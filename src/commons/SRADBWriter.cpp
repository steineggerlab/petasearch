#include "SRADBWriter.h"
#include "SRADBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "Concat.h"
#include "itoa.h"
#include "Timer.h"
#include "Parameters.h"


#define SIMDE_ENABLE_NATIVE_ALIASES

#include <simde/simde-common.h>

#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <unistd.h>

#ifdef OPENMP

#include <omp.h>

#endif

SRADBWriter::SRADBWriter(const char *dataFileName_, const char *indexFileName_, unsigned int threads, size_t mode, int
dbtype)
        : threads(threads), mode(mode), dbtype(dbtype) {
    dataFileName = strdup(dataFileName_);
    indexFileName = strdup(indexFileName_);

    dataFiles = new FILE *[threads];
    dataFilesBuffer = new char *[threads];

    dataFileNames = new char *[threads];

    indexFiles = new FILE *[threads];
    indexFileNames = new char *[threads];

    starts = new size_t[threads];
    std::fill(starts, starts + threads, 0);
    offsets = new size_t[threads];
    std::fill(offsets, offsets + threads, 0);
//    if ((mode & Parameters::WRITER_COMPRESSED_MODE) != 0) {
//        datafileMode = "wb+";
//    } else {
    datafileMode = "wb";
//    }

    closed = true;
}

SRADBWriter::~SRADBWriter() {
    delete[] offsets;
    delete[] starts;
    delete[] indexFileNames;
    delete[] indexFiles;
    delete[] dataFileNames;
    delete[] dataFilesBuffer;
    delete[] dataFiles;
    free(indexFileName);
    free(dataFileName);
}

// allocates heap memory, careful
char *SRADBWriter::makeResultFilename(const char *name, size_t split) {
    std::ostringstream ss;
    ss << name << "." << split;
    std::string s = ss.str();
    return strdup(s.c_str());
}

void SRADBWriter::open(size_t bufferSize) {
    if (bufferSize == SIZE_MAX) {
        if (Util::getTotalSystemMemory() < (8ull * 1024 * 1024 * 1024)) {
            // reduce this buffer if our system does not have much memory
            // createdb runs into trouble since it creates 2x32 splits with 64MB each (=4GB)
            // 8MB should be enough
            bufferSize = 8ull * 1024 * 1024;
        } else {
            bufferSize = 32ull * 1024 * 1024;
        }
    }
    for (unsigned int i = 0; i < threads; i++) {
        dataFileNames[i] = makeResultFilename(dataFileName, i);
        indexFileNames[i] = makeResultFilename(indexFileName, i);

        dataFiles[i] = FileUtil::openAndDelete(dataFileNames[i], datafileMode.c_str());
        int fd = fileno(dataFiles[i]);
        int flags;
        if ((flags = fcntl(fd, F_GETFL, 0)) < 0 || fcntl(fd, F_SETFD, flags | FD_CLOEXEC) == -1) {
            Debug(Debug::ERROR) << "Can not set mode for " << dataFileNames[i] << "!\n";
            EXIT(EXIT_FAILURE);
        }

        dataFilesBuffer[i] = new(std::nothrow) char[bufferSize];
        Util::checkAllocation(dataFilesBuffer[i], "Cannot allocate buffer for SRADBWriter");
        incrementMemory(bufferSize);
        this->bufferSize = bufferSize;

        // set buffer to 64
        if (setvbuf(dataFiles[i], dataFilesBuffer[i], _IOFBF, bufferSize) != 0) {
            Debug(Debug::WARNING) << "Write buffer could not be allocated (bufferSize=" << bufferSize << ")\n";
        }

        indexFiles[i] = FileUtil::openAndDelete(indexFileNames[i], "w");
        fd = fileno(indexFiles[i]);
        if ((flags = fcntl(fd, F_GETFL, 0)) < 0 || fcntl(fd, F_SETFD, flags | FD_CLOEXEC) == -1) {
            Debug(Debug::ERROR) << "Can not set mode for " << indexFileNames[i] << "!\n";
            EXIT(EXIT_FAILURE);
        }

        if (setvbuf(indexFiles[i], NULL, _IOFBF, bufferSize) != 0) {
            Debug(Debug::WARNING) << "Write buffer could not be allocated (bufferSize=" << bufferSize << ")\n";
        }

        if (dataFiles[i] == NULL) {
            perror(dataFileNames[i]);
            EXIT(EXIT_FAILURE);
        }

        if (indexFiles[i] == NULL) {
            perror(indexFileNames[i]);
            EXIT(EXIT_FAILURE);
        }
    }

    closed = false;
}

void SRADBWriter::writeDbtypeFile(const char *path, int dbtype, bool isCompressed) {
    if (dbtype == Parameters::DBTYPE_OMIT_FILE) {
        return;
    }

    std::string name = std::string(path) + ".dbtype";
    FILE *file = FileUtil::openAndDelete(name.c_str(), "wb");
    dbtype = isCompressed ? dbtype | (1 << 31) : dbtype & ~(1 << 31);
#if SIMDE_ENDIAN_ORDER == SIMDE_ENDIAN_BIG
    dbtype = __builtin_bswap32(dbtype);
#endif
    size_t written = fwrite(&dbtype, sizeof(int), 1, file);
    if (written != 1) {
        Debug(Debug::ERROR) << "Can not write to data file " << name << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(file) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << name << "\n";
        EXIT(EXIT_FAILURE);
    }
}


void SRADBWriter::close(bool merge) {
    // close all datafiles
    for (unsigned int i = 0; i < threads; i++) {
        if (fclose(dataFiles[i]) != 0) {
            Debug(Debug::ERROR) << "Cannot close data file " << dataFileNames[i] << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fclose(indexFiles[i]) != 0) {
            Debug(Debug::ERROR) << "Cannot close index file " << indexFileNames[i] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    merge = getenv("MMSEQS_FORCE_MERGE") != NULL || merge;
    mergeResults(dataFileName, indexFileName, (const char **) dataFileNames, (const char **) indexFileNames,
                 threads, merge);

    writeDbtypeFile(dataFileName, dbtype, (mode & Parameters::WRITER_COMPRESSED_MODE) != 0);

    for (unsigned int i = 0; i < threads; i++) {
        delete[] dataFilesBuffer[i];
        decrementMemory(bufferSize);
        free(dataFileNames[i]);
        free(indexFileNames[i]);
    }
    closed = true;
}

void SRADBWriter::writeStart(unsigned int thrIdx) {
    checkClosed();
    if (thrIdx >= threads) {
        Debug(Debug::ERROR) << "Thread index " << thrIdx << " > maximum thread number " << threads << "\n";
        EXIT(EXIT_FAILURE);
    }
    starts[thrIdx] = offsets[thrIdx];
}

size_t SRADBWriter::writeAdd(const char *data, size_t dataSize, unsigned int thrIdx) {
    checkClosed();
    if (thrIdx >= threads) {
        Debug(Debug::ERROR) << "Thread index " << thrIdx << " > maximum thread number " << threads << "\n";
        EXIT(EXIT_FAILURE);
    }
    size_t totalWriten = 0;
    size_t written;
    written = fwrite(data, sizeof(char), dataSize, dataFiles[thrIdx]);
    if (written != dataSize) {
        Debug(Debug::ERROR) << "Can not write to data file " << dataFileNames[thrIdx] << "\n";
        EXIT(EXIT_FAILURE);
    }
    offsets[thrIdx] += written;

    return totalWriten;
}

void SRADBWriter::writeEnd(unsigned int thrIdx, bool addNullByte, bool addIndexEntry) {
    size_t totalWritten = 0;
// entries are always separated by a null byte
    if (addNullByte == true) {
        char nullByte = '\0';
        const size_t written = fwrite(&nullByte, sizeof(char), 1, dataFiles[thrIdx]);
        if (written != 1) {
            Debug(Debug::ERROR) << "Can not write to data file " << dataFileNames[thrIdx] << "\n";
            EXIT(EXIT_FAILURE);
        }
        totalWritten += written;
        offsets[thrIdx] += 1;
    }

    if (addIndexEntry == true) {
//        size_t length = offsets[thrIdx] - starts[thrIdx];
// keep original size in index
        writeIndexEntry(starts[thrIdx], thrIdx);
    }
}

void SRADBWriter::writeIndexEntry(size_t offset, unsigned int thrIdx) {
    char buffer[1024];
    size_t len = indexToBuffer(buffer,offset);
    size_t written = fwrite(buffer, sizeof(char), len, indexFiles[thrIdx]);
    if (written != len) {
        Debug(Debug::ERROR) << "Can not write to data file " << dataFileName[thrIdx] << "\n";
        EXIT(EXIT_FAILURE);
    }
}


void SRADBWriter::writeData(const char *data, size_t dataSize, unsigned int key, unsigned int thrIdx, bool addNullByte,
                            bool addIndexEntry) {
    writeStart(thrIdx);
    writeAdd(data, dataSize, thrIdx);
    writeEnd(thrIdx, addNullByte, addIndexEntry);
}

size_t SRADBWriter::indexToBuffer(char *buff1, size_t offsetStart) {
    char *basePos = buff1;
    char *tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(offsetStart), buff1); //Itoa::u32toa_sse2
    // (static_cast<uint32_t>(key), buff1);
//    *(tmpBuff - 1) = '\t';
    // tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(offsetStart), tmpBuff);
//    *(tmpBuff - 1) = '\t';
//    tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(len), tmpBuff);
    *(tmpBuff - 1) = '\n';
    *(tmpBuff) = '\0';
    return tmpBuff - basePos;
}

void SRADBWriter::checkClosed() {
    if (closed == true) {
        Debug(Debug::ERROR) << "Trying to read a closed database. Datafile=" << dataFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void SRADBWriter::mergeResults(const std::string &outFileName, const std::string &outFileNameIndex,
                               const std::vector<std::pair<std::string, std::string >> &files) {
    const char **datafilesNames = new const char *[files.size()];
    const char **indexFilesNames = new const char *[files.size()];
    for (size_t i = 0; i < files.size(); i++) {
        datafilesNames[i] = files[i].first.c_str();
        indexFilesNames[i] = files[i].second.c_str();
    }
    mergeResults(outFileName.c_str(), outFileNameIndex.c_str(), datafilesNames, indexFilesNames, files.size(), true);
    delete[] datafilesNames;
    delete[] indexFilesNames;

    // leave only one dbtype file behind
    if (files.size() > 0) {
        std::string typeSrc = files[0].first + ".dbtype";
        std::string typeDest = outFileName + ".dbtype";
        if (FileUtil::fileExists(typeSrc.c_str())) {
            std::rename(typeSrc.c_str(), typeDest.c_str());
        }
        for (size_t i = 1; i < files.size(); i++) {
            std::string typeFile = files[i].first + ".dbtype";
            if (FileUtil::fileExists(typeFile.c_str())) {
                FileUtil::remove(typeFile.c_str());
            }
        }
    }
}
//
//template<>
//void SRADBWriter::writeIndexEntryToFile(FILE *outFile, char *buff1, SRADBReader<unsigned int>::Index &index) {
//    char *tmpBuff = Itoa::u32toa_sse2((uint32_t) index.id, buff1);
//    *(tmpBuff - 1) = '\t';
//    size_t currOffset = index.offset;
//    tmpBuff = Itoa::u64toa_sse2(currOffset, tmpBuff);
//    *(tmpBuff - 1) = '\t';
//    uint32_t sLen = index.length;
//    tmpBuff = Itoa::u32toa_sse2(sLen, tmpBuff);
//    *(tmpBuff - 1) = '\n';
//    *(tmpBuff) = '\0';
//    fwrite(buff1, sizeof(char), (tmpBuff - buff1), outFile);
//}
//
//template<>
//void SRADBWriter::writeIndexEntryToFile(FILE *outFile, char *buff1, SRADBReader<std::string>::Index &index) {
//    size_t keyLen = index.id.length();
//    char *tmpBuff = (char *) memcpy((void *) buff1, (void *) index.id.c_str(), keyLen);
//    tmpBuff += keyLen;
//    *(tmpBuff) = '\t';
//    tmpBuff++;
//    size_t currOffset = index.offset;
//    tmpBuff = Itoa::u64toa_sse2(currOffset, tmpBuff);
//    *(tmpBuff - 1) = '\t';
//    uint32_t sLen = index.length;
//    tmpBuff = Itoa::u32toa_sse2(sLen, tmpBuff);
//    *(tmpBuff - 1) = '\n';
//    *(tmpBuff) = '\0';
//    fwrite(buff1, sizeof(char), (tmpBuff - buff1), outFile);
//}
//
//template<>
//void SRADBWriter::writeIndex(FILE *outFile, size_t indexSize, SRADBReader<unsigned int>::Index *index) {
//    char buff1[1024];
//    for (size_t id = 0; id < indexSize; id++) {
//        writeIndexEntryToFile(outFile, buff1, index[id]);
//    }
//}
//
//template<>
//void SRADBWriter::writeIndex(FILE *outFile, size_t indexSize, SRADBReader<std::string>::Index *index) {
//    char buff1[1024];
//    for (size_t id = 0; id < indexSize; id++) {
//        writeIndexEntryToFile(outFile, buff1, index[id]);
//    }
//}
//

void SRADBWriter::mergeResults(const char *outFileName, const char *outFileNameIndex,
                               const char **dataFileNames, const char **indexFileNames,
                               unsigned long fileCount, bool mergeDatafiles) {
    Timer timer;
    std::vector<std::vector<std::string>> dataFilenames;
    for (unsigned int i = 0; i < fileCount; ++i) {
        dataFilenames.emplace_back(FileUtil::findDatafiles(dataFileNames[i]));
    }

    // merge results into one result file
    if (dataFilenames.size() > 1) {
        std::vector<FILE *> datafiles;
        std::vector<size_t> mergedSizes;
        for (unsigned int i = 0; i < dataFilenames.size(); i++) {
            std::vector<std::string> &filenames = dataFilenames[i];
            size_t cumulativeSize = 0;
            for (size_t j = 0; j < filenames.size(); ++j) {
                FILE *fh = fopen(filenames[j].c_str(), "r");
                if (fh == NULL) {
                    Debug(Debug::ERROR) << "Can not open result file " << filenames[j] << "!\n";
                    EXIT(EXIT_FAILURE);
                }
                struct stat sb;
                if (fstat(fileno(fh), &sb) < 0) {
                    int errsv = errno;
                    Debug(Debug::ERROR) << "Failed to fstat file " << filenames[j] << ". Error " << errsv << ".\n";
                    EXIT(EXIT_FAILURE);
                }
                datafiles.emplace_back(fh);
                cumulativeSize += sb.st_size;
            }
            mergedSizes.push_back(cumulativeSize);
        }

        if (mergeDatafiles) {
            FILE *outFh = FileUtil::openAndDelete(outFileName, "w");
            Concat::concatFiles(datafiles, outFh);
            if (fclose(outFh) != 0) {
                Debug(Debug::ERROR) << "Cannot close data file " << outFileName << "\n";
                EXIT(EXIT_FAILURE);
            }
        }

        for (unsigned int i = 0; i < datafiles.size(); ++i) {
            if (fclose(datafiles[i]) != 0) {
                Debug(Debug::ERROR) << "Cannot close data file in merge\n";
                EXIT(EXIT_FAILURE);
            }
        }

        if (mergeDatafiles) {
            for (unsigned int i = 0; i < dataFilenames.size(); i++) {
                std::vector<std::string> &filenames = dataFilenames[i];
                for (size_t j = 0; j < filenames.size(); ++j) {
                    FileUtil::remove(filenames[j].c_str());
                }
            }
        }

        // merge index
        mergeIndex(indexFileNames, dataFilenames.size(), mergedSizes);
    } else if (dataFilenames.size() == 1) {
        std::vector<std::string> &filenames = dataFilenames[0];
        if (filenames.size() == 1) {
            // In single thread SRADBReader mode it will create a .0
            // that should be moved to the final destination dest instead of dest.0
            FileUtil::move(filenames[0].c_str(), outFileName);
        } else {
            DBReader<unsigned int>::moveDatafiles(filenames, outFileName);
        }
    } else {
        FILE *outFh = FileUtil::openAndDelete(outFileName, "w");
        if (fclose(outFh) != 0) {
            Debug(Debug::ERROR) << "Cannot close data file " << outFileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        outFh = FileUtil::openAndDelete(outFileNameIndex, "w");
        if (fclose(outFh) != 0) {
            Debug(Debug::ERROR) << "Cannot close index file " << outFileNameIndex << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    if (dataFilenames.size() > 0) {
        FileUtil::move(indexFileNames[0], outFileNameIndex);

    }
    Debug(Debug::INFO) << "Time for merging to " << FileUtil::baseName(outFileName) << ": " << timer.lap() << "\n";
}

void SRADBWriter::mergeIndex(const char **indexFilenames, unsigned int fileCount, const std::vector<size_t>
&dataSizes) {
    FILE *index_file = fopen(indexFilenames[0], "a");
    if (index_file == NULL) {
        perror(indexFilenames[0]);
        EXIT(EXIT_FAILURE);
    }
    size_t globalOffset = dataSizes[0];
    for (unsigned int fileIdx = 1; fileIdx < fileCount; fileIdx++) {
        FILE *indexFD = FileUtil::openAndDelete(indexFilenames[fileIdx], "r");
        size_t size = FileUtil::getFileSize(std::string(indexFilenames[fileIdx]));
        char *string = static_cast<char *>(malloc(size + 1));
        char *l = static_cast<char *>(malloc(1024));
        fread(string, size, 1, indexFD);
        while (Util::getLine(string, size + 1, l, 1024)){
            if (Util::isNumber(std::string(l))) {
                std::string toWrite = std::to_string(
                        globalOffset + Util::fast_atoi<unsigned long>(l));
                fwrite(toWrite.c_str(), toWrite.size(), 1, index_file);
            }
        }


//        SRADBReader<unsigned int> reader(indexFilenames[fileIdx], indexFilenames[fileIdx], 1,
//                                      SRADBReader<unsigned int>::USE_INDEX);
//        reader.open(SRADBReader<unsigned int>::HARDNOSORT);
//        if (reader.getSize() > 0) {
//            SRADBReader<unsigned int>::Index *index = reader.getIndex();
//            for (size_t i = 0; i < reader.getSize(); i++) {
//                size_t currOffset = index[i].offset;
//                index[i].offset = globalOffset + currOffset;
//            }
//            writeIndex(index_file, reader.getSize(), index);
//        }
//        reader.close();
//        FileUtil::remove(indexFilenames[fileIdx]);

        globalOffset += dataSizes[fileIdx];
    }
    if (fclose(index_file) != 0) {
        Debug(Debug::ERROR) << "Cannot close index file " << indexFilenames[0] << "\n";
        EXIT(EXIT_FAILURE);
    }
}

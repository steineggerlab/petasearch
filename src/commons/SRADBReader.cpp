#include "SRADBReader.h"
#include "DBReader.h"
#include "BitManipulateMacros.h"
#include "FastSort.h"
#include <algorithm>
#include <cstring>
#include <cstddef>

#include <sys/mman.h>
#include <sys/stat.h>

#include <fcntl.h>

#include "MemoryMapped.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"

#ifdef OPENMP

#include <omp.h>

#endif


SRADBReader::SRADBReader(const char *dataFileName, const char *indexFileName, int threads, int mode) :
        threads(threads), dataMode(mode), dataFileName(strdup(dataFileName)), indexFileName(strdup(indexFileName)) {
}

void SRADBReader::open(int accessType) {
    this->accessType = accessType;
    if (dataFileName != NULL) {
        dbType = FileUtil::parseDbType(dataFileName);
    }
    if (dataMode & DBReader<unsigned int>::USE_DATA) {
        dataFileNames = FileUtil::findDatafiles(dataFileName);
        if (dataFileNames.empty()) {
            Debug(Debug::ERROR) << "No datafile could be found for " << dataFileName << "!\n";
            EXIT(EXIT_FAILURE);
        }
        totalDataSize = 0;
        dataFileCnt = dataFileNames.size();
        dataSizeOffset = new size_t[dataFileNames.size() + 1];
        dataFiles = new char *[dataFileNames.size()];
        for (size_t fileIdx = 0; fileIdx < dataFileNames.size(); fileIdx++) {
            FILE *dataFile = fopen(dataFileNames[fileIdx].c_str(), "r");
            if (dataFile == NULL) {
                Debug(Debug::ERROR) << "Cannot open data file " << dataFileName << "!\n";
                EXIT(EXIT_FAILURE);
            }
            size_t dataSize;
            dataFiles[fileIdx] = mmapData(dataFile, &dataSize);
            dataSizeOffset[fileIdx] = totalDataSize;
            totalDataSize += dataSize;
            if (fclose(dataFile) != 0) {
                Debug(Debug::ERROR) << "Cannot close file " << dataFileName << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        dataSizeOffset[dataFileNames.size()] = totalDataSize;
        dataMapped = true;
        if (accessType == DBReader<unsigned int>::LINEAR_ACCCESS ||
            accessType == DBReader<unsigned int>::SORT_BY_OFFSET) {
            setSequentialAdvice();
        }
        if (dataMode & DBReader<unsigned int>::USE_LOOKUP || dataMode & DBReader<unsigned int>::USE_LOOKUP_REV) {
            Debug(Debug::ERROR) << "SRAsearch Database does not have lookup file " << "!\n";
            EXIT(EXIT_FAILURE);
        }

        MemoryMapped indexData(indexFileName, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
        if (!indexData.isValid()) {
            Debug(Debug::ERROR) << "Cannot open index file " << indexFileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        char *indexDataChar = (char *) indexData.getData();
        size_t indexDataSize = indexData.size();
        size = Util::ompCountLines(indexDataChar, indexDataSize, threads);

        index = new unsigned long[size];
        Util::checkAllocation(index, "Cannot allocate index memory in DBReader");
        incrementMemory(sizeof(unsigned long) * size);

        readIndex(indexDataChar, indexDataSize, index);
//        bool isSortedById = readIndex(indexDataChar, indexDataSize, index, dataSize);
        indexData.close();
        // sortIndex also handles access modes that don't require sorting
    }

    closed = 0;
}

void SRADBReader::setSequentialAdvice() {
#ifdef HAVE_POSIX_MADVISE
    for (size_t i = 0; i < dataFileCnt; i++) {
        size_t dataSize = dataSizeOffset[i + 1] - dataSizeOffset[i];
        if (dataSize > 0 && posix_madvise(dataFiles[i], dataSize, POSIX_MADV_SEQUENTIAL) != 0) {
            Debug(Debug::ERROR) << "posix_madvise returned an error " << dataFileName << "\n";
        }
    }
#endif
}

void SRADBReader::readIndex(char *data, size_t indexDataSize, unsigned long *index) {
#ifdef OPENMP
    int threadCnt = 1;
    const int totalThreadCnt = threads;
    if (totalThreadCnt >= 4) {
        threadCnt = 4;
    }
#endif
    size_t globalIdOffset = 0;

    unsigned int localMaxSeqLen = 0;
    size_t localDataSize = 0;
    const unsigned int BATCH_SIZE = 1048576;
#pragma omp parallel num_threads(threadCnt) reduction(max: localMaxSeqLen) reduction(+: localDataSize)
    {
        size_t currPos = 0;
        char *indexDataChar = (char *) data;
        const char *cols[3];
        size_t lineStartId = __sync_fetch_and_add(&(globalIdOffset), BATCH_SIZE);
        size_t currLine = 0;
        size_t prev_offset = 0;

        while (currPos < indexDataSize) {
            if (currLine >= this->size) {
                Debug(Debug::ERROR) << "Corrupt memory, too many entries: " << currLine << " >= " << this->size << "\n";
                EXIT(EXIT_FAILURE);
            }
            if (currLine == lineStartId) {
                for (size_t startIndex = lineStartId;
                     startIndex < lineStartId + BATCH_SIZE && currPos < indexDataSize; startIndex++) {
                    Util::getWordsOfLine(indexDataChar, cols, 3);
                    size_t offset = Util::fast_atoi<size_t>(cols[0]);
                    size_t length = offset - prev_offset + 1;
                    localDataSize += length;
                    index[startIndex] = offset;
//                    index[startIndex].length = length;
                    localMaxSeqLen = std::max(static_cast<unsigned int>(length), localMaxSeqLen);
                    indexDataChar = Util::skipLine(indexDataChar);
                    currPos = indexDataChar - (char *) data;
                    currLine++;
                }
                lineStartId = __sync_fetch_and_add(&(globalIdOffset), BATCH_SIZE);
            } else {
                indexDataChar = Util::skipLine(indexDataChar);
                currPos = indexDataChar - (char *) data;
                currLine++;
            }

        }
    }
//    dataSize = localDataSize;
    maxSeqLen = localMaxSeqLen;
    seqBuffer = (char *) calloc(maxSeqLen * 3 / 2 + 1, sizeof(char));
}

void SRADBReader::close() {
    if (dataMode & DBReader<unsigned int>::USE_DATA) {
        unmapData();
    }
    closed = 1;
}

int SRADBReader::getDbtype() {
    return dbType;
}

char *SRADBReader::mmapData(FILE *file, size_t *dataSize) {
    struct stat sb;
    if (fstat(fileno(file), &sb) < 0) {
        int errsv = errno;
        Debug(Debug::ERROR) << "Failed to fstat File=" << dataFileName << ". Error " << errsv << ".\n";
        EXIT(EXIT_FAILURE);
    }

    *dataSize = sb.st_size;
    int fd = fileno(file);

    char *ret;
    if (*dataSize > 0) {
        if ((dataMode & DBReader<unsigned int>::USE_FREAD) == 0) {
            int mode;
            if (dataMode & DBReader<unsigned int>::USE_WRITABLE) {
                mode = PROT_READ | PROT_WRITE;
            } else {
                mode = PROT_READ;
            }
            ret = static_cast<char *>(mmap(NULL, *dataSize, mode, MAP_PRIVATE, fd, 0));
            if (ret == MAP_FAILED) {
                int errsv = errno;
                Debug(Debug::ERROR) << "Failed to mmap memory dataSize=" << *dataSize << " File=" << dataFileName
                                    << ". Error " << errsv << ".\n";
                EXIT(EXIT_FAILURE);
            }
        } else {
            ret = static_cast<char *>(malloc(*dataSize));
            Util::checkAllocation(ret, "Not enough system memory to read in the whole data file.");
            incrementMemory(*dataSize);

            size_t result = fread(ret, 1, *dataSize, file);
            if (result != *dataSize) {
                Debug(Debug::ERROR) << "Failed to read in datafile (" << dataFileName << "). Error " << errno << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        return ret;
    } else {
        return NULL;
    }
}

void SRADBReader::unmapData() {
    if (dataMapped == true) {
        for (size_t fileIdx = 0; fileIdx < dataFileNames.size(); fileIdx++) {
            size_t fileSize = dataSizeOffset[fileIdx + 1] - dataSizeOffset[fileIdx];
            if (fileSize > 0) {
                if ((dataMode & DBReader<unsigned int>::USE_FREAD) == 0) {
                    if (munmap(dataFiles[fileIdx], fileSize) < 0) {
                        Debug(Debug::ERROR) << "Failed to munmap memory dataSize=" << fileSize << " File="
                                            << dataFileName
                                            << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                } else {
                    free(dataFiles[fileIdx]);
                }
            }
        }
    }
    dataMapped = false;
}

char *SRADBReader::getData(size_t id, int thread_idx) {
    char *rawString = getDataUncompressed(id);
    unsigned short *packedArray = reinterpret_cast<unsigned short *>(rawString);
    // FIXME: seqLen calculation is WRONG!
    int idex = 0;
    int j = 0;
    while (!IS_LAST_15_BITS(packedArray[idex])) {
        seqBuffer[j] = GET_HIGH_CHAR(packedArray[idex]);
        seqBuffer[j + 1] = GET_MID_CHAR(packedArray[idex]);
        seqBuffer[j + 2] = GET_LOW_CHAR(packedArray[idex]);
        j += 3;
        idex++;
    }
    char p1 = GET_HIGH_CHAR(packedArray[idex]);
    char p2 = GET_MID_CHAR(packedArray[idex]);
    char p3 = GET_LOW_CHAR(packedArray[idex]);
    seqBuffer[j] = p1;
    seqBuffer[j + 1] = '\0';
    if (p2 != '@') {
        seqBuffer[j + 1] = p2;
        seqBuffer[j + 2] = '\0';
    }
    if (p3 != '@') {
        seqBuffer[j + 2] = p3;
        seqBuffer[j + 3] = '\0';
    }
    return seqBuffer;
}

char *SRADBReader::getDataUncompressed(size_t id) {
    checkClosed();
    if (!(dataMode & DBReader<unsigned int>::USE_DATA)) {
        Debug(Debug::ERROR) << "DBReader is just open in INDEXONLY mode. Call of getData is not allowed" << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (id >= size) {
        Debug(Debug::ERROR) << "Invalid database read for database data file=" << dataFileName << ", database index="
                            << indexFileName << "\n";
        Debug(Debug::ERROR) << "getData: local id (" << id << ") >= db size (" << size << ")\n";
        EXIT(EXIT_FAILURE);
    }
    return getDataByOffset(index[id]);
}

char *SRADBReader::getDataByOffset(size_t offset) {
    if (offset >= totalDataSize) {
        Debug(Debug::ERROR) << "Invalid database read for database data file=" << dataFileName << ", database index="
                            << indexFileName << "\n";
        Debug(Debug::ERROR) << "Size of data: " << totalDataSize << "\n";
        Debug(Debug::ERROR) << "Requested offset: " << offset << "\n";
        EXIT(EXIT_FAILURE);
    }
    size_t cnt = 0;
    while (!(offset >= dataSizeOffset[cnt] && offset < dataSizeOffset[cnt + 1])) {
        cnt++;
    }
    size_t fileOffset = offset - dataSizeOffset[cnt];
    return dataFiles[cnt] + fileOffset;
}

size_t SRADBReader::getSize() {
    checkClosed();
    return size;
}

size_t SRADBReader::getSeqLen(size_t id) {
//    return id;
    if (id < size - 1) {
        return (index[id + 1] - index[id]) / 2 * 3;
    } else if (id == size - 1) { // this is the last element
        return (totalDataSize - index[id]) / 2 * 3;
    } else { // invalid id
        Debug(Debug::ERROR) << "Invalid database read for database data file=" << dataFileName << ", database index="
                            << indexFileName << "\n";
        Debug(Debug::ERROR) << "local id (" << id << ") >= db size (" << size << ")\n";
        EXIT(EXIT_FAILURE);
    }
}

unsigned int SRADBReader::getDbKey(size_t id) {
    return id;
}

size_t SRADBReader::getAminoAcidDBSize() {
    checkClosed();
//    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE) || Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_PROFILE_STATE_PROFILE)) {
//        // Get the actual profile column without the null byte per entry
//        return (dataSize / Sequence::PROFILE_READIN_SIZE) - size;
//    } else {
    // Get the actual number of residues witout \n\0 per entry
    return totalDataSize * 3 / 2 - (2 * size);
//    }
//    return _dbreader->getAminoAcidDBSize() * 3 / 2;
}

void SRADBReader::checkClosed() {
    if (closed == 1) {
        Debug(Debug::ERROR) << "Trying to read a closed database.\n";
        EXIT(EXIT_FAILURE);
    }
}

SRADBReader::~SRADBReader() {
    delete[] index;
    free(seqBuffer);

    if (dataFileName != NULL) {
        free(dataFileName);
    }

    if (indexFileName != NULL) {
        free(indexFileName);
    }

    if (dataSizeOffset != NULL) {
        delete[] dataSizeOffset;
    }

    if (dataFiles != NULL) {
        delete[] dataFiles;
    }
}

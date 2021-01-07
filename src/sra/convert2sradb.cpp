//
// Created by matchy233
//

#include "FileUtil.h"
#include "DBReader.h"
#include "BitManipulateMacros.h"

#include "LocalParameters.h"
#include "DBWriter.h"
#include "KSeqWrapper.h"

int convert2sradb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    char newline = '\n';

    std::vector<std::string> filenames(par.filenames);
    std::string dataFile = filenames.back();
    filenames.pop_back();

    for (size_t i = 0; i < filenames.size(); i++) {
        if (FileUtil::directoryExists(filenames[i].c_str())) {
            Debug(Debug::ERROR) << "File " << filenames[i] << " is a directory\n";
            EXIT(EXIT_FAILURE);
        }
    }

    bool dbInput = false;
    if (FileUtil::fileExists(par.db1dbtype.c_str())) {
        if (filenames.size() > 1) {
            Debug(Debug::ERROR) << "Only one database can be used with database input\n";
            EXIT(EXIT_FAILURE);
        }
        dbInput = true;
    }

    int dbType = LocalParameters::DBTYPE_SRA_DB;

    std::string indexFile = dataFile + ".index";
    const unsigned int shuffleSplits = par.shuffleDatabase ? 32 : 1;
    std::string hdrDataFile = dataFile + "_h";
    std::string hdrIndexFile = dataFile + "_h.index";

    unsigned int entries_num = 0;

    Debug::Progress progress;
    std::vector<unsigned short> *sourceLookup = new std::vector<unsigned short>[shuffleSplits]();
    for (size_t i = 0; i < shuffleSplits; ++i) {
        sourceLookup[i].reserve(16384);
    }

    Debug(Debug::INFO) << "Converting sequences\n";

    std::string sourceFile = dataFile + ".source";

    FILE *source = fopen(sourceFile.c_str(), "w");
    if (source == NULL) {
        Debug(Debug::ERROR) << "Cannot open " << sourceFile << " for writing\n";
        EXIT(EXIT_FAILURE);
    }

    DBWriter hdrWriter(hdrDataFile.c_str(), hdrIndexFile.c_str(), shuffleSplits, par.compressed,
                       Parameters::DBTYPE_GENERIC_DB);
    hdrWriter.open();
    DBWriter seqWriter(dataFile.c_str(), indexFile.c_str(), shuffleSplits, par.compressed, dbType);
    seqWriter.open();

    size_t fileCount = filenames.size();
    DBReader<unsigned int> *reader = NULL;
    if (dbInput) {
        reader = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(), 1,
                                            DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX |
                                            DBReader<unsigned int>::USE_LOOKUP);
        reader->open(DBReader<unsigned int>::LINEAR_ACCCESS);
        fileCount = reader->getSize();
    }

    for (size_t fileIdx = 0; fileIdx < fileCount; fileIdx++) {
        unsigned int numEntriesInCurrFile = 0;
        std::string header;
        header.reserve(1024);

        std::string sourceName;
        if (dbInput) {
            sourceName = reader->getLookupEntryName(fileIdx);
        } else {
            sourceName = FileUtil::baseName(filenames[fileIdx]);
        }
        char buffer[4096];
        size_t len = snprintf(buffer, sizeof(buffer), "%zu\t%s\n", fileIdx, sourceName.c_str());
        int written = fwrite(buffer, sizeof(char), len, source);
        if (written != (int) len) {
            Debug(Debug::ERROR) << "Cannot write to source file " << sourceFile << "\n";
            EXIT(EXIT_FAILURE);
        }

        KSeqWrapper *kseq = NULL;
        if (dbInput) {
            kseq = new KSeqBuffer(reader->getData(fileIdx, 0), reader->getEntryLen(fileIdx) - 1);
        } else {
            kseq = KSeqFactory(filenames[fileIdx].c_str());
        }

        while (kseq->ReadEntry()) {
            progress.updateProgress();
            const KSeqWrapper::KSeqEntry &e = kseq->entry;

            if (e.name.l == 0) {
                Debug(Debug::ERROR) << "Fasta entry " << entries_num << " is invalid\n";
                EXIT(EXIT_FAILURE);
            }

            // Header craetion
            header.append(e.name.s, e.name.l);
            if (e.comment.l > 0) {
                header.append(" ", 1);
                header.append(e.comment.s, e.comment.l);
            }

            std::string headerId = Util::parseFastaHeader(header.c_str());
            if (headerId.empty()) {
                // An identifier is necessary for these two cases, so we should just give up
                Debug(Debug::WARNING) << "Cannot extract identifier from entry " << entries_num << "\n";
            }
            header.push_back('\n');

            unsigned int id = entries_num;

            // Writie down entry in a compact way!
            unsigned int splitIdx = id % shuffleSplits;
            sourceLookup[splitIdx].emplace_back(fileIdx);

            hdrWriter.writeData(header.c_str(), header.length(), id, splitIdx);

            int rem = e.sequence.l % 3;
            int padding = rem == 0 ? 0 : 1;
            unsigned short result[e.sequence.l / 3 + padding];

            int i = 0;
            for (; i < e.sequence.l - 3; i = i + 3) {
                result[i / 3] = PACK_TO_SHORT(e.sequence.s[i], e.sequence.s[i + 1], e.sequence.s[i + 2]);
            }
            result[i] = 0;
            for (int j = 0; j < rem; j++) {
                result[i] <<= 5U;
                result[i] |= GET_LAST_5_BITS(e.sequence.s[i + j]);
            }
            result[i] |= 0x1000U; // Set last bit
            const char *packedSeq = reinterpret_cast<const char *>(result);

            seqWriter.writeStart(splitIdx);
            seqWriter.writeAdd(packedSeq, sizeof(packedSeq), splitIdx);
            seqWriter.writeAdd(&newline, 1, splitIdx);
            seqWriter.writeEnd(id, splitIdx, true);

            entries_num++;
            numEntriesInCurrFile++;
            header.clear();
        }
        delete kseq;
    }
    Debug(Debug::INFO) << "\n";
    if (fclose(source) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << sourceFile << "\n";
        EXIT(EXIT_FAILURE);
    }
    hdrWriter.close(true, false);
    seqWriter.close(true, false);
    Debug(Debug::INFO) << "Database type: SRA Database" << "\n";

    if (dbInput) {
        reader->close();
        delete reader;
    }

    if (entries_num == 0) {
        Debug(Debug::ERROR) << "The input files have no entry: ";
        for (size_t fileIdx = 0; fileIdx < filenames.size(); fileIdx++) {
            Debug(Debug::ERROR) << " - " << filenames[fileIdx] << "\n";
        }
        Debug(Debug::ERROR) << "Please check your input files. Only files in fasta/fastq[.gz|bz2] are supported\n";
        EXIT(EXIT_FAILURE);
    }

    if (par.writeLookup == true) {
        DBReader<unsigned int> readerHeader(hdrDataFile.c_str(), hdrIndexFile.c_str(), 1, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        readerHeader.open(DBReader<unsigned int>::NOSORT);
        // create lookup file
        std::string lookupFile = dataFile + ".lookup";
        FILE* file = FileUtil::openAndDelete(lookupFile.c_str(), "w");
        std::string buffer;
        buffer.reserve(2048);
        unsigned int splitIdx = 0;
        unsigned int splitCounter = 0;
        DBReader<unsigned int>::LookupEntry entry;
        for (unsigned int id = 0; id < readerHeader.getSize(); id++) {
            size_t splitSize = sourceLookup[splitIdx].size();
            if (splitSize == 0 || splitCounter > sourceLookup[splitIdx].size() - 1) {
                splitIdx++;
                splitCounter = 0;
            }
            char *header = readerHeader.getData(id, 0);
            entry.id = id;
            entry.entryName = Util::parseFastaHeader(header);
            if (entry.entryName.empty()) {
                Debug(Debug::WARNING) << "Cannot extract identifier from entry " << entries_num << "\n";
            }
            entry.fileNumber = sourceLookup[splitIdx][splitCounter];
            readerHeader.lookupEntryToBuffer(buffer, entry);
            int written = fwrite(buffer.c_str(), sizeof(char), buffer.size(), file);
            if (written != (int)buffer.size()) {
                Debug(Debug::ERROR) << "Cannot write to lookup file " << lookupFile << "\n";
                EXIT(EXIT_FAILURE);
            }
            buffer.clear();
            splitCounter++;
        }
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << lookupFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        readerHeader.close();
    }
    delete[] sourceLookup;
}

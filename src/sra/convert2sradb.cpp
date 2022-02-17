
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "BitManipulateMacros.h"

#include "LocalParameters.h"
#include "SRADBWriter.h"
#include "KSeqWrapper.h"

#include "SRAUtil.h"

int convert2sradb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    char newline = '\n';

    /* Retrieve all input files*/
    std::vector<std::string> filenames(par.filenames);

    /* Get the last input file as the output file */
    std::string outputDataFile = filenames.back();
    filenames.pop_back();

    /* Determine whether any input file is a directory */
    for (size_t i = 0; i < filenames.size(); i++) {
        if (FileUtil::directoryExists(filenames[i].c_str())) {
            Debug(Debug::ERROR) << "File " << filenames[i] << " is a directory" << newline;
            EXIT(EXIT_FAILURE);
        }
    }

    /* Determine whether it is fasta input or database input
     *    Here we assume that the input is a database if it has a
     *    corresponding ".dbtype" file
     */
    bool isDbInput = FileUtil::fileExists(par.db1dbtype.c_str());

    /* If it is database input, we only accept 1 db each time*/
    if (isDbInput && (filenames.size() > 1)) {
        Debug(Debug::ERROR) << "Only one database can be used with database input" << newline;
        EXIT(EXIT_FAILURE);
    }


    /* Name output files and database type */
    int outputDbType = Parameters::DBTYPE_AMINO_ACIDS;

    const unsigned int shuffleSplits = 1;

    std::string outputIndexFile = outputDataFile + ".index";
    std::string outputHdrDataFile = outputDataFile + "_h";
    std::string outputHdrIndexFile = outputDataFile + "_h.index";
    std::string sourceFile = outputDataFile + ".source";

    unsigned int entries_num = 0;

    Debug::Progress progress;

    Debug(Debug::INFO) << "Converting sequences" << newline;

    FILE *source = fopen(sourceFile.c_str(), "w");
    if (source == NULL) {
        Debug(Debug::ERROR) << "Cannot open " << sourceFile << " for writing\n";
        EXIT(EXIT_FAILURE);
    }

    SRADBWriter hdrWriter(outputHdrDataFile.c_str(), outputHdrIndexFile.c_str(),
                          shuffleSplits, par.compressed,
                          Parameters::DBTYPE_GENERIC_DB);
    hdrWriter.open();
    SRADBWriter seqWriter(outputDataFile.c_str(), outputIndexFile.c_str(),
                          shuffleSplits, par.compressed,
                          Parameters::DBTYPE_AMINO_ACIDS);
    seqWriter.open();

    size_t fileCount = filenames.size();
    DBReader<unsigned int> *reader = NULL;
    DBReader<unsigned int> *hdrReader = NULL;
    if (isDbInput) {
        reader = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(), 1,
                                            DBReader<unsigned int>::USE_DATA |
                                            DBReader<unsigned int>::USE_INDEX |
                                            DBReader<unsigned int>::USE_LOOKUP);
        reader->open(DBReader<unsigned int>::NOSORT);
        hdrReader = new DBReader<unsigned int>((par.db1 + "_h").c_str(), (par.db1 + "_h.index").c_str(), 1,
                                               DBReader<unsigned int>::USE_DATA |
                                               DBReader<unsigned int>::USE_INDEX);
        hdrReader->open(DBReader<unsigned int>::NOSORT);
        fileCount = reader->getSize();
    }

    Debug(Debug::INFO) << "max sequence length: " << par.maxSeqLen << newline;

    /* Process all inputs */
    for (size_t fileIdx = 0; fileIdx < fileCount; fileIdx++) {
        unsigned int numEntriesInCurrFile = 0;
        std::string header;
        header.reserve(1024);

        std::string sourceName;
        if (isDbInput) {
            unsigned int dbKey = reader->getDbKey(fileIdx);
            size_t lookupId = reader->getLookupIdByKey(dbKey);
            sourceName = reader->getLookupEntryName(lookupId);
        } else {
            sourceName = FileUtil::baseName(filenames[fileIdx]);
        }
        char buffer[4096];
        size_t len = snprintf(buffer, sizeof(buffer), "%zu\t%s\n", fileIdx, sourceName.c_str());
        int written = fwrite(buffer, sizeof(char), len, source);
        if (written != (int) len) {
            Debug(Debug::ERROR) << "Cannot write to source file " << sourceFile << newline;
            EXIT(EXIT_FAILURE);
        }

        // TODO: fix kseqbuffer so that we can use the same API
//        KSeqWrapper *kseq = NULL;
//        if (isDbInput) {
//            std::string fakeString = " >sp|Q8I6R7|ACN2_ACAGO Acanthoscurrin-2 (Fragment) OS=Acanthoscurria gomesiana "
//                                     "OX=115339 GN=acantho2 PE=1 SV=1\n" + std::string(reader->getData(fileIdx, 0));
//            kseq = new KSeqBuffer(fakeString.c_str(), reader->getEntryLen(fileIdx));
//        } else {
//            kseq = KSeqFactory(filenames[fileIdx].c_str());
//        }

        if (isDbInput) {
            // DB INPUT CASE
            progress.updateProgress();
            const size_t s_upperlimit = reader->getSeqLen(fileIdx);
            char *kseq = static_cast<char *>(calloc(s_upperlimit + 1, sizeof(char)));
            stripInvalidChars(reader->getData(fileIdx, 0), kseq);
            const size_t s = strlen(kseq);
            if (s == 0) {
                Debug(Debug::ERROR) << "Fasta entry " << entries_num << " is invalid\n";
                EXIT(EXIT_FAILURE);
            }

            const char *hdr = hdrReader->getData(fileIdx, 0);
            const size_t hdr_s = hdrReader->getEntryLen(fileIdx);

            // Header creation
            header.append(hdr, hdr_s);

            std::string headerId = Util::parseFastaHeader(header.c_str());
            if (headerId.empty()) {
                // An identifier is necessary for these two cases, so we should just give up
                Debug(Debug::WARNING) << "Cannot extract identifier from entry " << entries_num << newline;
            }

            unsigned int id = par.identifierOffset + entries_num;
            unsigned int splitIdx = id % shuffleSplits;

            /* Write header */
            hdrWriter.writeData(header.c_str(), header.length(), splitIdx, true, true);

            /* e.g. a string of len 12 will be packed into buffer of len 4
                    a string of len 13 need to be packed into buff of len 5
             */
            // calculate remainder
            size_t rem = s % 3;
            // need one more unsigned short if s cannot divide by 3
            int padding = rem == 0 ? 0 : 1;
            rem = (rem == 0) ? 3 : rem;
            unsigned short *resultBuffer = (unsigned short *) calloc((s / 3 + padding),
                                                                     sizeof(unsigned short));
            size_t i = 0;
            // Only go into the for loop if len > 3
            if (s > 3) {
                for (; i < s - 3; i = i + 3) {
                    resultBuffer[i / 3] = PACK_TO_SHORT(kseq[i], kseq[i + 1], kseq[i + 2]);
                }
            }
            /* Write to the last unsigned short */
            // initialize to 0
            // i >= s - 3; i % 3 = 0
            resultBuffer[i / 3] = 0U;

            for (size_t j = 0; j < 3; j++) {
                resultBuffer[i / 3] <<= 5U;
                if (j < rem && kseq[i + j] != '\n') {
                    resultBuffer[i / 3] |= GET_LAST_5_BITS(kseq[i + j]);
                }
            }
            resultBuffer[i / 3] |= 0x8000U; // Set most significant bit
            const char *packedSeq = reinterpret_cast<const char *>(resultBuffer);

            seqWriter.writeStart(splitIdx);
            seqWriter.writeAdd(packedSeq, sizeof(unsigned short) * (s / 3 + padding), splitIdx);
            seqWriter.writeEnd(splitIdx, false);

            entries_num++;
            numEntriesInCurrFile++;
            header.clear();
            free(kseq);
            free(resultBuffer);
        } else {
            // fasta / fastq case
            KSeqWrapper *kseq = KSeqFactory(filenames[fileIdx].c_str());
            while (kseq->ReadEntry()) {
                progress.updateProgress();
                const KSeqWrapper::KSeqEntry &e = kseq->entry;

                if (e.name.l == 0) {
                    Debug(Debug::ERROR) << "Fasta entry " << entries_num << " is invalid\n";
                    EXIT(EXIT_FAILURE);
                }
                // Header creation
                header.append(e.name.s, e.name.l);
                if (e.comment.l > 0) {
                    header.append(" ", 1);
                    header.append(e.comment.s, e.comment.l);
                }

                std::string headerId = Util::parseFastaHeader(header.c_str());
                if (headerId.empty()) {
                    // An identifier is necessary for these two cases, so we should just give up
                    Debug(Debug::WARNING) << "Cannot extract identifier from entry " << entries_num << newline;
                }
                header.push_back(newline);

                unsigned int id = par.identifierOffset + entries_num;
                unsigned int splitIdx = id % shuffleSplits;

                /* Write header */
                hdrWriter.writeData(header.c_str(), header.length(), splitIdx, true, true);

                unsigned long rem = e.sequence.l % 3;
                int padding = rem == 0 ? 0 : 1;
                rem = (rem == 0) ? 3 : rem;

                unsigned short *resultBuffer = (unsigned short *) calloc((e.sequence.l / 3 + padding),
                                                                         sizeof(unsigned short));

                size_t i = 0;
                if (e.sequence.l > 3) {
                    for (; i < e.sequence.l - 3; i = i + 3) {
                        resultBuffer[i / 3] = PACK_TO_SHORT(e.sequence.s[i],
                                                            e.sequence.s[i + 1],
                                                            e.sequence.s[i + 2]);
                    }
                }
                resultBuffer[i / 3] = 0U;
                for (unsigned int j = 0; j < 3; j++) {
                    resultBuffer[i / 3] <<= 5U;
                    if (j < rem && e.sequence.s[i + j] != '\n') {
                        resultBuffer[i / 3] |= GET_LAST_5_BITS(e.sequence.s[i + j]);
                    }
                }
                resultBuffer[i / 3] |= 0x8000U; // Set last bit
                const char *packedSeq = reinterpret_cast<const char *>(resultBuffer);

                seqWriter.writeStart(splitIdx);
                seqWriter.writeAdd(packedSeq, sizeof(unsigned short) * (e.sequence.l / 3 + padding), splitIdx);
                seqWriter.writeEnd(splitIdx, false);

                entries_num++;
                numEntriesInCurrFile++;
                header.clear();
                free(resultBuffer);
            }
            delete kseq;
            kseq = nullptr;
        }
    }

    Debug(Debug::INFO) << newline;
    if (fclose(source) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << sourceFile << newline;
        EXIT(EXIT_FAILURE);
    }

    SRADBWriter::writeDbtypeFile(seqWriter.getDataFileName(), outputDbType, par.compressed);

    hdrWriter.close(true);
    seqWriter.close(true);
    Debug(Debug::INFO) << "Database type: SRA Database" << newline;

    if (isDbInput) {
        reader->close();
        hdrReader->close();
        delete reader;
        reader = nullptr;
        delete hdrReader;
        hdrReader = nullptr;
    }

    if (entries_num == 0) {
        Debug(Debug::ERROR) << "The input files have no entry: ";
        for (size_t fileIdx = 0; fileIdx < filenames.size(); fileIdx++) {
            Debug(Debug::ERROR) << " - " << filenames[fileIdx] << newline;
        }
        Debug(Debug::ERROR) <<
                            "Please check your input files. Only files in fasta/fastq[.gz|bz2] are supported"
                            << newline;
        EXIT(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}

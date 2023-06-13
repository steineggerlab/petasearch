#include "LocalParameters.h"
#include "Debug.h"
#include "NucleotideMatrix.h"
#include "ExtendedSubstitutionMatrix.h"
#include "Indexer.h"
#include "KmerGenerator.h"
#include "MathUtil.h"
#include "QueryTableEntry.h"
#include "DBWriter.h"
#include "BitManipulateMacros.h"
#include "SRAUtil.h"
#include "FastSort.h"

#include <fcntl.h>  // open, read
#include <unistd.h>
#include <cstdlib> // aligned_alloc

#ifdef OPENMP
#include <omp.h>
#endif

#define MEM_SIZE_16MB ((size_t) (16 * 1024 * 1024))
#define MEM_SIZE_32MB ((size_t) (32 * 1024 * 1024))

QueryTableEntry *removeNotHitSequences(
    QueryTableEntry *startPos, QueryTableEntry *endPos, QueryTableEntry *resultTable, unsigned int requiredKmerMatches
) {
    QueryTableEntry *currentReadPos = startPos;
    QueryTableEntry *currentWritePos = resultTable;
    while (currentReadPos < endPos) {
        size_t count = 1;
        while (currentReadPos < endPos
               && currentReadPos->targetSequenceID != UINT_MAX
               && currentReadPos->targetSequenceID == (currentReadPos + 1)->targetSequenceID
               && currentReadPos->querySequenceId == (currentReadPos + 1)->querySequenceId) {
            ++count;
            ++currentReadPos;
        }
        if (count > requiredKmerMatches) {
            memcpy(currentWritePos, currentReadPos - (count - 1), sizeof(QueryTableEntry) * count);
            currentWritePos += count;
        }
        ++currentReadPos;
    }
    return currentWritePos;
}

inline void parallelReadIntoVec(
    int fd,
    std::vector<void *> &destBlocks,
    std::vector<ssize_t> &destBlockSize,
    size_t blockSize,
    bool allocateNewSpace = true,
    size_t offsetBlock = 0
) {
    size_t end = destBlocks.size();
// #pragma omp parallel for schedule(dynamic, 10)
    for (size_t j = 0; j < end; j++) {
        // TODO: determine the alignment dynamically instead of using hard-coded 512
        if (allocateNewSpace) {
            destBlocks[j] = aligned_alloc(512, blockSize);
            if (destBlocks[j] == nullptr) {
                Debug(Debug::ERROR) << "Cannot allocate memory for target table\n";
                EXIT(EXIT_FAILURE);
            }
        }
        off_t offset = (off_t) ((offsetBlock * end + j) * blockSize);
        if ((destBlockSize[j] = pread(fd, destBlocks[j], blockSize, offset)) < 0) {
            Debug(Debug::ERROR) << "Cannot read chunk #" << (j + 1) << " from target table\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

int resultTableSort(const QueryTableEntry &first, const QueryTableEntry &second) {
    if (first.targetSequenceID != second.targetSequenceID) {
        return first.targetSequenceID < second.targetSequenceID;
    }
    if (first.querySequenceId != second.querySequenceId) {
        return first.querySequenceId < second.querySequenceId;
    }
    if (first.Query.kmerPosInQuery != second.Query.kmerPosInQuery) {
        return first.Query.kmerPosInQuery < second.Query.kmerPosInQuery;
    }
    if (first.Query.kmer != second.Query.kmer) {
        return first.Query.kmer < second.Query.kmer;
    }
    return false;
}

int queryTableSort(const QueryTableEntry &first, const QueryTableEntry &second) {
    if (first.Query.kmer != second.Query.kmer) {
        return first.Query.kmer < second.Query.kmer;
    }
    if (first.querySequenceId != second.querySequenceId) {
        return first.querySequenceId > second.querySequenceId;  // Note the '>' operator here
    }
    if (first.Query.kmerPosInQuery != second.Query.kmerPosInQuery) {
        return first.Query.kmerPosInQuery < second.Query.kmerPosInQuery;
    }
    return false;
}

void createQueryTable(LocalParameters &par, std::vector<QueryTableEntry> &queryTable) {
    Timer timer;

    int seqType = FileUtil::parseDbType(par.db1.c_str());

    const bool useProfileSearch = Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_HMM_PROFILE);

    BaseMatrix *subMat;
    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.seedScoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    } else if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_AMINO_ACIDS)) {
        subMat = new SubstitutionMatrix(par.seedScoringMatrixFile.values.aminoacid().c_str(), 8.0, -0.2f);
    } else if (useProfileSearch) {
        subMat = new SubstitutionMatrix(par.seedScoringMatrixFile.values.aminoacid().c_str(), 2.0f, 0.0);
    } else {
        Debug(Debug::ERROR) << "Invalid input type (Support: nucleotide, amino acid, profile)\n";
        EXIT(EXIT_FAILURE);
    }

    DBReader<unsigned int> reader(
        par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA
    );
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    Debug(Debug::INFO) << "Input preparation time: " << timer.lap() << "\n";

    const unsigned int kmerSize = par.kmerSize;
    Debug(Debug::INFO) << "Number of sequences: " << reader.getSize() << "\n";

    // FIXME: should not be hard coded. This is only good for k = 9
    double similarKmerFactor = 1.5 * (useProfileSearch ? 5 : 1);
    size_t tableCapacity = (size_t) (similarKmerFactor * (double) (kmerCount + 1));
    queryTable.reserve(tableCapacity);

    const int xIndex = subMat->aa2num[(int) 'X'];

    ScoreMatrix twoMatrix, threeMatrix;
    if (!useProfileSearch) {
        twoMatrix = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
        threeMatrix = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 3);
    }

    Debug::Progress progress(reader.getSize());

#pragma omp parallel default(none) \
shared(par, reader, subMat, progress, seqType, twoMatrix, threeMatrix, tableCapacity, queryTable, useProfileSearch, xIndex)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        Indexer idx(subMat->alphabetSize - 1, par.kmerSize);
        Sequence sequence(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false,
                          useProfileSearch ? false : true);

        const int kmerThr = useProfileSearch ? par.kmerScore.values.profile() : par.kmerScore.values.sequence();
        KmerGenerator kmerGenerator(kmerSize, subMat->alphabetSize - 1, kmerThr);

        if (useProfileSearch && sequence.profile_matrix != nullptr) {
            kmerGenerator.setDivideStrategy(sequence.profile_matrix);
        } else {
            kmerGenerator.setDivideStrategy(&threeMatrix, &twoMatrix);
        }

        std::vector<QueryTableEntry> localTable;
        localTable.reserve(tableCapacity / total_threads);

#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();
            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, (int) thread_idx);
            unsigned int seqLen = reader.getSeqLen(i);
            sequence.mapSequence(i, key, data, seqLen);

            while (sequence.hasNextKmer()) {
                const unsigned char *kmer = sequence.nextKmer();

                unsigned int xCount = 0;
                for (unsigned int pos = 0; pos < kmerSize; ++pos) {
                    xCount += (kmer[pos] == xIndex);
                }

                if (xCount) {
                    continue;
                }

                if (useProfileSearch) {
                    std::pair<size_t *, size_t> similarKmerList = kmerGenerator.generateKmerList(kmer);
                    size_t lim = similarKmerList.second;
                    for (size_t j = 0; j < lim; ++j) {
                        QueryTableEntry entry{};
                        entry.querySequenceId = key;
                        entry.targetSequenceID = UINT_MAX;
                        entry.Query.kmer = similarKmerList.first[j];
                        entry.Query.kmerPosInQuery = sequence.getCurrentPosition();
                        localTable.emplace_back(entry);
                    }
                    QueryTableEntry entry{};
                    entry.querySequenceId = key;
                    entry.targetSequenceID = UINT_MAX;
                    entry.Query.kmer = idx.int2index(kmer, 0, par.kmerSize);
                    entry.Query.kmerPosInQuery = sequence.getCurrentPosition();
                    localTable.emplace_back(entry);
                } else if (par.exactKmerMatching) {
                    QueryTableEntry entry{};
                    entry.querySequenceId = key;
                    entry.targetSequenceID = UINT_MAX;
                    entry.Query.kmer = idx.int2index(kmer, 0, par.kmerSize);
                    entry.Query.kmerPosInQuery = sequence.getCurrentPosition();
                    localTable.emplace_back(entry);
                } else {
                    // FIXME: too memory consuming when k = 11, need to adjust
                    //  (at least make the program does not terminate with an bad_alloc() error)
                    std::pair<size_t *, size_t> similarKmerList = kmerGenerator.generateKmerList(kmer);
                    size_t lim = similarKmerList.second;
                    for (size_t j = 0; j < lim; ++j) {
                        QueryTableEntry entry{};
                        entry.querySequenceId = key;
                        entry.targetSequenceID = UINT_MAX;
                        entry.Query.kmer = similarKmerList.first[j];
                        entry.Query.kmerPosInQuery = sequence.getCurrentPosition();
                        localTable.emplace_back(entry);
                    }
                }
            }
        }

#pragma omp critical
        queryTable.insert(queryTable.end(), localTable.begin(), localTable.end());
    }
    Debug(Debug::INFO) << "\nk-mers: " << queryTable.size()
                       << "\nk-mers per pos: " << (double) queryTable.size() / (double) reader.getAminoAcidDBSize()
                       << "\nRequired Memory: " << queryTable.size() * sizeof(QueryTableEntry) / 1024 / 1024 << " MB"
                       << "\ntime: " << timer.lap() << "\n";

    timer.reset();
    SORT_PARALLEL(queryTable.begin(), queryTable.end(), queryTableSort);
    Debug(Debug::INFO) << "Sorting time: " << timer.lap() << "\n";

    delete subMat;
    subMat = nullptr;

    if (!useProfileSearch) {
        ExtendedSubstitutionMatrix::freeScoreMatrix(twoMatrix);
        ExtendedSubstitutionMatrix::freeScoreMatrix(threeMatrix);
    }
    reader.close();
}

int comparekmertables(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.spacedKmer = false;

    par.parseParameters(argc, argv, command, true, 0, LocalParameters::PARSE_VARIADIC);

    std::vector<QueryTableEntry> qTable;
    createQueryTable(par, qTable);

    // FIXME: accept single file input also
    std::vector<std::string> targetTables = SRAUtil::getFileNamesFromFile(par.db2);
    std::vector<std::string> resultFiles = SRAUtil::getFileNamesFromFile(par.db3);
    if (targetTables.empty()) {
        Debug(Debug::ERROR) << "Expected at least one targetTable entry in the target table file\n";
        EXIT(EXIT_FAILURE);
    }
    if (targetTables.size() != resultFiles.size()) {
        Debug(Debug::ERROR) << "Number of targetTable and result table is not equal\n";
        EXIT(EXIT_FAILURE);
    }

    const unsigned long queryTableSize = qTable.size() * sizeof(QueryTableEntry);
    const size_t numQTableAvailInMem = Util::getTotalSystemMemory() / 3 / queryTableSize;
    const size_t chunkSize = numQTableAvailInMem >= targetTables.size() ? 1 : targetTables.size() / numQTableAvailInMem;

    const unsigned long MAXIMUM_NUM_OF_BLOCKS =
            (Util::getTotalSystemMemory() - numQTableAvailInMem * queryTableSize) / (MEM_SIZE_16MB + MEM_SIZE_32MB);
    const unsigned long maximumNumOfBlocksPerDB = MAXIMUM_NUM_OF_BLOCKS / targetTables.size();
    size_t localThreads = par.threads;
#ifdef OPENMP
    localThreads = std::max(std::min(localThreads, targetTables.size()), (size_t)1);
#endif

#pragma omp parallel num_threads(localThreads) default(none) shared(par, resultFiles, qTable, targetTables, std::cerr, std::cout, maximumNumOfBlocksPerDB, chunkSize)
    {
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < targetTables.size(); ++i) {
            Timer timer;
            std::vector<QueryTableEntry> localQTable(qTable); //creates a deep copy of the queryTable
            Debug(Debug::INFO) << "Deep copy time: " << timer.lap() << "\n";
            QueryTableEntry *startPosQueryTable = localQTable.data();
            QueryTableEntry *endQueryPos = startPosQueryTable + localQTable.size();

            const std::string& targetName = targetTables[i];

            timer.reset();
#if !defined(O_DIRECT)
            const int mode = (O_RDONLY | O_SYNC);
#else
            const int mode = (O_RDONLY | O_DIRECT | O_SYNC);
#endif
            /* Open target table in direct mode */
            int fdTargetTable = open(targetName.c_str(), mode);
            if (fdTargetTable < 0) {
                Debug(Debug::ERROR) << "Open target table " << targetName << " failed\n";
                EXIT(EXIT_FAILURE);
            }
#if !defined(O_DIRECT) && defined(F_NOCACHE)
            fcntl(fdTargetTable, F_NOCACHE, 1);
#endif

            /* Open ID table in direct mode */
            int fdIDTable = open((targetName + "_ids").c_str(), mode);
            if (fdIDTable < 0) {
                Debug(Debug::ERROR) << "Open ID table " << targetName << "_ids" << "failed\n";
                EXIT(EXIT_FAILURE);
            }
#if !defined(O_DIRECT) && defined(F_NOCACHE)
            fcntl(fdIDTable, F_NOCACHE, 1);
#endif

            /* Get file size in bytes */
            size_t targetTableSize = FileUtil::getFileSize(targetName);
            size_t idTableSize = FileUtil::getFileSize((targetName + "_ids"));

            size_t totalNumOfTargetBlocks = targetTableSize / MEM_SIZE_16MB + (targetTableSize % MEM_SIZE_16MB == 0 ? 0 : 1);
            size_t totalNumOfIDBlocks = idTableSize / MEM_SIZE_32MB + (idTableSize % MEM_SIZE_32MB == 0 ? 0 : 1);

            size_t numOfTargetBlocks = std::min(maximumNumOfBlocksPerDB, totalNumOfTargetBlocks);
            size_t numOfIDBlocks = std::min(maximumNumOfBlocksPerDB, totalNumOfIDBlocks);

            std::vector<void *> targetTableBlocks(numOfTargetBlocks);
            std::vector<ssize_t> targetTableBlockSize(numOfTargetBlocks, -1);
            std::vector<void *> IDTableBlocks(numOfIDBlocks);
            std::vector<ssize_t> IDTableBlockSize(numOfIDBlocks, -1);

            /* Read in 16MB chunks for target table */
            parallelReadIntoVec(fdTargetTable, targetTableBlocks, targetTableBlockSize, MEM_SIZE_16MB);

            /* Read in 32MB chunks for ID table */
            parallelReadIntoVec(fdIDTable, IDTableBlocks, IDTableBlockSize, MEM_SIZE_32MB);

            size_t IDTableIndex = 0;
            size_t IDReadGroup = 0;

            unsigned int *startPosIDTable, *currentIDPos, *endIDPos;
            startPosIDTable = (unsigned int *) IDTableBlocks[IDTableIndex];
            currentIDPos = startPosIDTable;
            endIDPos = startPosIDTable + (MEM_SIZE_32MB / sizeof(unsigned int));

            QueryTableEntry *currentQueryPos = startPosQueryTable;
            QueryTableEntry *endPosQueryTable = startPosQueryTable;

            size_t equalKmers = 0;
            unsigned long long currentKmer = 0;
            bool first = true;
            uint64_t currDiffIndex = 0;

            bool breakOut = false;

            unsigned long long totalBlocksRead = numOfTargetBlocks;
            size_t targetReadGroup = 0;

            while (numOfTargetBlocks == totalNumOfTargetBlocks || totalBlocksRead < totalNumOfTargetBlocks) {
                for (size_t j = 0; j < numOfTargetBlocks; j++) {
                    unsigned short *startPosTargetTable, *endTargetPos, *currentTargetPos;
                    startPosTargetTable = (unsigned short *) targetTableBlocks[j];
                    endTargetPos = startPosTargetTable + (targetTableBlockSize[j] / sizeof(unsigned short));
                    currentTargetPos = startPosTargetTable;
                    // cover the rare case that the first (real) target entry is larger than USHRT_MAX
                    while (currentTargetPos < endTargetPos && !IS_LAST_15_BITS(*currentTargetPos)) {
                        currDiffIndex = DECODE_15_BITS(currDiffIndex, *currentTargetPos);
                        currDiffIndex <<= 15U;
                        ++currentTargetPos;
                    }
                    currDiffIndex = DECODE_15_BITS(currDiffIndex, *currentTargetPos);
                    currentKmer += currDiffIndex;
                    currDiffIndex = 0;

                    while (LIKELY(currentTargetPos < endTargetPos) && currentQueryPos < endQueryPos) {
                        if (currentKmer == currentQueryPos->Query.kmer) {
                            if (first) {
                                startPosQueryTable = currentQueryPos;
                                first = false;
                            }
                            ++equalKmers;
                            currentQueryPos->targetSequenceID = *currentIDPos;
                            ++currentQueryPos;
                            while (LIKELY(currentQueryPos < endQueryPos) &&
                                   currentQueryPos->Query.kmer == currentKmer) {
                                currentQueryPos->targetSequenceID = *currentIDPos;
                                ++currentQueryPos;
                            }
                            endPosQueryTable = currentQueryPos;
                            ++currentTargetPos;
                            ++currentIDPos;
                            if (UNLIKELY(currentIDPos >= endIDPos)) {
                                ++IDTableIndex;
                                if (UNLIKELY(IDTableIndex >= numOfIDBlocks)) {
                                    // parallel read
                                    parallelReadIntoVec(
                                        fdIDTable, IDTableBlocks, IDTableBlockSize, MEM_SIZE_32MB,  false, ++IDReadGroup
                                    );
                                    IDTableIndex = 0;
                                }
                                startPosIDTable = (unsigned int *) IDTableBlocks[IDTableIndex];
                                currentIDPos = startPosIDTable;
                                endIDPos = startPosIDTable + (MEM_SIZE_32MB / sizeof(unsigned int));
                            }
                            while (UNLIKELY(currentTargetPos < endTargetPos && !IS_LAST_15_BITS(*currentTargetPos))) {
                                currDiffIndex = DECODE_15_BITS(currDiffIndex, *currentTargetPos);
                                currDiffIndex <<= 15U;
                                ++currentTargetPos;
                            }
                            if (UNLIKELY(currentTargetPos >= endTargetPos)) {
                                break;
                            }
                            currDiffIndex = DECODE_15_BITS(currDiffIndex, *currentTargetPos);
                            currentKmer += currDiffIndex;
                            currDiffIndex = 0;
                        }

                        while (LIKELY(currentQueryPos < endQueryPos) &&
                               currentQueryPos->Query.kmer < currentKmer) {
                            ++currentQueryPos;
                        }

                        while (currentQueryPos < endQueryPos &&
                               currentTargetPos < endTargetPos &&
                               currentKmer < currentQueryPos->Query.kmer) {
                            ++currentTargetPos;
                            ++currentIDPos;
                            if (UNLIKELY(currentIDPos >= endIDPos)) {
                                ++IDTableIndex;
                                if (UNLIKELY(IDTableIndex >= numOfIDBlocks)) {
                                    // parallel read
                                    parallelReadIntoVec(
                                        fdIDTable, IDTableBlocks, IDTableBlockSize, MEM_SIZE_32MB, false, ++IDReadGroup
                                    );
                                    IDTableIndex = 0;
                                }
                                startPosIDTable = (unsigned int *) IDTableBlocks[IDTableIndex];
                                currentIDPos = startPosIDTable;
                                endIDPos = startPosIDTable + (MEM_SIZE_32MB / sizeof(unsigned int));
                            }
                            while (UNLIKELY(currentTargetPos < endTargetPos && !IS_LAST_15_BITS(*currentTargetPos))) {
                                currDiffIndex = DECODE_15_BITS(currDiffIndex, *currentTargetPos);
                                currDiffIndex <<= 15U;
                                ++currentTargetPos;
                            }
                            if (UNLIKELY(currentTargetPos >= endTargetPos)) {
                                breakOut = true;
                                break;
                            }
                            currDiffIndex = DECODE_15_BITS(currDiffIndex, *currentTargetPos);
                            currentKmer += currDiffIndex;
                            currDiffIndex = 0;
                        }
                        if (UNLIKELY(breakOut)) {
                            breakOut = false;
                            break;
                        }
                    }
                }

                if (numOfTargetBlocks == totalNumOfTargetBlocks) {
                    break;
                }
                parallelReadIntoVec(
                    fdTargetTable, targetTableBlocks, targetTableBlockSize, MEM_SIZE_16MB, false, ++targetReadGroup
                );
                totalBlocksRead += numOfTargetBlocks;
            }

            double timediff = timer.getTimediff();
            Debug(Debug::INFO) << timediff << " s; Rate "
                               << ((double) (targetTableSize + idTableSize) / 1e+9) / timediff << " GB/s \n";
            Debug(Debug::INFO) << "Number of equal k-mers: " << equalKmers << "\n";

            for (size_t j = 0; j < numOfTargetBlocks; j++) {
                free(targetTableBlocks[j]);
            }
            for (size_t j = 0; j < numOfIDBlocks; j++) {
                free(IDTableBlocks[j]);
            }

            if (close(fdIDTable) < 0) {
                Debug(Debug::ERROR) << "Cannot close ID table\n";
                EXIT(EXIT_FAILURE);
            }

            if (close(fdTargetTable) < 0) {
                Debug(Debug::ERROR) << "Cannot close target table\n";
                EXIT(EXIT_FAILURE);
            }

            timer.reset();
            SORT_SERIAL(startPosQueryTable, endQueryPos, resultTableSort);
            Debug(Debug::INFO) << "Result table sort time: " << timer.lap() << "\n";
            timer.reset();

            QueryTableEntry *resultTable = new QueryTableEntry[endPosQueryTable - startPosQueryTable + 1];
            QueryTableEntry *truncatedResultEndPos = removeNotHitSequences(startPosQueryTable, endQueryPos, resultTable, par.requiredKmerMatches);
            Debug(Debug::INFO) << "Duplicate elimination time: " << timer.lap() << "\n";
            timer.reset();

            std::string resultDB = resultFiles[i];
            DBWriter writer(
                resultDB.c_str(), (resultDB + ".index").c_str(), 1, par.compressed, Parameters::DBTYPE_PREFILTER_RES
            );
            writer.open();

            QueryTableEntry *startPos = resultTable;
            QueryTableEntry *endPos = truncatedResultEndPos;
            std::string result;
            result.reserve(10 * 1024);
            char buffer[1024];
            for (QueryTableEntry *currentPos = startPos; currentPos < endPos - 1; ++currentPos) {
//        bool didAppend = false;
                while (currentPos < endPos - 1 && currentPos->targetSequenceID == (currentPos + 1)->targetSequenceID) {
                    size_t len = QueryTableEntry::queryEntryToBuffer(buffer, *currentPos);
//            if (didAppend == false) {
                    result.append(buffer, len);
//                didAppend = true;
//            }
                    ++currentPos;
                }
                writer.writeData(result.c_str(), result.length(), currentPos->targetSequenceID, 0);
                result.clear();
            }
            writer.close();
            Debug(Debug::INFO) << "Result write time: " << timer.lap() << "\n";

            delete[] resultTable;
            resultTable = nullptr;
        }
    }

    return EXIT_SUCCESS;
}

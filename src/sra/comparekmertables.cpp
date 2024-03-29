#include "LocalParameters.h"
#include "Debug.h"
#include "NucleotideMatrix.h"
#include "ExtendedSubstitutionMatrix.h"
#include "Indexer.h"
#include "KmerGenerator.h"
#include "FixedKmerGenerator.h"
#include "MathUtil.h"
#include "QueryTableEntry.h"
#include "DBWriter.h"
#include "BitManipulateMacros.h"
#include "SRAUtil.h"
#include "FastSort.h"
#include "tantan.h"

#include <map>
#include <fcntl.h>  // open, read
#include <unistd.h>
#include <cstdlib> // aligned_alloc
#include <sys/types.h>
#include <sys/stat.h>

#ifdef OPENMP
#include <omp.h>
#endif

#define MEM_SIZE_16MB ((size_t) (16 * 1024 * 1024))
#define MEM_SIZE_32MB ((size_t) (32 * 1024 * 1024))

QueryTableEntry *copyHitSequences(
    QueryTableEntry *startPos, QueryTableEntry *endPos, QueryTableEntry *destPos
) {
    QueryTableEntry *currentReadPos = startPos;
    QueryTableEntry *currentWritePos = destPos;
    while (currentReadPos < endPos) {
        if (currentReadPos->targetSequenceID != UINT_MAX) {
            // If this entry had a hit, we copy it to the destination location
            memcpy(currentWritePos, currentReadPos, sizeof(QueryTableEntry));
            ++currentWritePos;
        }
        ++currentReadPos;
    }
    return currentWritePos;
}

QueryTableEntry *removeNotHitSequences(
    QueryTableEntry *startPos, QueryTableEntry *endPos, unsigned int requiredKmerMatches
) {
    QueryTableEntry *currentReadPos = startPos;
    QueryTableEntry *currentWritePos = startPos;
    while (currentReadPos < endPos) {
        size_t count = 1;
        while (currentReadPos < endPos - 1
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

    // FIXME: incorporate previous
    // double similarKmerFactor = 1.5 * (useProfileSearch ? 5 : 1);
    // size_t tableCapacity = (size_t) (similarKmerFactor * (double) (kmerCount + 1));
    const size_t kmerCount = reader.getAminoAcidDBSize() - (reader.getSize() * (par.kmerSize + 1));
    const size_t tableCapacity = (size_t) (par.maxKmerPerPos * (kmerCount + 1));
    queryTable.reserve(tableCapacity);

    const int xIndex = subMat->aa2num[(int) 'X'];

    ScoreMatrix twoMatrix, threeMatrix;
    if (!useProfileSearch) {
        twoMatrix = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
        threeMatrix = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 3);
    }

    ProbabilityMatrix *probMatrix = NULL;
    if (par.maskMode == 1) {
        probMatrix = new ProbabilityMatrix(*subMat);
    }
    const unsigned int total_threads = par.threads;

    Debug::Progress progress(reader.getSize());
#pragma omp parallel default(none) shared(par, reader, subMat, progress, seqType, twoMatrix, threeMatrix, tableCapacity, queryTable, useProfileSearch, probMatrix) firstprivate(xIndex, kmerSize, total_threads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        Indexer idx(subMat->alphabetSize - 1, kmerSize);
        Sequence sequence(par.maxSeqLen, seqType, subMat, kmerSize, par.spacedKmer, par.compBiasCorrection, useProfileSearch ? false : true, par.spacedKmerPattern);

        const int kmerThr = useProfileSearch ? par.kmerScore.values.profile() : par.kmerScore.values.sequence();
        FixedKmerGenerator kmerGenerator(kmerSize, subMat->alphabetSize - 1, kmerThr, par.maxKmerPerPos);

        if (useProfileSearch && sequence.profile_matrix != nullptr) {
            kmerGenerator.setDivideStrategy(sequence.profile_matrix);
        } else {
            kmerGenerator.setDivideStrategy(&threeMatrix, &twoMatrix);
        }

        std::vector<QueryTableEntry> localTable;
        localTable.reserve(tableCapacity / total_threads);

        float *compositionBias = nullptr;
        if (par.compBiasCorrection == 1) {
            compositionBias = new float[par.maxSeqLen];
        }

#pragma omp for schedule(dynamic, 1) nowait
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();
            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, (int) thread_idx);
            unsigned int seqLen = reader.getSeqLen(i);
            sequence.mapSequence(i, key, data, seqLen);

            if (par.compBiasCorrection == 1) {
                SubstitutionMatrix::calcLocalAaBiasCorrection(subMat, sequence.numSequence, sequence.L, compositionBias, par.compBiasCorrectionScale);
            }

            if (par.maskMode == 1) {
                tantan::maskSequences(
                    (char*)sequence.numSequence,
                    (char*)(sequence.numSequence + sequence.L),
                    50 /*options.maxCycleLength*/,
                    probMatrix->probMatrixPointers,
                    0.005 /*options.repeatProb*/,
                    0.05 /*options.repeatEndProb*/,
                    0.5 /*options.repeatOffsetProbDecay*/,
                    0, 0,
                    par.maskProb /*options.minMaskProb*/,
                    probMatrix->hardMaskTable
                );
                const char* charSeq = sequence.getSeqData();
                for (int i = 0; i < sequence.L; i++) {
                    sequence.numSequence[i] = (islower(charSeq[i])) ? xIndex : sequence.numSequence[i];
                }
            }

            while (sequence.hasNextKmer()) {
                const unsigned char *kmer = sequence.nextKmer();
                if (sequence.kmerContainsX()) {
                    continue;
                }

                if (par.compBiasCorrection == 1) {
                    const unsigned char *pos = sequence.getAAPosInSpacedPattern();
                    const unsigned short current_i = sequence.getCurrentPosition();
                    float biasCorrection = 0;
                    for (unsigned int i = 0; i < kmerSize; i++) {
                        biasCorrection += compositionBias[current_i + static_cast<short>(pos[i])];
                    }
                    short bias = std::min((short)0, static_cast<short>((biasCorrection < 0.0) ? biasCorrection - 0.5 : biasCorrection + 0.5));
                    short kmerMatchScore = std::max(kmerThr - bias, 0);

                    // Debug(Debug::ERROR) << "bias: " << bias << " kmerMatchScore: " << kmerMatchScore << "\n";

                    // adjust kmer threshold based on composition bias
                    kmerGenerator.setThreshold(kmerMatchScore);
                }

                QueryTableEntry entry{};
                entry.querySequenceId = key;
                entry.targetSequenceID = UINT_MAX;
                entry.Query.kmer = idx.int2index(kmer, 0, kmerSize);
                // idx.printKmer(entry.Query.kmer, kmerSize, subMat->num2aa);
                // Debug(Debug::INFO) << "\n";
                entry.Query.kmerPosInQuery = sequence.getCurrentPosition();
                localTable.emplace_back(entry);
                if (par.exactKmerMatching == false) {
                    std::pair<size_t *, size_t> similarKmerList = kmerGenerator.generateKmerList(kmer); // , false, par.maxKmerPerPos);
                    for (size_t j = 0; j < similarKmerList.second; ++j) {
                        entry.querySequenceId = key;
                        entry.targetSequenceID = UINT_MAX;
                        entry.Query.kmer = similarKmerList.first[j];
                        // idx.printKmer(entry.Query.kmer, kmerSize, subMat->num2aa);
                        // Debug(Debug::INFO) << "\n";
                        entry.Query.kmerPosInQuery = sequence.getCurrentPosition();
                        localTable.emplace_back(entry);
                    }
                }
            }
        }

        if (par.compBiasCorrection == 1) {
            delete[] compositionBias;
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

std::vector<size_t> roundRobinOrder(const std::vector<std::string>& filenames) {
    std::multimap<long, std::pair<size_t, std::string>> deviceToFileNames;
    std::vector<size_t> result;

    // Fill the multimap with filenames keyed by device identifier
    for (size_t i = 0; i < filenames.size(); i++) {
        struct stat sb;
        if (stat(filenames[i].c_str(), &sb) == 0) {
            deviceToFileNames.emplace(sb.st_dev, std::make_pair(i, filenames[i]));
        }
    }

    // Iterate over the multimap in round-robin order
    while (!deviceToFileNames.empty()) {
        for (auto it = deviceToFileNames.begin(); it != deviceToFileNames.end();) {
            // Add the index to the result
            result.emplace_back(it->second.first);

            // Store the iterator pointing to the next device
            auto nextIt = it;
            ++nextIt;

            // Remove the current filename from the multimap
            deviceToFileNames.erase(it);

            // Move to the next device
            it = nextIt;
        }
    }

    return result;
}

template<typename T>
void reorderVectorInPlace(std::vector<T>& original, const std::vector<size_t>& indices) {
    std::vector<T> temp = original;
    for (size_t i = 0; i < indices.size(); i++) {
        original[i] = temp[indices[i]];
    }
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

    std::vector<size_t> indices = roundRobinOrder(targetTables);
    reorderVectorInPlace(targetTables, indices);
    reorderVectorInPlace(resultFiles, indices);

    const unsigned long queryTableSize = qTable.size() * sizeof(QueryTableEntry);
    const size_t numQTableAvailInMem = Util::getTotalSystemMemory() / 3 / queryTableSize;
    // const size_t chunkSize = numQTableAvailInMem >= targetTables.size() ? 1 : targetTables.size() / numQTableAvailInMem;

    const unsigned long MAXIMUM_NUM_OF_BLOCKS =
            (Util::getTotalSystemMemory() - numQTableAvailInMem * queryTableSize) / (MEM_SIZE_16MB + MEM_SIZE_32MB);
    const unsigned long maximumNumOfBlocksPerDB = MAXIMUM_NUM_OF_BLOCKS / targetTables.size();

    size_t localThreads = par.threads;
#ifdef OPENMP
    localThreads = std::max(std::min(localThreads, targetTables.size()), (size_t)1);
#endif

#pragma omp parallel num_threads(localThreads) default(none) shared(par, resultFiles, qTable, targetTables, std::cerr, std::cout, maximumNumOfBlocksPerDB)
    {
        Timer timer;
        QueryTableEntry* localQTable = new QueryTableEntry[qTable.size()];
        std::memcpy(localQTable, qTable.data(), qTable.size() * sizeof(QueryTableEntry));
        Debug(Debug::INFO) << "Deep copy time: " << timer.lap() << "\n";

        QueryTableEntry* resultTable = new QueryTableEntry[qTable.size()];

        std::string result;
        result.reserve(10 * 1024 * 1024);
        char buffer[1024];

#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < targetTables.size(); ++i) {
            // std::vector<QueryTableEntry> localQTable(qTable); //creates a deep copy of the queryTable
            QueryTableEntry *startPosQueryTable = localQTable;
            QueryTableEntry *endQueryPos = startPosQueryTable + qTable.size();

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
            QueryTableEntry *resultTableEndPos = copyHitSequences(startPosQueryTable, endPosQueryTable, resultTable);
            Debug(Debug::INFO) << "Not hit kmer elimination time: " << timer.lap() << "\n";
            timer.reset();
            
            timer.reset();
            SORT_SERIAL(resultTable, resultTableEndPos, resultTableSort);
            Debug(Debug::INFO) << "Result table sort time: " << timer.lap() << "\n";
            timer.reset();

            QueryTableEntry *truncatedResultEndPos = removeNotHitSequences(resultTable, resultTableEndPos, par.requiredKmerMatches);
            Debug(Debug::INFO) << "Duplicate elimination time: " << timer.lap() << "\n";
            timer.reset();
            Debug(Debug::INFO) << "Reduced k-mers " << (endPosQueryTable - startPosQueryTable) << " -> " << (resultTableEndPos - resultTable) << " -> " << (truncatedResultEndPos - resultTable) << "\n";

            std::string resultDB = resultFiles[i];
            DBWriter writer(
                resultDB.c_str(), (resultDB + ".index").c_str(), 1, par.compressed, Parameters::DBTYPE_PREFILTER_RES
            );
            writer.open();

            for (QueryTableEntry *currentPos = resultTable; currentPos < truncatedResultEndPos - 1; ++currentPos) {
//        bool didAppend = false;
                while (currentPos < truncatedResultEndPos - 1 && currentPos->targetSequenceID == (currentPos + 1)->targetSequenceID) {
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
            result.clear();
            Debug(Debug::INFO) << "Result write time: " << timer.lap() << "\n";

        }
        delete[] localQTable;
        delete[] resultTable;
    }

    return EXIT_SUCCESS;
}

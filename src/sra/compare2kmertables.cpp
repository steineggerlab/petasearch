#include "LocalParameters.h"
#include "Command.h"
#include "Debug.h"
#include "NucleotideMatrix.h"
#include "ExtendedSubstitutionMatrix.h"
#include "Indexer.h"
#include "KmerGenerator.h"
#include "MathUtil.h"
#include "QueryTableEntry.h"
#include "DBWriter.h"
#include "MemoryMapped.h"

#include "omptl/omptl_algorithm"

#ifdef OPENMP
#include <omp.h>
#endif

int resultTableSort(const QueryTableEntry &first, const QueryTableEntry &second);
int truncatedResultTableSort(const QueryTableEntry &first, const QueryTableEntry &second);
QueryTableEntry *removeNotHitSequences(QueryTableEntry *startPos, QueryTableEntry *endPos, QueryTableEntry *resultTable, LocalParameters &par);

int queryTableSort(const QueryTableEntry &first, const QueryTableEntry &second) {
    if (first.Query.kmer < second.Query.kmer) {
        return true;
    }
    if (second.Query.kmer < first.Query.kmer) {
        return false;
    }
    if (first.querySequenceId > second.querySequenceId) {
        return true;
    }
    if (second.querySequenceId > first.querySequenceId) {
        return false;
    }
    if (first.Query.kmerPosInQuery < second.Query.kmerPosInQuery) {
        return true;
    }
    if (second.Query.kmerPosInQuery < first.Query.kmerPosInQuery) {
        return false;
    }
    return false;
}


void createQueryTable(LocalParameters &par, std::vector<QueryTableEntry> &queryTable){
    Timer timer;
    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    BaseMatrix *subMat;
    int seqType = reader.getDbtype();
    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.seedScoringMatrixFile.nucleotides, 1.0, 0.0);
    } else {
        subMat = new SubstitutionMatrix(par.seedScoringMatrixFile.aminoacids, 8.0, -0.2f);
    }
    Debug(Debug::INFO) << "input prepared, time spent: " << timer.lap() << "\n";
    size_t kmerCount = 0;
#pragma omp parallel for reduction(+:kmerCount)
    for (size_t i = 0; i < reader.getSize(); ++i) {
        size_t currentSequenceLength = reader.getSeqLen(i);
        //number of ungapped k-mers per sequence = seq.length-k-mer.size+1
        kmerCount += currentSequenceLength >= (size_t)par.kmerSize ? currentSequenceLength - par.kmerSize + 1 : 0;
    }
    Debug(Debug::INFO) << "Number of sequences: " << reader.getSize() << "\n";
    float similarKmerFactor = 1.5 ;
    size_t tableCapacity = (size_t) (similarKmerFactor * (kmerCount + 1));
    queryTable.reserve(tableCapacity);

    const int xIndex = subMat->aa2int[(int)'X'];

    Debug::Progress progress(reader.getSize());

    ScoreMatrix twoMatrix = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
    ScoreMatrix threeMatrix = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 3);

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
        unsigned int total_threads = 1;
#ifdef OPENMP
        thread_idx = (unsigned int)omp_get_thread_num();
        total_threads = (unsigned int)omp_get_num_threads();
#endif
        Indexer idx(subMat->alphabetSize - 1, par.kmerSize);
        Sequence sequence(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false);
        KmerGenerator kmerGenerator(par.kmerSize, subMat->alphabetSize - 1, par.kmerScore);
        kmerGenerator.setDivideStrategy(&threeMatrix, &twoMatrix);

        std::vector<QueryTableEntry> localTable;
        localTable.reserve(tableCapacity / total_threads);

#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();
            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, thread_idx);
            unsigned int seqLen = reader.getSeqLen(i);
            sequence.mapSequence(i, 0, data, seqLen);

            while (sequence.hasNextKmer()) {
                const int *kmer = sequence.nextKmer();

                int xCount = 0;
                for (int pos = 0; pos < par.kmerSize; ++pos) {
                    xCount += (kmer[pos] == xIndex);
                }
                if (xCount) {
                    continue;
                }
                if (par.exactKmerMatching) {
                    QueryTableEntry entry;
                    entry.querySequenceId = key;
                    entry.targetSequenceID = UINT_MAX;
                    entry.Query.kmer = idx.int2index(kmer, 0, par.kmerSize);
                    entry.Query.kmerPosInQuery = sequence.getCurrentPosition();
                    localTable.emplace_back(entry);
                } else {
                    std::pair<size_t*, size_t> similarKmerList = kmerGenerator.generateKmerList(kmer);
                    for (size_t j = 0; j < similarKmerList.second; ++j) {
                        QueryTableEntry entry;
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
                       << "\nRequired Memory: " << queryTable.size()  * sizeof(QueryTableEntry) / 1024 / 1024 << " MB"
                       << "\ntime: " << timer.lap() << "\n";

    Debug(Debug::INFO) << "start sorting \n";
    omptl::sort(queryTable.begin(), queryTable.end(), queryTableSort);
    Debug(Debug::INFO) << "Required time for sorting" << timer.lap() << "\n";

    reader.close();
}


int compare2kmertables(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, true, 0, LocalParameters::PARSE_VARIADIC);

    Debug(Debug::INFO) << "mapping query and target files \n";
    std::vector<QueryTableEntry> qTable;
    createQueryTable(par, qTable);

    QueryTableEntry *startPosQueryTable = qTable.data();
    QueryTableEntry *endQueryPos = startPosQueryTable + qTable.size();

    MemoryMapped targetTable(par.db2, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    MemoryMapped targetIds(std::string(par.db2 + "_ids"), MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    if (targetTable.isValid() == false || targetIds.isValid() == false) {
        Debug(Debug::ERROR) << "Could not open target database " << par.db2 << "\n";
        EXIT(EXIT_FAILURE);
    }
    unsigned short *startPosTargetTable = (unsigned short *)targetTable.getData();
    unsigned int *startPosIDTable = (unsigned int *)targetIds.getData();

    QueryTableEntry *currentQueryPos = startPosQueryTable;
    unsigned short *currentTargetPos = startPosTargetTable;
    unsigned short *endTargetPos = startPosTargetTable + targetTable.size() / sizeof(unsigned short);
    unsigned int *currentIDPos = startPosIDTable;

    size_t equalKmers = 0;
    size_t currentKmer = 0;


    Debug(Debug::INFO) << "start comparing \n";

    Timer timer;
    // cover the rare case that the first (real) target entry is larger than USHRT_MAX
    while (currentTargetPos <= endTargetPos && *currentTargetPos == USHRT_MAX) {
        currentKmer += USHRT_MAX;
        ++currentTargetPos;
        ++currentIDPos;
    }
    currentKmer += *currentTargetPos;

    while (__builtin_expect(currentTargetPos < endTargetPos, 1) && currentQueryPos < endQueryPos) {
        if (currentKmer == currentQueryPos->Query.kmer) {
            ++equalKmers;
            currentQueryPos->targetSequenceID = *currentIDPos;
            ++currentQueryPos;
            while (__builtin_expect(currentQueryPos < endQueryPos, 1) && currentQueryPos->Query.kmer == currentKmer){
                currentQueryPos->targetSequenceID = *currentIDPos;
                ++currentQueryPos;
            }
            ++currentTargetPos;
            ++currentIDPos;
            while (__builtin_expect(*currentTargetPos == USHRT_MAX && currentTargetPos < endTargetPos, 0)) {
                currentKmer += USHRT_MAX;
                ++currentTargetPos;
                ++currentIDPos;
            }
            currentKmer += *currentTargetPos;
        }

        while (currentQueryPos->Query.kmer < currentKmer && __builtin_expect(currentQueryPos < endQueryPos, 1)) {
            ++currentQueryPos;
        }

        while (currentKmer < currentQueryPos->Query.kmer && currentTargetPos < endTargetPos) {
            ++currentTargetPos;
            ++currentIDPos;
            while (__builtin_expect(*currentTargetPos == USHRT_MAX && currentTargetPos < endTargetPos, 0)) {
                currentKmer += USHRT_MAX;
                ++currentTargetPos;
                ++currentIDPos;
            }
            currentKmer += *currentTargetPos;
        }
    }


    double timediff = timer.getTimediff();
    Debug(Debug::INFO) << timediff << " s; Rate " << ((targetTable.size() + targetIds.size()) / 1e+9) / timediff << " GB/s \n";
    Debug(Debug::INFO) << "Number of equal k-mers: " << equalKmers << "\n";

    Debug(Debug::INFO) << "Sorting result table\n";
    omptl::sort(startPosQueryTable, endQueryPos, resultTableSort);

    Debug(Debug::INFO) << "Removing sequences with less than two hits\n";
    QueryTableEntry *resultTable = new QueryTableEntry[qTable.size()];
    QueryTableEntry *truncatedResultEndPos = removeNotHitSequences(startPosQueryTable, endQueryPos, resultTable, par);

    Debug(Debug::INFO) << "Sorting result table after target id \n";
    omptl::sort(resultTable, truncatedResultEndPos, truncatedResultTableSort);

    Debug(Debug::INFO) << "Writing result files\n";
    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed, Parameters::DBTYPE_PREFILTER_RES);
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

    delete[] resultTable;

    return EXIT_SUCCESS;
}

QueryTableEntry *removeNotHitSequences(QueryTableEntry *startPos, QueryTableEntry *endPos, QueryTableEntry *resultTable, LocalParameters &par) {
    QueryTableEntry *currentReadPos = startPos;
    QueryTableEntry *currentWritePos = resultTable;
    while (currentReadPos < endPos) {
        size_t count = 1;
        while (currentReadPos < endPos
                && currentReadPos->targetSequenceID != UINT_MAX
                && currentReadPos->targetSequenceID == (currentReadPos + 1)->targetSequenceID
                && currentReadPos->querySequenceId  == (currentReadPos + 1)->querySequenceId) {
            ++count;
            ++currentReadPos;
        }
        if (count > par.requiredKmerMatches) {
            memcpy(currentWritePos, currentReadPos - (count - 1), sizeof(QueryTableEntry) * count);
            currentWritePos += count;
        }
        ++currentReadPos;
    }
    return currentWritePos;
}

int resultTableSort(const QueryTableEntry &first, const QueryTableEntry &second) {
    if (first.querySequenceId < second.querySequenceId) {
        return true;
    }
    if (second.querySequenceId < first.querySequenceId) {
        return false;
    }
    if (first.targetSequenceID < second.targetSequenceID) {
        return true;
    }
    if (second.targetSequenceID < first.targetSequenceID) {
        return false;
    }
    if (first.Query.kmerPosInQuery < second.Query.kmerPosInQuery) {
        return true;
    }
    if (second.Query.kmerPosInQuery < first.Query.kmerPosInQuery) {
        return false;
    }
    if (first.Query.kmer < second.Query.kmer) {
        return true;
    }
    if (second.Query.kmer < first.Query.kmer) {
        return false;
    }
    return false;
}

int truncatedResultTableSort(const QueryTableEntry &first, const QueryTableEntry &second) {
    if (first.targetSequenceID < second.targetSequenceID) {
        return true;
    }
    if (second.targetSequenceID < first.targetSequenceID) {
        return false;
    }
    if (first.querySequenceId < second.querySequenceId) {
        return true;
    }
    if (second.querySequenceId < first.querySequenceId) {
        return false;
    }
    if (first.Query.kmerPosInQuery < second.Query.kmerPosInQuery) {
        return true;
    }
    if (second.Query.kmerPosInQuery < first.Query.kmerPosInQuery) {
        return false;
    }
    if (first.Query.kmer < second.Query.kmer) {
        return true;
    }
    if (second.Query.kmer < first.Query.kmer) {
        return false;
    }
    return false;
}
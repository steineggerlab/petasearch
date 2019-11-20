#include "LocalParameters.h"
#include "Command.h"
#include "Debug.h"
#include <sys/stat.h>
#include <sys/mman.h>
#include "MathUtil.h"
#include "QueryTableEntry.h"
#include "omptl/omptl_algorithm"
#include "DBWriter.h"

int resultTableSort(const QueryTableEntry &first, const QueryTableEntry &second);
void writeResultTable(QueryTableEntry *startPos, QueryTableEntry *endPos, Parameters &par);
int truncatedResultTableSort(const QueryTableEntry &first, const QueryTableEntry &second);
QueryTableEntry *removeNotHitSequences(QueryTableEntry *startPos, QueryTableEntry *endPos, QueryTableEntry *resultTable, LocalParameters &par);

int compare2kmertables(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, true, 0, 0);

    Debug(Debug::INFO) << "mapping query and target files \n";
    FILE *handleQueryKmerTable = fopen(par.db1.c_str(), "rb");
    int fdQueryTable = fileno(handleQueryKmerTable);
    struct stat fileStatsQueryTable;
    fstat(fdQueryTable, &fileStatsQueryTable);
    size_t fileSizeQueryTable = fileStatsQueryTable.st_size;

    FILE *handleTargetKmerTable = fopen(par.db2.c_str(), "rb");
    int fdTargetTable = fileno(handleTargetKmerTable);
    struct stat fileStatsTargetTable;
    fstat(fdTargetTable, &fileStatsTargetTable);
    size_t fileSizeTargetTable = fileStatsTargetTable.st_size;

    FILE *handleTargetIDTable = fopen(par.db3.c_str(), "rb");
    int fdTargetIDTable = fileno(handleTargetIDTable);
    struct stat fileStatsTargetIDTable;
    fstat(fdTargetIDTable, &fileStatsTargetIDTable);
    size_t fileSizeTargetIDTable = fileStatsTargetIDTable.st_size;

    QueryTableEntry *startPosQueryTable = (QueryTableEntry *)mmap(NULL, fileSizeQueryTable, PROT_READ | PROT_WRITE, MAP_PRIVATE, fdQueryTable, 0);
    unsigned short *startPosTargetTable = (unsigned short *)mmap(NULL, fileSizeTargetTable, PROT_READ, MAP_PRIVATE, fdTargetTable, 0);
    unsigned int *startPosIDTable = (unsigned int *)mmap(NULL, fileSizeTargetIDTable, PROT_READ, MAP_PRIVATE, fdTargetIDTable, 0);
    if (posix_madvise(startPosQueryTable, fileSizeQueryTable, POSIX_MADV_SEQUENTIAL | POSIX_MADV_WILLNEED) != 0) {
        Debug(Debug::ERROR) << "posix_madvise returned an error for the query k-mer table\n";
    }
    if (posix_madvise(startPosTargetTable, fileSizeTargetTable, POSIX_MADV_SEQUENTIAL | POSIX_MADV_WILLNEED) != 0) {
        Debug(Debug::ERROR) << "posix_madvise returned an error for the target k-mer table\n";
    }
    if (posix_madvise(startPosIDTable, fileSizeTargetIDTable, POSIX_MADV_SEQUENTIAL | POSIX_MADV_WILLNEED) != 0) {
        Debug(Debug::ERROR) << "posix_madvise returned an error for the target id table\n";
    }

    QueryTableEntry *currentQueryPos = startPosQueryTable;
    QueryTableEntry *endQueryPos = startPosQueryTable + fileSizeQueryTable / sizeof(QueryTableEntry);
    unsigned short *currentTargetPos = startPosTargetTable;
    unsigned short *endTargetPos = startPosTargetTable + fileSizeTargetTable / sizeof(unsigned short);
    unsigned int *currentIDPos = startPosIDTable;
    size_t equalKmers = 0;
    size_t currentKmer = 0;
    struct timeval startTime;
    struct timeval endTime;
    gettimeofday(&startTime, NULL);

    Debug(Debug::INFO) << "start comparing \n";

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

    gettimeofday(&endTime, NULL);
    double timediff = (endTime.tv_sec - startTime.tv_sec) + 1e-6 * (endTime.tv_usec - startTime.tv_usec);
    Debug(Debug::INFO) << timediff << " s; Rate " << ((fileSizeTargetTable + fileSizeQueryTable + fileSizeTargetIDTable) / 1e+9) / timediff
                       << " GB/s \n";
    Debug(Debug::INFO) << "Number of equal k-mers: " << equalKmers << "\n";

    Debug(Debug::INFO) << "Sorting result table\n";
    omptl::sort(startPosQueryTable, endQueryPos, resultTableSort);

    Debug(Debug::INFO) << "Removing sequences with less than two hits \n";
    QueryTableEntry *resultTable = (QueryTableEntry *)malloc(fileSizeQueryTable);

    QueryTableEntry *truncatedResultEndPos = removeNotHitSequences(startPosQueryTable, endQueryPos, resultTable, par);

    Debug(Debug::INFO) << "Sorting result table after target id  \n";
    omptl::sort(resultTable, truncatedResultEndPos, truncatedResultTableSort);
    Debug(Debug::INFO) << "Writing result files \n";
    writeResultTable(resultTable, truncatedResultEndPos, par);

    munmap(startPosQueryTable, fileSizeQueryTable);
    munmap(startPosTargetTable, fileSizeTargetTable);
    munmap(startPosIDTable, fileSizeTargetIDTable);
    fclose(handleQueryKmerTable);
    fclose(handleTargetKmerTable);
    fclose(handleTargetIDTable);

    return 0;
}

void writeResultTable(QueryTableEntry *startPos, QueryTableEntry *endPos, Parameters &par) {
    DBWriter *writer = new DBWriter(par.db4.c_str(), par.db4Index.c_str(), 1, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    writer->open();
    for (QueryTableEntry *currentPos = startPos; currentPos < endPos; ++currentPos) {
        size_t blockSize = 0;
        while (currentPos < endPos && currentPos->targetSequenceID == (currentPos + 1)->targetSequenceID) {
            ++blockSize;
            ++currentPos;
        }
        writer->writeData((char *)(currentPos - blockSize), blockSize * sizeof(QueryTableEntry), currentPos->targetSequenceID, 0U, false);
    }
    writer->close();
}

QueryTableEntry *removeNotHitSequences(QueryTableEntry *startPos, QueryTableEntry *endPos, QueryTableEntry *resultTable, LocalParameters &par) {
    QueryTableEntry *currentReadPos = startPos;
    QueryTableEntry *currentWritePos = resultTable;
    while (currentReadPos < endPos) {
        size_t count = 0;
        while (currentReadPos < endPos
                && currentReadPos->targetSequenceID != 0
                && currentReadPos->targetSequenceID == (currentReadPos + 1)->targetSequenceID
                && currentReadPos->querySequenceId  == (currentReadPos + 1)->querySequenceId) {
            count++;
            ++currentReadPos;
        }
        if (count >= par.requiredKmerMatches) {
            memcpy(currentWritePos, currentReadPos - count, sizeof(QueryTableEntry) * count);
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
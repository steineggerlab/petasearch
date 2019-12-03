#include "LocalParameters.h"
#include "Command.h"
#include "Debug.h"
#include "IndexReader.h"
#include "DBReader.h"
#include "NucleotideMatrix.h"
#include "MathUtil.h"
#include "QueryTableEntry.h"
#include "TargetTableEntry.h"
#include "omptl/omptl_algorithm"
#include <sys/mman.h>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

#define BUFFERSIZE Util::getPageSize()

void writeQueryTable(QueryTableEntry *queryTable, size_t kmerCount, std::string queryID);
void writeTargetTables(TargetTableEntry *targetTable, size_t kmerCount, std::string blockID);
int queryTableSort(const QueryTableEntry &first, const QueryTableEntry &second);
int targetTableSort(const TargetTableEntry &first, const TargetTableEntry &second);
size_t countKmer(DBReader<unsigned int> *reader, unsigned int kmerSize);
int xCountInSequence(const int *kmer, size_t kmerSize, const int xIndex);
int createQueryTable(Parameters &par, DBReader<unsigned int> *reader, BaseMatrix *subMat);
int createTargetTable(Parameters &par, DBReader<unsigned int> *reader, BaseMatrix *subMat);
void writeKmerDiff(size_t lastKmer, TargetTableEntry *entryToWrite, FILE *handleKmerTable, FILE *handleIDTable);

size_t  diffLargerThenUShortMax = 0;
size_t entryDiffLargerUShortMax =0;

int createkmertable(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, true, 0, 0);
    Timer timer;
    Debug(Debug::INFO) << "Preparing input database\n";

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    BaseMatrix *subMat;
    int seqType = reader.getDbtype();
    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.nucleotides, 1.0, 0.0);
    }
    else {
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.aminoacids, 2.0, 0.0);
    }
    Debug(Debug::INFO) << "input prepared, time spent: " << timer.lap() << "\n";
    int result = EXIT_FAILURE;
    if (par.createTargetTable) {
        result = createTargetTable(par, &reader, subMat);
    }
    else {
        result = createQueryTable(par, &reader, subMat);
    }

    delete subMat;
    reader.close();
    return result;
}

int createTargetTable(Parameters &par, DBReader<unsigned int> *reader, BaseMatrix *subMat) {
    Timer timer;
    size_t kmerCount = countKmer(reader, par.kmerSize);
    TargetTableEntry *targetTable = NULL;
    Debug(Debug::INFO) << "Number of sequences: " << reader->getSize() << "\n"
                       << "Number of all overall kmers: " << kmerCount << "\n"
                       << "Creating TargetTable. Requiring " << ((kmerCount + 1) * sizeof(TargetTableEntry)) / 1024 / 1024 << " MB of memory for it\n";
    size_t page_size = Util::getPageSize();
    targetTable = (TargetTableEntry *)mem_align(page_size, (kmerCount + 1) * sizeof(TargetTableEntry));
    if (madvise(targetTable, (kmerCount + 1) * sizeof(TargetTableEntry), MADV_HUGEPAGE | MADV_SEQUENTIAL) != 0) {
        Debug(Debug::WARNING) << "madvise returned an error\n";
    }
    Debug(Debug::INFO) << "Memory allocated \n"
                       << timer.lap() << "\n"
                       << "Extracting k-mers\n";

    int seqType = reader->getDbtype();
    size_t tableIndex = 0;
    Debug::Progress progress(reader->getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int)omp_get_thread_num();
#endif
        Indexer idx(subMat->alphabetSize - 1, par.kmerSize);
        Sequence s(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false);
        TargetTableEntry *localBuffer = new TargetTableEntry[BUFFERSIZE];
        size_t localTableIndex = 0;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < reader->getSize(); ++i) {
            progress.updateProgress();
            char *data = reader->getData(i, thread_idx);
            unsigned int seqLen = reader->getSeqLen(i);
            s.mapSequence(i, 0, data, seqLen);
            const int xIndex = s.subMat->aa2int[(int)'X'];
            while (s.hasNextKmer()) {
                const int *kmer = s.nextKmer();
                if (xCountInSequence(kmer, par.kmerSize, xIndex)) {
                    continue;
                }

                localBuffer[localTableIndex].sequenceID = i;
                localBuffer[localTableIndex].kmerAsLong = idx.int2index(kmer, 0, par.kmerSize);
                localBuffer[localTableIndex].sequenceLength = reader->getSeqLen(i);
                ++localTableIndex;
                if (localTableIndex >= BUFFERSIZE) {
                    size_t writeOffset = __sync_fetch_and_add(&tableIndex, localTableIndex);
                    memcpy(targetTable + writeOffset, localBuffer, sizeof(TargetTableEntry) * localTableIndex);
                    localTableIndex = 0;
                }
            }
        }
        if (localTableIndex > 0) {
            size_t writeOffset = __sync_fetch_and_add(&tableIndex, localTableIndex);
            memcpy(targetTable + writeOffset, localBuffer, sizeof(TargetTableEntry) * localTableIndex);
        }
        delete[] localBuffer;
    }

    Debug(Debug::INFO) << "k-mers: " << tableIndex << " time: " << timer.lap() << "\n";
    Debug(Debug::INFO) << "start sorting \n";
    omptl::sort(targetTable, targetTable + tableIndex, targetTableSort);
    Debug(Debug::INFO) << timer.lap() << "\n";
    writeTargetTables(targetTable, tableIndex, par.db2);
    Debug(Debug::INFO) << timer.lap() << "\n";
    free(targetTable);
    return EXIT_SUCCESS;
}

int createQueryTable(Parameters &par, DBReader<unsigned int> *reader, BaseMatrix *subMat) {
    Timer timer;
    size_t kmerCount = countKmer(reader, par.kmerSize);
    QueryTableEntry *queryTable = NULL;
    Debug(Debug::INFO) << "Number of sequences: " << reader->getSize() << "\n"
                       << "Number of all overall k-mers: " << kmerCount << "\n"
                       << "Creating QueryTable. Requiring " << ((kmerCount + 1) * sizeof(QueryTableEntry)) / 1024 / 1024 << " MB of memory for it\n";
    size_t page_size = Util::getPageSize();
    queryTable = (QueryTableEntry *)mem_align(page_size, (kmerCount + 1) * sizeof(QueryTableEntry));
    if (madvise(queryTable, (kmerCount + 1) * sizeof(QueryTableEntry), MADV_HUGEPAGE | MADV_SEQUENTIAL) != 0) {
        Debug(Debug::WARNING) << "madvise returned an error\n";
    }
    Debug(Debug::INFO) << "Memory allocated \n"
                       << timer.lap() << "\n"
                       << "Extracting k-mers\n";

    const int xIndex = subMat->aa2int[(int)'X'];

    size_t tableIndex = 0;
    int seqType = reader->getDbtype();
    Debug::Progress progress(reader->getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int)omp_get_thread_num();
#endif
        Indexer idx(subMat->alphabetSize - 1, par.kmerSize);
        Sequence s(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false);
        QueryTableEntry *localBuffer = new QueryTableEntry[BUFFERSIZE];
        size_t localTableIndex = 0;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < reader->getSize(); ++i) {
            progress.updateProgress();
            char *data = reader->getData(i, thread_idx);
            unsigned int seqLen = reader->getSeqLen(i);
            s.mapSequence(i, 0, data, seqLen);
            short kmerPosInSequence = 0;
            while (s.hasNextKmer()) {
                const int *kmer = s.nextKmer();
                if (xCountInSequence(kmer, par.kmerSize, xIndex)) {
                    continue;
                }

                localBuffer[localTableIndex].querySequenceId = i;
                localBuffer[localTableIndex].targetSequenceID = UINT_MAX;
                localBuffer[localTableIndex].Query.kmer = idx.int2index(kmer, 0, par.kmerSize);
                localBuffer[localTableIndex].Query.kmerPosInQuery = kmerPosInSequence;
                ++localTableIndex;
                ++kmerPosInSequence;

                if (localTableIndex >= BUFFERSIZE) {
                    size_t writeOffset = __sync_fetch_and_add(&tableIndex, localTableIndex);
                    memcpy(queryTable + writeOffset, localBuffer, sizeof(QueryTableEntry) * localTableIndex);
                    localTableIndex = 0;
                }
            }
        }
        if (localTableIndex > 0) {
            size_t writeOffset = __sync_fetch_and_add(&tableIndex, localTableIndex);
            memcpy(queryTable + writeOffset, localBuffer, sizeof(QueryTableEntry) * localTableIndex);
        }
        delete[] localBuffer;
    }

    Debug(Debug::INFO) << "\nkmers: " << tableIndex << " time: " << timer.lap() << "\n";
    Debug(Debug::INFO) << "start sorting \n";

    omptl::sort(queryTable, queryTable + tableIndex, queryTableSort);
    Debug(Debug::INFO) << timer.lap() << "\n";
    writeQueryTable(queryTable, tableIndex, par.db2);
    Debug(Debug::INFO) << timer.lap() << "\n";
    free(queryTable);
    return EXIT_SUCCESS;
}

size_t countKmer(DBReader<unsigned int> *reader, unsigned int kmerSize) {
    size_t kmerCount = 0;
    for (size_t i = 0; i < reader->getSize(); i++) {
        size_t currentSequenceLength = reader->getSeqLen(i);
        //number of ungapped k-mers per sequence = seq.length-k-mer.size+1
        kmerCount += currentSequenceLength >= kmerSize ? currentSequenceLength - kmerSize + 1 : 0;
    }
    return kmerCount;
}

int xCountInSequence(const int *kmer, size_t kmerSize, const int xIndex) {
    int xCount = 0;
    for (size_t pos = 0; pos < kmerSize; pos++) {
        xCount += (kmer[pos] == xIndex);
    }
    return xCount;
}

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

int targetTableSort(const TargetTableEntry &first, const TargetTableEntry &second) {
    if (first.kmerAsLong < second.kmerAsLong) {
        return true;
    }
    if (second.kmerAsLong < first.kmerAsLong) {
        return false;
    }
    if (first.sequenceLength > second.sequenceLength) {
        return true;
    }
    if (second.sequenceLength > first.sequenceLength) {
        return false;
    }
    if (first.sequenceID < second.sequenceID) {
        return true;
    }
    if (second.sequenceID < first.sequenceID) {  
        return false;
    }
    return false;
}

void writeTargetTables(TargetTableEntry *targetTable, size_t kmerCount, std::string blockID) {
    std::string kmerTableFileName = blockID + "_k-merTable";
    std::string idTablefileName = blockID + "_IDTable";
    Debug(Debug::INFO) << "Writing k-mer target table to file: " << kmerTableFileName << "\n";
    Debug(Debug::INFO) << "Writing target ID table to file:  " << idTablefileName << "\n";
    FILE *handleKmerTable = fopen(kmerTableFileName.c_str(), "wb");
    FILE *handleIDTable = fopen(idTablefileName.c_str(), "wb");
    TargetTableEntry *entryToWrite = targetTable;
    TargetTableEntry *posInTable = targetTable;
    size_t uniqueKmerCount = 0;
    size_t lastKmer = 0;
    Debug::Progress progress(kmerCount);
    
    for (size_t i = 0; i < kmerCount; ++i, ++posInTable) {
        progress.updateProgress();
        if (posInTable->kmerAsLong != entryToWrite->kmerAsLong) {
            writeKmerDiff(lastKmer, entryToWrite, handleKmerTable, handleIDTable);
            lastKmer = entryToWrite->kmerAsLong;
            entryToWrite = posInTable;
            ++uniqueKmerCount;
        }
    }
    //write last one
    writeKmerDiff(lastKmer, entryToWrite, handleKmerTable, handleIDTable);
    ++uniqueKmerCount;

    fclose(handleKmerTable);
    fclose(handleIDTable);
    Debug(Debug::INFO) << "Wrote " << uniqueKmerCount << " unique k-mers.\n";
    Debug(Debug::INFO) << "For "<< entryDiffLargerUShortMax  << " entries the difference between the previous were larger than maxshort.\n";
    Debug(Debug::INFO) << "Created " << diffLargerThenUShortMax << " extra unsigned short max entries to store the diff.\n";
}

void writeKmerDiff(size_t lastKmer, TargetTableEntry *entryToWrite, FILE *handleKmerTable, FILE *handleIDTable) {
    size_t kmerdiff = entryToWrite->kmerAsLong - lastKmer;
    bool first = true;
    unsigned short maxshort = USHRT_MAX;
    while (kmerdiff > maxshort) {
        if(first) {
            ++entryDiffLargerUShortMax;
            first = false;
        }
        ++diffLargerThenUShortMax;
        fwrite(&(maxshort), sizeof(unsigned short), 1, handleKmerTable);
        fwrite(&(entryToWrite->sequenceID), sizeof(unsigned int), 1, handleIDTable);
        kmerdiff -= maxshort;
    }
    fwrite((short *)&(kmerdiff), sizeof(unsigned short), 1, handleKmerTable);
    fwrite(&(entryToWrite->sequenceID), sizeof(unsigned int), 1, handleIDTable);
}

void writeQueryTable(QueryTableEntry *queryTable, size_t kmerCount, std::string queryID) {
    std::string fileName = queryID + "_queryTable";
    Debug(Debug::INFO) << "Writing query table to file: " << fileName << "\n";
    FILE *handleQueryTable = fopen(fileName.c_str(), "wb");
    fwrite(queryTable, sizeof(QueryTableEntry), kmerCount, handleQueryTable);
    fclose(handleQueryTable);
}
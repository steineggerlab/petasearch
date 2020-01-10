#include "LocalParameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "NucleotideMatrix.h"
#include "MathUtil.h"
#include "QueryTableEntry.h"
#include "TargetTableEntry.h"
#include "ExtendedSubstitutionMatrix.h"
#include "KmerGenerator.h"

#include "omptl/omptl_algorithm"
#include <sys/mman.h>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

#define BUFFERSIZE Util::getPageSize()

void writeQueryTable(const std::vector<QueryTableEntry>& queryTable, const std::string& queryID);
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
        subMat = new NucleotideMatrix(par.seedScoringMatrixFile.nucleotides, 1.0, 0.0);
    } else {
        subMat = new SubstitutionMatrix(par.seedScoringMatrixFile.aminoacids, 8.0, -0.2f);
    }
    Debug(Debug::INFO) << "input prepared, time spent: " << timer.lap() << "\n";
    int result;
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
    Debug(Debug::INFO) << "Number of sequences: " << reader->getSize() << "\n"
                       << "Number of all overall k-mers: " << kmerCount << "\n"
                       << "Creating QueryTable. Requiring " << ((kmerCount + 1) * sizeof(QueryTableEntry)) / 1024 / 1024 << " MB of memory for it\n";

    float similarKmerFactor = 1.5 ;
    size_t tableCapacity = (size_t) (similarKmerFactor * (kmerCount + 1));
    std::vector<QueryTableEntry> queryTable;
    queryTable.reserve(tableCapacity);

    const int xIndex = subMat->aa2int[(int)'X'];

    int seqType = reader->getDbtype();
    Debug::Progress progress(reader->getSize());

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
        for (size_t i = 0; i < reader->getSize(); ++i) {
            progress.updateProgress();
            char *data = reader->getData(i, thread_idx);
            unsigned int seqLen = reader->getSeqLen(i);
            sequence.mapSequence(i, 0, data, seqLen);

            while (sequence.hasNextKmer()) {
                const int *kmer = sequence.nextKmer();
                if (xCountInSequence(kmer, par.kmerSize, xIndex)) {
                    continue;
                }

                if (par.exactKmerMatching) {
                    QueryTableEntry entry;
                    entry.querySequenceId = i;
                    entry.targetSequenceID = UINT_MAX;
                    entry.Query.kmer = idx.int2index(kmer, 0, par.kmerSize);
                    entry.Query.kmerPosInQuery = sequence.getCurrentPosition();
                    localTable.emplace_back(entry);
                } else {
                    std::pair<size_t*, size_t> similarKmerList = kmerGenerator.generateKmerList(kmer);
//                    size_t *end = similarKmerList.first+similarKmerList.second;
//                    size_t  * originalkmer = std::find(similarKmerList.first, end ,idx.int2index(kmer, 0, par.kmerSize));
//                    if(originalkmer == end && similarKmerList.second != 0){
//                        Debug(Debug::ERROR) << "original k-mer not in List\n";
//                        exit(-1);
//                    }
                    for (size_t j = 0; j < similarKmerList.second; ++j) {
                        QueryTableEntry entry;
                        entry.querySequenceId = i;
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
        <<"\nRequired Memory: " << queryTable.size()  * sizeof(QueryTableEntry) / 1024 / 1024 << " MB \n"
        << " time: " << timer.lap() << "\n";
    Debug(Debug::INFO) << "start sorting \n";

    omptl::sort(queryTable.begin(), queryTable.end(), queryTableSort);
    Debug(Debug::INFO) << timer.lap() << "\n";
    writeQueryTable(queryTable, par.db2);
    Debug(Debug::INFO) << timer.lap() << "\n";
    return EXIT_SUCCESS;
}

short getKmerThreshold() {
    //TODO add function
    return 100;
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
    std::string kmerTableFileName = blockID;
    std::string idTableFileName = blockID + "_ids";
    Debug(Debug::INFO) << "Writing k-mer target table to file: " << kmerTableFileName << "\n";
    Debug(Debug::INFO) << "Writing target ID table to file:  " << idTableFileName << "\n";
    FILE *handleKmerTable = fopen(kmerTableFileName.c_str(), "wb");
    FILE *handleIDTable = fopen(idTableFileName.c_str(), "wb");
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
    Debug(Debug::INFO) << "For "<< entryDiffLargerUShortMax  << " entries the difference between the previous were larger than max short.\n";
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

void writeQueryTable(const std::vector<QueryTableEntry>& queryTable, const std::string& queryID) {
    Debug(Debug::INFO) << "Writing query table to file: " << queryID << "\n";
    FILE *handleQueryTable = fopen(queryID.c_str(), "wb");
    // C++11 ??
    fwrite(queryTable.data(), sizeof(QueryTableEntry), queryTable.size(), handleQueryTable);
    fclose(handleQueryTable);
}
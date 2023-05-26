#include "LocalParameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "SRADBReader.h"
#include "NucleotideMatrix.h"
#include "QueryTableEntry.h"
#include "TargetTableEntry.h"
#include "ExtendedSubstitutionMatrix.h"
#include "KmerGenerator.h"
#include "BitManipulateMacros.h"
#include "FastSort.h"

#include <sys/mman.h>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

#define KMER_BUFSIZ 500000000
#define ID_BUFSIZ 250000000

void writeTargetTables(TargetTableEntry *targetTable, size_t kmerCount, const std::string &blockID);

int queryTableSort(const QueryTableEntry &first, const QueryTableEntry &second);

int targetTableSort(const TargetTableEntry &first, const TargetTableEntry &second);

void writeKmerDiff(size_t lastKmer, TargetTableEntry *entryToWrite, FILE *handleKmerTable, FILE *handleIDTable,
                   uint16_t *kmerBuf, unsigned int *IDBuf);

static inline void writeKmer(uint16_t *buffer, FILE *handleKmerTable, uint16_t *toWrite, size_t size);

static inline void writeID(unsigned int *buffer, FILE *handleIDTable, unsigned int toWrite);

static inline void flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable);

static inline void flushIDBuf(unsigned int *buffer, FILE *handleIDTable);

static unsigned int kmerBufIdx = 0;
static unsigned int IDBufIdx = 0;

int createkmertable(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, true, 0, 0);
    Timer timer;

    SRADBReader reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    BaseMatrix *subMat;
    int seqType = reader.getDbtype();
    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.seedScoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    } else {
        subMat = new SubstitutionMatrix(par.seedScoringMatrixFile.values.aminoacid().c_str(), 8.0, -0.2f);
    }
    Debug(Debug::INFO) << "input prepared, time spent: " << timer.lap() << "\n";
    size_t kmerCount = 0;
    const unsigned int kmerSize = par.kmerSize;
#pragma omp parallel for default(none) shared(par, reader) firstprivate(kmerSize) reduction(+:kmerCount)
    for (size_t i = 0; i < reader.getSize(); ++i) {
        size_t currentSequenceLength = reader.getSeqLen(i);
        //number of ungapped k-mers per sequence = seq.length-k-mer.size+1
        kmerCount += currentSequenceLength >= (size_t) kmerSize ? currentSequenceLength - kmerSize + 1 : 0;
    }
    TargetTableEntry *targetTable = NULL;
    Debug(Debug::INFO) << "Number of sequences: " << reader.getSize() << "\n"
                       << "Number of all overall kmers: " << kmerCount << "\n"
                       << "Target table requires "
                       << ((kmerCount + 1) * sizeof(TargetTableEntry)) / 1024 / 1024 << " MB memory\n";

    size_t targetTableSize = std::min(Util::getTotalSystemMemory() - 32UL * 1024UL * 1024UL * 1024UL, (kmerCount + 1) * sizeof(TargetTableEntry));
    // TODO: check if overflow with target maximum index
    targetTable = (TargetTableEntry *) calloc(targetTableSize, sizeof(char));
    // (kmerCount + 1), sizeof(TargetTableEntry));

    if (targetTable == NULL) {
        Debug(Debug::ERROR) << "Could not allocate memory for target table\n";
        EXIT(EXIT_FAILURE);
    }

    const size_t pageSize = Util::getPageSize();
    const size_t threadBufferSize = 16 * pageSize;

    const int xIndex = subMat->aa2num[(int) 'X'];

    size_t tableIndex = 0;
    Debug::Progress progress(reader.getSize());
#pragma omp parallel default(none) shared(par, subMat, seqType, reader, tableIndex, targetTable, pageSize, threadBufferSize) firstprivate(xIndex)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        Indexer idx(subMat->alphabetSize - 1, par.kmerSize);
        Sequence s(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false);
        TargetTableEntry *localBuffer = (TargetTableEntry *) mem_align(pageSize, threadBufferSize * sizeof(TargetTableEntry));
        size_t localTableIndex = 0;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < reader.getSize(); ++i) {
//            progress.updateProgress();
            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, thread_idx);
            unsigned int seqLen = reader.getSeqLen(i);
            s.mapSequence(i, key, data, seqLen);
            while (s.hasNextKmer()) {
                const unsigned char *kmer = s.nextKmer();
                int xCount = 0;
                for (size_t pos = 0; pos < (unsigned) par.kmerSize; pos++) {
                    xCount += (kmer[pos] == xIndex);
                }
                if (xCount) {
                    continue;
                }

                localBuffer[localTableIndex].sequenceID = s.getId(); // for debug purposes: s.getDbKey();
                localBuffer[localTableIndex].kmerAsLong = idx.int2index(kmer, 0, par.kmerSize);
                localBuffer[localTableIndex].sequenceLength = s.L;
                ++localTableIndex;
                if (localTableIndex >= threadBufferSize) {
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
        free(localBuffer);
    }

    Debug(Debug::INFO) << "k-mers: " << tableIndex << " time: " << timer.lap() << "\n";
    SORT_PARALLEL(targetTable, targetTable + tableIndex, targetTableSort);
    Debug(Debug::INFO) << "Sorting time: " << timer.lap() << "\n";
    writeTargetTables(targetTable, tableIndex, par.db2);
    Debug(Debug::INFO) << "Writing time: " << timer.lap() << "\n";
    free(targetTable);

    delete subMat;
    subMat = nullptr;
    reader.close();
    return EXIT_SUCCESS;
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

void writeTargetTables(TargetTableEntry *targetTable, size_t kmerCount, const std::string &blockID) {
    const std::string &kmerTableFileName = blockID;
    std::string idTableFileName = blockID + "_ids";
    Debug(Debug::INFO) << "Writing k-mer target table to file: " << kmerTableFileName << "\n";
    Debug(Debug::INFO) << "Writing target ID table to file:  " << idTableFileName << "\n";
    FILE *handleKmerTable = fopen(kmerTableFileName.c_str(), "wb");
    FILE *handleIDTable = fopen(idTableFileName.c_str(), "wb");
    TargetTableEntry *entryToWrite = targetTable;
    TargetTableEntry *posInTable = targetTable;
    size_t uniqueKmerCount = 0;
    size_t lastKmer = 0;
//    Debug::Progress progress(kmerCount);

    uint16_t *kmerLocalBuf = (uint16_t *) malloc(sizeof(uint16_t) * KMER_BUFSIZ);
    unsigned int *IDLocalBuf = (unsigned int *) malloc(sizeof(unsigned int) * ID_BUFSIZ);
    for (size_t i = 0; i < kmerCount; ++i, ++posInTable) {
//        progress.updateProgress();
        if (posInTable->kmerAsLong != entryToWrite->kmerAsLong) {
            writeKmerDiff(lastKmer, entryToWrite, handleKmerTable, handleIDTable, kmerLocalBuf, IDLocalBuf);
            lastKmer = entryToWrite->kmerAsLong;
            entryToWrite = posInTable;
            ++uniqueKmerCount;
        }
    }
    // write last one
    writeKmerDiff(lastKmer, entryToWrite, handleKmerTable, handleIDTable, kmerLocalBuf, IDLocalBuf);
    ++uniqueKmerCount;

    flushKmerBuf(kmerLocalBuf, handleKmerTable);
    flushIDBuf(IDLocalBuf, handleIDTable);
    free(kmerLocalBuf);
    free(IDLocalBuf);
    fclose(handleKmerTable);
    fclose(handleIDTable);
    Debug(Debug::INFO) << "Wrote " << uniqueKmerCount << " unique k-mers\n";
//    Debug(Debug::INFO) << "For "<< entryDiffLargerUShortMax  << " entries the difference between the previous were larger than max short.\n";
//    Debug(Debug::INFO) << "Created " << diffLargerThenUShortMax << " extra unsigned short max entries to store the diff.\n";
}

static inline void flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable) {
    fwrite(buffer, sizeof(uint16_t), kmerBufIdx, handleKmerTable);
    kmerBufIdx = 0;
}

static inline void flushIDBuf(unsigned int *buffer, FILE *handleIDTable) {
    fwrite(buffer, sizeof(unsigned int), IDBufIdx, handleIDTable);
    IDBufIdx = 0;
}

static inline void writeKmer(uint16_t *buffer, FILE *handleKmerTable, uint16_t *toWrite, size_t size) {
    if (kmerBufIdx + size >= KMER_BUFSIZ) {
        flushKmerBuf(buffer, handleKmerTable);
    }
    memcpy(buffer + kmerBufIdx, toWrite, sizeof(uint16_t) * size);
    kmerBufIdx += size;
}

static inline void writeID(unsigned int *buffer, FILE *handleIDTable, unsigned int toWrite) {
    if (IDBufIdx == ID_BUFSIZ) {
        flushIDBuf(buffer, handleIDTable);
    }
    buffer[IDBufIdx] = toWrite;
    IDBufIdx++;
}

void writeKmerDiff(size_t lastKmer, TargetTableEntry *entryToWrite, FILE *handleKmerTable, FILE *handleIDTable,
                   uint16_t *kmerBuf, unsigned int *IDBuf) {
    uint64_t kmerdiff = entryToWrite->kmerAsLong - lastKmer;
    // Consecutively store 15 bits of information into a short, until kmer diff is all
    uint16_t buffer[5]; // 15*5 = 75 > 64
    buffer[4] = SET_END_FLAG(GET_15_BITS(kmerdiff));
    kmerdiff >>= 15U;
    int idx = 3;
    while (kmerdiff) {
        uint16_t toWrite = GET_15_BITS(kmerdiff);
        kmerdiff >>= 15U;
        buffer[idx] = toWrite;
        idx--;
    }
    writeKmer(kmerBuf, handleKmerTable, (buffer + idx + 1), (4 - idx));
    writeID(IDBuf, handleIDTable, entryToWrite->sequenceID);
}

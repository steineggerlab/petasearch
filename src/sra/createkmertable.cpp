#include "LocalParameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "NucleotideMatrix.h"
#include "MathUtil.h"
#include "QueryTableEntry.h"
#include "TargetTableEntry.h"
#include "ExtendedSubstitutionMatrix.h"
#include "KmerGenerator.h"
#include "BitManipulateMacros.h"

#include <sys/mman.h>
#include <algorithm>

#include "ips4o/ips4o.hpp"

#ifdef OPENMP
#include <omp.h>
#endif

void writeTargetTables(TargetTableEntry *targetTable, size_t kmerCount, std::string blockID);
int queryTableSort(const QueryTableEntry &first, const QueryTableEntry &second);
int targetTableSort(const TargetTableEntry &first, const TargetTableEntry &second);
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
    size_t kmerCount = 0;
#pragma omp parallel for reduction(+:kmerCount)
    for (size_t i = 0; i < reader.getSize(); ++i) {
        size_t currentSequenceLength = reader.getSeqLen(i);
        //number of ungapped k-mers per sequence = seq.length-k-mer.size+1
        kmerCount += currentSequenceLength >= (unsigned)par.kmerSize ? currentSequenceLength - par.kmerSize + 1 : 0;
    }
    TargetTableEntry *targetTable = NULL;
    Debug(Debug::INFO) << "Number of sequences: " << reader.getSize() << "\n"
                       << "Number of all overall kmers: " << kmerCount << "\n"
                       << "Creating TargetTable. Requiring " << ((kmerCount + 1) * sizeof(TargetTableEntry)) / 1024 / 1024 << " MB of memory for it\n";
    targetTable = (TargetTableEntry *)calloc((kmerCount + 1), sizeof(TargetTableEntry));
    Debug(Debug::INFO) << "Memory allocated \n"
                       << timer.lap() << "\n"
                       << "Extracting k-mers\n";

    const size_t pageSize = Util::getPageSize();
    const size_t threadBufferSize = 16 * pageSize;

    size_t tableIndex = 0;
    Debug::Progress progress(reader.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int)omp_get_thread_num();
#endif
        Indexer idx(subMat->alphabetSize - 1, par.kmerSize);
        Sequence s(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false);
        TargetTableEntry *localBuffer = (TargetTableEntry *)mem_align(pageSize, threadBufferSize * sizeof(TargetTableEntry));
        size_t localTableIndex = 0;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();
            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, thread_idx);
            unsigned int seqLen = reader.getSeqLen(i);
            s.mapSequence(i, 0, data, seqLen);
            const int xIndex = s.subMat->aa2num[(int)'X'];
            while (s.hasNextKmer()) {
                const unsigned char*kmer = s.nextKmer();
                int xCount = 0;
                for (size_t pos = 0; pos < (unsigned)par.kmerSize; pos++) {
                    xCount += (kmer[pos] == xIndex);
                }
                if (xCount) {
                    continue;
                }

                localBuffer[localTableIndex].sequenceID = key;
                localBuffer[localTableIndex].kmerAsLong = idx.int2index(kmer, 0, par.kmerSize);
                localBuffer[localTableIndex].sequenceLength = reader.getSeqLen(i);
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
    Debug(Debug::INFO) << "start sorting \n";
    ips4o::parallel::sort(targetTable, targetTable + tableIndex, targetTableSort);
    Debug(Debug::INFO) << timer.lap() << "\n";
    writeTargetTables(targetTable, tableIndex, par.db2);
    Debug(Debug::INFO) << timer.lap() << "\n";
    free(targetTable);


    delete subMat;
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

    // TODO: add a buffer array to save the overhead of writing
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
//    Debug(Debug::INFO) << "For "<< entryDiffLargerUShortMax  << " entries the difference between the previous were larger than max short.\n";
//    Debug(Debug::INFO) << "Created " << diffLargerThenUShortMax << " extra unsigned short max entries to store the diff.\n";
}

void writeKmerDiff(size_t lastKmer, TargetTableEntry *entryToWrite, FILE *handleKmerTable, FILE *handleIDTable) {
    uint64_t kmerdiff = entryToWrite->kmerAsLong - lastKmer;
    // Consecutively store 15 bits of information into a short, until kmer diff is all
    uint16_t buffer[5] = { 0 }; // 15*5 = 75 > 64
    buffer[4] = SET_END_FLAG(GET_15_BITS(kmerdiff));
    kmerdiff >>= 15U;
    int idx = 3;
    while (kmerdiff) {
        uint16_t toWrite = GET_15_BITS(kmerdiff);
        kmerdiff >>= 15U;
        buffer[idx] = toWrite;
        idx--;
    }
    int start = -1;
    for (int i = 0; i < 5; i++) {
        uint16_t bits = buffer[i];
        if (bits) {
            start = i;
            break;
        }
    }
    for (; start < 5; start++) {
        uint16_t bits = buffer[start];
        fwrite(&(bits), sizeof(uint16_t), 1, handleKmerTable);
        if (IS_LAST_15_BITS(bits)) {
            fwrite(&(entryToWrite->sequenceID), sizeof(unsigned int), 1, handleIDTable);
        }
    }
}

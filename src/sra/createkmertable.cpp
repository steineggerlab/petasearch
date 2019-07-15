#include "LocalParameters.h"
#include "Command.h"
#include "Debug.h"
#include "IndexReader.h"
#include "DBReader.h"
#include "NucleotideMatrix.h"
#include "fstream"
#include "MathUtil.h"
#include "QueryTableEntry.h"
#include "TargetTableEntry.h"
#include "omptl/omptl_algorithm"
#include <sys/mman.h>

#include <algorithm>
#include <cassert>

#ifdef OPENMP
#include <omp.h>
#endif

#define KMER_SIZE 9
#define SPACED_KMER true


long kmer2long(const int* index, size_t kmerSize, long **aminoAcidValueAtPosition );
void writeQueryTable(QueryTableEntry* queryTable, size_t kmerCount, std::string queryID);
void writeTargetTables(TargetTableEntry* targetTable, size_t kmerCount, std::string blockID);
int queryTableSort(const QueryTableEntry &first, const QueryTableEntry &second);
int targetTableSort(const TargetTableEntry &first, const TargetTableEntry &second);
size_t countKmer(DBReader<unsigned int> *reader, Parameters & par);
int xCountInSequence(const int* kmer, size_t kmerSize, const int xIndex);
int createQueryTable(Parameters& par, DBReader<unsigned int> *reader,  BaseMatrix * subMat, long** aminoAcidValueAtPosition);
int createTargetTable(Parameters& par, DBReader<unsigned int> *reader,  BaseMatrix * subMat, long** aminoAcidValueAtPosition);

int isQuery(){
    return false;
}


int createkmertable(int argc, const char ** argv, const Command& command){
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.kmerSize = KMER_SIZE;
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, 2, true);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    // if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        // reader.readMmapedDataInMemory();
    // }

    BaseMatrix * subMat;
    int seqType = reader.getDbtype();
    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.c_str(), 1.0, 0.0);
    } else {
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);
    }

    long** aminoAcidValueAtPosition;
    aminoAcidValueAtPosition = new long *[par.kmerSize];
    for(int i = 0; i < par.kmerSize; i++){
        long posBaseValue = MathUtil::ipow<long>(subMat->alphabetSize, i);
        aminoAcidValueAtPosition[i] = new long[subMat->alphabetSize];
        for (int aa = 0; aa < subMat->alphabetSize; aa++){
            aminoAcidValueAtPosition[i][aa] = aa * posBaseValue;
        }  
    }

    int result = EXIT_FAILURE;
    if(par.createTargetTable){
        result = createTargetTable(par, &reader, subMat, aminoAcidValueAtPosition);
    }else{
        result = createQueryTable(par, &reader, subMat, aminoAcidValueAtPosition);
    }

    delete subMat;
    for(int i = 0; i < par.kmerSize; ++i){
        delete [] aminoAcidValueAtPosition[i];
    }
    delete [] aminoAcidValueAtPosition;
    reader.close();
    return result;
}

int createTargetTable(Parameters& par, DBReader<unsigned int> *reader,  BaseMatrix * subMat, long** aminoAcidValueAtPosition){
    Timer timer;
    size_t kmerCount = countKmer(reader,par);
    TargetTableEntry* targetTable = NULL;
    Debug(Debug::INFO) << "Number of sequences: " << reader->getSize() << "\n"
                       << "Number of all overall kmers: " << kmerCount << "\n"
                       << "Creating TargetTable. Requiring " << ((kmerCount+1)*sizeof(TargetTableEntry))/1024/1024 << " MB of memory for it\n";
    size_t page_size = Util::getPageSize();
    targetTable = (TargetTableEntry*) mem_align(page_size, (kmerCount + 1) * sizeof(TargetTableEntry));
    if (madvise(targetTable,  (kmerCount + 1) * sizeof(TargetTableEntry), MADV_HUGEPAGE|MADV_SEQUENTIAL) != 0){
        Debug(Debug::WARNING) << "madvise returned an error\n";
    }
    Debug(Debug::INFO) << "Memory allocated \n"
                       << timer.lap() << "\n"
                       << "Extracting k-mers\n";

    int seqType = reader->getDbtype();
    size_t tableIndex = 0;    
    #pragma omp parallel 
    {   
        unsigned int thread_idx = 0;
        #ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
        #endif
        Indexer idx(subMat->alphabetSize-1, par.kmerSize);
        Sequence s(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false);

        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < reader->getSize(); ++i) {
            char *data = reader->getData(i, thread_idx);
            s.mapSequence(i, 0, data);
            short kmerPosInSequence = 0;
            const int xIndex = s.subMat->aa2int[(int) 'X'];
            while (s.hasNextKmer()) {
                const int *kmer = s.nextKmer();
                ++kmerPosInSequence;
               if(xCountInSequence(kmer,par.kmerSize, xIndex)){
                    continue;
                }
                size_t localTableIndex = __sync_fetch_and_add(&tableIndex, 1);
                targetTable[localTableIndex].sequenceID = i;
                // targetTable[localTableIndex].kmerAsLong = kmer2long(kmer,par.kmerSize,aminoAcidValueAtPosition);
                targetTable[localTableIndex].kmerAsLong = idx.int2index(kmer, 0, par.kmerSize);
                targetTable[localTableIndex].sequenceLength = reader->getSeqLens(i) - 2;
            }
        }
    }

    Debug(Debug::INFO) << "kmers: " << tableIndex << " time: "<< timer.lap() << "\n";
    Debug(Debug::INFO) << "start sorting \n";
    omptl::sort(targetTable, targetTable+kmerCount, targetTableSort);
    Debug(Debug::INFO) << timer.lap() << "\n";
    writeTargetTables(targetTable,kmerCount, par.db2);
    Debug(Debug::INFO) << timer.lap() << "\n";
    free(targetTable);
    return EXIT_SUCCESS;
}

int createQueryTable(Parameters& par, DBReader<unsigned int> *reader,  BaseMatrix * subMat, long** aminoAcidValueAtPosition){
    Timer timer;
    size_t kmerCount = countKmer(reader,par);
    QueryTableEntry* queryTable = NULL;
    Debug(Debug::INFO) << "Number of sequences: " << reader->getSize() << "\n"
                       << "Number of all overall kmers: " << kmerCount << "\n"
                       << "Creating QueryTable. Requiring " << ((kmerCount+1)*sizeof(QueryTableEntry))/1024/1024 << " MB of memory for it\n";
    size_t page_size = Util::getPageSize();
    queryTable = (QueryTableEntry*) mem_align(page_size, (kmerCount + 1) * sizeof(QueryTableEntry));
    if (madvise(queryTable,  (kmerCount + 1) * sizeof(QueryTableEntry), MADV_HUGEPAGE|MADV_SEQUENTIAL) != 0){
        Debug(Debug::WARNING) << "madvise returned an error\n";
    }
    Debug(Debug::INFO) << "Memory allocated \n"
                       << timer.lap() << "\n"
                       << "Extracting k-mers\n";

    const int xIndex = subMat->aa2int[(int) 'X'];

    size_t tableIndex = 0;    
    int seqType = reader->getDbtype();
    #pragma omp parallel 
    {   
        unsigned int thread_idx = 0;
        #ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
        #endif
        Indexer idx(subMat->alphabetSize-1, par.kmerSize);
        Sequence s(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false);

        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < reader->getSize(); ++i) {
            char *data = reader->getData(i, thread_idx);
            s.mapSequence(i, 0, data);
            short kmerPosInSequence = 0;
            // bool hadKmer = false;
            while (s.hasNextKmer()) {
                // hadKmer = true;
                const int *kmer = s.nextKmer();
                ++kmerPosInSequence;
                if(xCountInSequence(kmer,par.kmerSize, xIndex)){
                    continue;
                }

                size_t localTableIndex = __sync_fetch_and_add(&tableIndex, 1);
                queryTable[localTableIndex].querySequenceId = i;
                queryTable[localTableIndex].targetSequenceID = 0;
                // queryTable[localTableIndex].Query.kmer = kmer2long(kmer, par.kmerSize, aminoAcidValueAtPosition);
                queryTable[localTableIndex].Query.kmer = idx.int2index(kmer, 0, par.kmerSize);
                queryTable[localTableIndex].Query.kmerPosInQuery = kmerPosInSequence;

                // Debug(Debug::INFO) << kmer2long(kmer, par.kmerSize, aminoAcidValueAtPosition) << "\n";
                // idx.printKmer(kmer, par.kmerSize, subMat->int2aa);
                // Debug(Debug::INFO) << "\n" << idx.int2index(kmer, 0, par.kmerSize) << "\n";
            }
            // if (hadKmer)
            // break;
        }
    }

    Debug(Debug::INFO) << "kmers: " << tableIndex << " time: "<< timer.lap() << "\n";
    Debug(Debug::INFO) << "start sorting \n";

    omptl::sort(queryTable,queryTable+tableIndex, queryTableSort);
    Debug(Debug::INFO) << timer.lap() << "\n";
    writeQueryTable(queryTable,tableIndex,par.db2);
    Debug(Debug::INFO) << timer.lap() << "\n";
    free(queryTable);
    return EXIT_SUCCESS;
}

size_t countKmer(DBReader<unsigned int> *reader, Parameters & par){
    size_t kmerCount = 0;
    for (size_t i = 0; i < reader->getSize(); i++) {
        size_t currentSequenceLength = reader->getSeqLens(i) - 2;
        //number of ungapped k-mers per sequence = seq.length-k-mer.size+1
        kmerCount += currentSequenceLength >= (unsigned int )par.kmerSize ? currentSequenceLength - par.kmerSize + 1 : 0;
    }
    return kmerCount; 
}

int xCountInSequence(const int* kmer, size_t kmerSize, const int xIndex){
    int xCount = 0;
    for (size_t pos = 0; pos < kmerSize; pos++) {
        xCount += (kmer[pos] == xIndex);
    }
    return xCount;
}

long kmer2long(const int* index, size_t kmerSize, long** aminoAcidValueAtPosition){
    long kmerAsLong = 0;
    for(size_t i = 0; i < kmerSize; ++i){
        kmerAsLong += aminoAcidValueAtPosition[i][index[i]];
    }
    return kmerAsLong;
}

int queryTableSort(const QueryTableEntry &first, const QueryTableEntry &second){
    if (first.Query.kmer < second.Query.kmer)
        return true;
    if (second.Query.kmer < first.Query.kmer)
        return false;
    if (first.querySequenceId > second.querySequenceId)
        return true;
    if (second.querySequenceId > first.querySequenceId)
        return false;
    if (first.Query.kmerPosInQuery < second.Query.kmerPosInQuery)
        return true;
    if (second.Query.kmerPosInQuery < first.Query.kmerPosInQuery)
        return false;
    return false;
}

int targetTableSort(const TargetTableEntry &first, const TargetTableEntry &second){
    if (first.kmerAsLong < second.kmerAsLong)
        return true;
    if (second.kmerAsLong < first.kmerAsLong)
        return false;
    if (first.sequenceLength > second.sequenceLength)
        return true;
    if (second.sequenceLength > first.sequenceLength)
        return false;
    if (first.sequenceID < second.sequenceID)
        return true;
    if (second.sequenceID < first.sequenceID)
        return false;
    return false;
}

void writeTargetTables(TargetTableEntry* targetTable, size_t kmerCount, std::string blockID){
    std::string kmerTableFileName = blockID + "_k-merTable";
    std::string idTablefileName = blockID + "_IDTable";
    Debug(Debug::INFO) << "Writing k-mer target table to file: " << kmerTableFileName <<"\n";
    Debug(Debug::INFO) << "Writing target ID table to file:  " << idTablefileName <<"\n";
    FILE* handleKmerTable = fopen(kmerTableFileName.c_str(), "wb");
    FILE* handleIDTable = fopen(idTablefileName.c_str(),"wb");
    TargetTableEntry* entryToWrite = targetTable;
    TargetTableEntry* posInTable = targetTable;
    size_t uniqueKmerCount = 0;
    for(size_t i = 0; i < kmerCount; ++i, ++posInTable){
        if(posInTable->kmerAsLong == entryToWrite->kmerAsLong){
            if(posInTable->sequenceLength > entryToWrite->sequenceLength){
                entryToWrite = posInTable;
            }
        }else{
            fwrite(&(entryToWrite->kmerAsLong),sizeof(long),1,handleKmerTable);
            fwrite(&(entryToWrite->sequenceID),sizeof(unsigned int),1,handleIDTable);
            entryToWrite = posInTable;
            ++uniqueKmerCount;
        }
        //TODO Test
        // if(posInTable->kmerAsLong != entryToWrite->kmerAsLong){
        //     fwrite(&(entryToWrite->kmerAsLong),sizeof(entryToWrite->kmerAsLong),1,handleKmerTable);
        //     fwrite(&(entryToWrite->sequenceID),sizeof(entryToWrite->sequenceID),1,handleIDTable);
        //     entryToWrite=posInTable;
        //     ++uniqueKmerCount;
        // }
    }
    //write last one
    fwrite(&(entryToWrite->kmerAsLong),sizeof(long),1,handleKmerTable);
    fwrite(&(entryToWrite->sequenceID),sizeof(unsigned int),1,handleIDTable);
    ++uniqueKmerCount;

    fclose(handleKmerTable);
    fclose(handleIDTable);
    Debug(Debug::INFO)<<"Wrote "<< uniqueKmerCount << " unique k-mers.\n";
}

void writeQueryTable(QueryTableEntry* queryTable, size_t kmerCount, std::string queryID){
    std::string fileName = "queryTable_"+queryID;
    Debug(Debug::INFO) << "Writing query table to file: " << fileName <<"\n";
    FILE* handleQueryTable = fopen(fileName.c_str(),"wb");
    fwrite(queryTable,sizeof(QueryTableEntry),kmerCount,handleQueryTable);
    fclose(handleQueryTable);
}
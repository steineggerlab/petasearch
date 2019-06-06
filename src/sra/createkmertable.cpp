#include "Parameters.h"
#include "Command.h"
#include "Debug.h"
#include "IndexReader.h"
#include "DBReader.h"
#include "NucleotideMatrix.h"
#include "fstream"
#include "MathUtil.h"
#include "QueryTableEntry.h"
#include "TargetTableEntry.h"
#include <algorithm>

#define KMER_SIZE 9
#define SPACED_KMER true


long kmer2long(const int* index, size_t kmerSize, size_t alphabetSize);
void writeQueryTable(QueryTableEntry* querryTable, size_t kmerCount, std::string querryID);
void writeTargetTables(TargetTableEntry* targetTable, size_t kmerCount, std::string blockID);
int querryTableSort(const QueryTableEntry &entryOne, const QueryTableEntry &entry2);
int targetTableSort(const TargetTableEntry &entryOne, const TargetTableEntry &entryTwo);


int isQuery(){
    return false;
}


int createkmertable(int argc, const char ** argv, const Command& command){
    Parameters& par = Parameters::getInstance();
    par.kmerSize = KMER_SIZE;
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, 2, false);
    int indexSrcType = IndexReader::SEQUENCES;

    IndexReader reader(par.db1, par.threads, indexSrcType, 0);
    int seqType = reader.sequenceReader->getDbtype();
    BaseMatrix * subMat;

    size_t isNucl=Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES);
    if (isNucl) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.c_str(), 1.0, 0.0);
    } else {
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);
    }

    size_t maxLen = 0;
    size_t kmerCount=0;
    for(size_t i = 0; i < reader.sequenceReader->getSize(); i++){
        size_t currentSequenceLenth = reader.sequenceReader->getSeqLens(i) - 2;
        //number of ungapped k-mers per sequence = seq.length-k-mer.size+1
        kmerCount+=reader.sequenceReader->getSeqLens(i) - 2 - par.kmerSize+1;
        maxLen = std::max(maxLen, currentSequenceLenth);
    }

    if(isQuery()){
        Debug(Debug::INFO) << "creating QueryTable with " <<kmerCount <<" entries. Requiring " 
                           << (kmerCount*sizeof(QueryTableEntry))/1024/1024 << " MB of memory for it\n";
        QueryTableEntry* querryTable = (QueryTableEntry*) calloc(kmerCount,sizeof(QueryTableEntry));
        Debug(Debug::INFO)<< "Memory allocated \n"
        <<"Number of sequences: " << reader.sequenceReader->getSize()<<"\n";
        #pragma omp parallel 
        {
            Indexer idx(subMat->alphabetSize-1, par.kmerSize);
            Sequence s(maxLen, seqType, subMat,
                        par.kmerSize, par.spacedKmer, false);
            #pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < reader.sequenceReader->getSize(); ++i) {
                char *data = reader.sequenceReader->getData(i, 0);
                s.mapSequence(i, 0, data);
                size_t queryIndex = 0;
                const int xIndex = s.subMat->aa2int[(int) 'X'];
                short kmerPosInSequence = 0;
                while (s.hasNextKmer()) {
                    const int *kmer = s.nextKmer();
                    kmerPosInSequence++;
                    int xCount = 0;
                    for (int pos = 0; pos < par.kmerSize; pos++) {
                        xCount += (kmer[pos] == xIndex);
                    }
                    if (xCount > 0) {
                        continue;
                    }
                
                    queryIndex = __sync_fetch_and_add(&queryIndex, 1);
                    querryTable[queryIndex].querySequenceId = s.getId();;
                    querryTable[queryIndex].Query.kmer = kmer2long(kmer,par.kmerSize,subMat->alphabetSize-1);
                    querryTable[queryIndex].Query.kmerPosInQuerry = kmerPosInSequence;
                }
            }
        }
        Debug(Debug::INFO) << "start sorting \n";
        std::stable_sort(querryTable,querryTable+kmerCount,querryTableSort);
        writeQueryTable(querryTable,kmerCount,"Test");
        free(querryTable);
            
    }else{
        Timer timer;
        Debug(Debug::INFO) << "Number of all overall kmers: " << kmerCount <<" Requiring " 
                           << (kmerCount*sizeof(TargetTableEntry))/1024/1024 << " MB of memory while creating the target table\n";
        TargetTableEntry* targetTable = (TargetTableEntry*) calloc(kmerCount,sizeof(TargetTableEntry));
        size_t targetIndex = 0;
        Debug(Debug::INFO) << "Memory allocated \n"
                           << "Number of sequences: " << reader.sequenceReader->getSize()<<"\n";
        #pragma omp parallel 
        {
            Indexer idx(subMat->alphabetSize-1, par.kmerSize);
            Sequence s(maxLen, seqType, subMat,
                        par.kmerSize, par.spacedKmer, false);

            #pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < reader.sequenceReader->getSize(); ++i) {
                char *data = reader.sequenceReader->getData(i, 0);
                s.mapSequence(i, 0, data);
                
                const int xIndex = s.subMat->aa2int[(int) 'X'];
                while (s.hasNextKmer()) {
                    const int *kmer = s.nextKmer();
                    int xCount = 0;
                    for (int pos = 0; pos < par.kmerSize; pos++) {
                        xCount += (kmer[pos] == xIndex);
                    }
                    if (xCount > 0) {
                        continue;
                    }

                    targetIndex = __sync_fetch_and_add(&targetIndex, 1);
                    targetTable[targetIndex].sequenceID = s.getId();
                    targetTable[targetIndex].kmerAsLong = kmer2long(kmer,par.kmerSize,subMat->alphabetSize-1);
                    targetTable[targetIndex].sequenceLength = reader.sequenceReader->getSeqLens(i) - 2;
                }
            }
        }

        Debug(Debug::INFO) << timer.lap() << "\n";
        Debug(Debug::INFO) << "start sorting \n";
        std::stable_sort(targetTable,targetTable+kmerCount,targetTableSort);
        Debug(Debug::INFO) << timer.lap() << "\n";
        writeTargetTables(targetTable,kmerCount,"test");
        Debug(Debug::INFO) << timer.lap() << "\n";
        Debug(Debug::INFO)<<(targetTable)->kmerAsLong << "\n";
        free(targetTable);
       
    }
    return EXIT_SUCCESS;
}


int querryTableSort(const QueryTableEntry &entryOne, const QueryTableEntry &entryTwo){
    return entryOne.Query.kmer < entryTwo.Query.kmer;
}

int targetTableSort(const TargetTableEntry &entryOne, const TargetTableEntry &entryTwo){
    if(entryOne.kmerAsLong != entryTwo.kmerAsLong){
        return entryOne.kmerAsLong < entryTwo.kmerAsLong;
    }
    return entryOne.sequenceLength < entryTwo.sequenceLength;
    
}


void writeTargetTables(TargetTableEntry* targetTable, size_t kmerCount, std::string blockID){
    std::string kmerTableFileName = blockID + "_k-merTable";
    std::string idTablefileName = blockID + "_IDTable";
    Debug(Debug::INFO) << "Writing k-mer target table to file: " << kmerTableFileName <<"\n";
    Debug(Debug::INFO) << "Writing target ID table to file:  " << idTablefileName <<"\n";
    FILE* handleKmerTable = fopen(kmerTableFileName.c_str(), "ab+");
    FILE* handleIDTable = fopen(idTablefileName.c_str(),"ab+");
    TargetTableEntry* entryToWrite = targetTable;
    TargetTableEntry* posInTable = targetTable;
    size_t uniqueKmerCount = 0;
    for(size_t i = 0; i < kmerCount; ++i, ++posInTable){
        if(posInTable->kmerAsLong == entryToWrite->kmerAsLong){
            if(posInTable->sequenceLength > entryToWrite->sequenceLength){
                entryToWrite=posInTable;
            }
        }else{
            fwrite(&(entryToWrite->kmerAsLong),sizeof(entryToWrite->kmerAsLong),1,handleKmerTable);
            fwrite(&(entryToWrite->sequenceID),sizeof(entryToWrite->sequenceID),1,handleIDTable);
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
    fwrite(&(entryToWrite->kmerAsLong),sizeof(entryToWrite->kmerAsLong),1,handleKmerTable);
    fwrite(&(entryToWrite->sequenceID),sizeof(entryToWrite->sequenceID),1,handleIDTable);
    ++uniqueKmerCount;

    fclose(handleKmerTable);
    fclose(handleIDTable);
    Debug(Debug::INFO)<<"Wrote "<< uniqueKmerCount << " unique k-mers.\n";
}

void writeQueryTable(QueryTableEntry* querryTable, size_t kmerCount, std::string querryID){
    std::string fileName = "querryTable_"+querryID;
    Debug(Debug::INFO) << "Writing query table to file: " << fileName <<"\n";
    FILE* handleQueryTable = fopen(fileName.c_str(),"ab+");
    fwrite(querryTable,sizeof(querryTable),kmerCount,handleQueryTable);
    fclose(handleQueryTable);
}



long kmer2long(const int* index, size_t kmerSize, size_t alphabetSize){
    long kmerAsLong = 0;
    for(size_t i = 0; i < kmerSize; ++i){
        kmerAsLong += index[i]*MathUtil::ipow<long>(alphabetSize, i);
    }
    return kmerAsLong;

}





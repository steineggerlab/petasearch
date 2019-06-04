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

#define KMER_SIZE 5
#define SPACED_KMER true


long index2long(size_t index[], size_t kmerSize, size_t alphabetSize);
void writeKmerTableUsingOfStream(char * filename, unsigned int kmerCountTable[], BaseMatrix *subMat);
void writeKmerTableFWrite(char * filename, unsigned char kmerCountTable[], BaseMatrix *subMat );
void writeQueryTable(QueryTableEntry* querryTable, size_t kmerCount, std::string querryID);

int isQuery(){
    return true;
}
int querryTableSort(QueryTableEntry entryOne, QueryTableEntry entry2);
int targetTableSort(TargetTableEntry entryOne, TargetTableEntry entryTwo);


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
        size_t currentSequenceLenth = reader.sequenceReader->getSeqLens(i);
        kmerCount+=reader.sequenceReader->getSeqLens(i)-par.kmerSize+1;
        maxLen = std::max(maxLen, currentSequenceLenth);
    }

    if(isQuery()){
        Debug(Debug::INFO) << "creating QueryTable with " <<kmerCount <<" entries. Requiring " << (kmerCount*sizeof(QueryTableEntry))/1024/1024 << " MB of memory for it\n";
        QueryTableEntry* querryTable = (QueryTableEntry*) calloc(kmerCount,sizeof(QueryTableEntry));
        QueryTableEntry* posInQueryTable  = querryTable;
        Debug(Debug::INFO)<< "Memory allocated \n"
        <<"Number of sequences: " << reader.sequenceReader->getSize()<<"\n";
        #pragma omp parallel 
        {
            Indexer idx(subMat->alphabetSize-1, par.kmerSize);
            Sequence s(maxLen, seqType, subMat,
                        par.kmerSize, par.spacedKmer, false);
            #pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < reader.sequenceReader->getSize(); ++i) {
                QueryTableEntry * localpos;
                char *data = reader.sequenceReader->getData(i, 0);
                s.mapSequence(i, 0, data);
                
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
                
                    #pragma omp atomic capture
                    localpos = posInQueryTable++;                    
            
                    localpos->querySequenceId = i;
                    localpos->Query.kmer = index2long(idx.workspace,KMER_SIZE,subMat->alphabetSize-1);
                    localpos->Query.kmerPosInQuerry = kmerPosInSequence;
                }
            }
        }
        Debug(Debug::INFO) << "start sorting \n";
        std::sort(querryTable,querryTable+kmerCount,querryTableSort);
        writeQueryTable(querryTable,kmerCount,"Test");
            
    }else{
        Debug(Debug::INFO) << "creating TargetTable. Requiring " << (kmerCount*sizeof(TargetTableEntry))/1024/1024 << " MB\n";
        TargetTableEntry* targetTable = (TargetTableEntry*) calloc(kmerCount,sizeof(TargetTableEntry));
        TargetTableEntry* posInTargetTable = targetTable;
        Debug(Debug::INFO)<< "Memory allocated \n";
        #pragma omp parallel 
        {
            Indexer idx(subMat->alphabetSize-1, par.kmerSize);
            Sequence s(maxLen, seqType, subMat,
                        par.kmerSize, par.spacedKmer, false);

            #pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < reader.sequenceReader->getSize(); ++i) {
                TargetTableEntry * localpos;
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

                    #pragma omp atomic capture
                    localpos = posInTargetTable++;

                    localpos->sequenceID = i;
                    localpos->kmerAsLong = index2long(idx.workspace,KMER_SIZE,subMat->alphabetSize-1);
                    localpos->sequenceLength = reader.sequenceReader->getSeqLens(i);
                }
            }
        }
        std::sort(targetTable,targetTable+kmerCount,targetTableSort);
        //TODO write to target-kmer table and target-id table.
       
    }
    return 0;
}


int querryTableSort(QueryTableEntry entryOne, QueryTableEntry entryTwo){
    return entryOne.Query.kmer > entryTwo.Query.kmer;
}

int targetTableSort(TargetTableEntry entryOne, TargetTableEntry entryTwo){
    if(entryOne.kmerAsLong != entryTwo.kmerAsLong){
        return entryOne.kmerAsLong > entryTwo.kmerAsLong;
    }
    return entryOne.sequenceLength > entryTwo.sequenceLength;
    
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
    for(size_t i = 0; i < kmerCount; ++i, ++posInTable){
        if(posInTable->kmerAsLong == entryToWrite->kmerAsLong){
            if(posInTable->sequenceLength > entryToWrite->sequenceLength){
                entryToWrite=posInTable;
            }
        }else{
            fwrite(&(entryToWrite->kmerAsLong),sizeof(entryToWrite->kmerAsLong),1,handleKmerTable);
            fwrite(&(entryToWrite->sequenceID),sizeof(entryToWrite->sequenceID),1,handleIDTable);
            entryToWrite = posInTable;
        }
        //TODO Test
        // if(posInTable->kmerAsLong != entryToWrite->kmerAsLong){
        //     fwrite(&(entryToWrite->kmerAsLong),sizeof(entryToWrite->kmerAsLong),1,handleKmerTable);
        //     fwrite(&(entryToWrite->sequenceID),sizeof(entryToWrite->sequenceID),1,handleIDTable);
        //     entryToWrite=posInTable;
        // }
    }
    fclose(handleKmerTable);
    fclose(handleKmerTable);
}

void writeQueryTable(QueryTableEntry* querryTable, size_t kmerCount, std::string querryID){
    std::string fileName = "querryTable_"+querryID;
    Debug(Debug::INFO) << "Writing query table to file: " << fileName <<"\n";
    FILE* handleQueryTable = fopen(fileName.c_str(),"ab+");
    fwrite(querryTable,sizeof(querryTable),kmerCount,handleQueryTable);
    fclose(handleQueryTable);
}

// int createkmertable(int argc, const char **argv, const Command& command){
//     Parameters& par = Parameters::getInstance();
//     par.kmerSize = KMER_SIZE;
//     par.spacedKmer = false;
//     par.parseParameters(argc, argv, command, 2, false);
//     int indexSrcType = IndexReader::SEQUENCES;
//     Debug(Debug::INFO) << par.db1 << "\n";
    
//     IndexReader reader(par.db1, par.threads, indexSrcType, 0);
//     int seqType = reader.sequenceReader->getDbtype();
//     BaseMatrix * subMat;

//     size_t isNucl=Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES);
//     if (isNucl) {
//         subMat = new NucleotideMatrix(par.scoringMatrixFile.c_str(), 1.0, 0.0);
//     } else {
//         subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);
//     }

//     size_t maxLen = 0;
//     for(size_t i = 0; i < reader.sequenceReader->getSize(); i++){
//         maxLen = std::max(maxLen, reader.sequenceReader->getSeqLens(i));
//     }

//     size_t idxSize = MathUtil::ipow<size_t>(subMat->alphabetSize-1, par.kmerSize);
//     Debug(Debug::INFO) << "Index Size: " << idxSize << "\n";
//     struct timeval startTime;
//     struct timeval endTime; 
//     Debug(Debug::INFO) << "Starting zeroing memory" << "\n";
//     gettimeofday(&startTime, NULL);
//     //unsigned char* kmerCountTable=new unsigned char[idxSize];
//     //memset((long*)kmerCountTable,  0l, sizeof(unsigned char)*idxSize);
//     unsigned char * kmerCountTable = (unsigned char *) calloc(idxSize,sizeof(unsigned char));

//     gettimeofday(&endTime, NULL);
//     double timediff = (endTime.tv_sec - startTime.tv_sec) + 1e-6 * (endTime.tv_usec - startTime.tv_usec);
//     Debug(Debug::INFO) << "memory zerod. Required time: " << timediff << " seconds \n";
// #pragma omp parallel
//     {
//         Indexer idx(subMat->alphabetSize-1, par.kmerSize);
//         Sequence s(maxLen, seqType, subMat,
//                           par.kmerSize, par.spacedKmer, false);

// #pragma omp for schedule(dynamic, 1)
//         for (size_t i = 0; i < reader.sequenceReader->getSize(); ++i) {
//             char *data = reader.sequenceReader->getData(i, 0);
//             s.mapSequence(i, 0, data);
//             const int xIndex = s.subMat->aa2int[(int) 'X'];
//             while (s.hasNextKmer()) {
//                 const int *kmer = s.nextKmer();
//                 int xCount = 0;
//                 for (int pos = 0; pos < par.kmerSize; pos++) {
//                     xCount += (kmer[pos] == xIndex);
//                 }
//                 if (xCount > 0) {
//                     continue;
//                 }
//                 size_t kmerIdx = (isNucl) ? Indexer::computeKmerIdx(kmer, par.kmerSize) : idx.int2index(kmer, 0, par.kmerSize);
//                 //no need for syncronized, in all cases the index get increased at least by 1 which is sufficent
//                 // kmerCountTable[kmerIdx]++;
//         #pragma omp atomic write
// 		kmerCountTable[kmerIdx] = 1;
//                 //__sync_fetch_and_add(&kmerCountTable[kmerIdx], 1);
//             }
//         }
//     }

// 	gettimeofday(&startTime, NULL);
//     writeKmerTableFWrite((char*)par.db2.c_str(),kmerCountTable,subMat);
//     // writeKmerTableUsingOfStream(KMERTABLEFILE,kmerCountTable,subMat);
//     gettimeofday(&endTime, NULL);
//     timediff = (endTime.tv_sec - startTime.tv_sec) + 1e-6 * (endTime.tv_usec - startTime.tv_usec);
//     Debug(Debug::INFO)<<"Time required to write k-mer table to file (wall clock writing time): " << timediff << " s\n";

//     delete [] kmerCountTable;

//     return EXIT_SUCCESS;
// }

long index2long(size_t index[], size_t kmerSize, size_t alphabetSize){
    long kmerAsLong = 0;
    for(size_t i = 0; i < kmerSize; ++i){
        kmerAsLong += index[i]*MathUtil::ipow<long>(alphabetSize, i);
    }
    return kmerAsLong;

}

// void writeKmerTableFWrite(char* filename, unsigned char kmerCountTable[], BaseMatrix *subMat ){
//     FILE* handle = fopen(filename, "ab+");
//     Indexer idx(subMat->alphabetSize-1, KMER_SIZE);
//     size_t idxSize = MathUtil::ipow<size_t>(subMat->alphabetSize-1, KMER_SIZE);
//     long sum =0;
//     size_t count = 0;
//     for(size_t i = 0; i < idxSize; ++i){
//         int kmerCount = kmerCountTable[i];
//         if(kmerCount){
//             idx.index2int(idx.workspace, i, KMER_SIZE);
//             long kmerAsLong= index2long(idx.workspace,KMER_SIZE,subMat->alphabetSize-1);
//             fwrite(&kmerAsLong, sizeof(kmerAsLong),1,handle);
//             sum += kmerAsLong;
//             ++count;
//         }
//     }

//     Debug(Debug::INFO) << "sum: " << sum <<"\n" << "count: "<< count << "\n" << "last written: " 
//         << index2long(idx.workspace,KMER_SIZE,subMat->alphabetSize-1) << "\n";
//     fclose(handle);

// }

// void writeKmerTableUsingOfStream(char * filename, unsigned int kmerCountTable[], BaseMatrix *subMat ){
//     Indexer idx(subMat->alphabetSize-1, KMER_SIZE);
//     size_t idxSize = MathUtil::ipow<size_t>(subMat->alphabetSize-1, KMER_SIZE);
//     std::ofstream tableFile (filename, std::ios::out | std::ios::binary);
//     long sum =0;
//     size_t count = 0;
//     for(size_t i = 0; i < idxSize; ++i){
//         int kmerCount = kmerCountTable[i];
//         if(kmerCount){
//             idx.index2int(idx.workspace, i, KMER_SIZE);
//             long kmerAslong = index2long(idx.workspace,KMER_SIZE,subMat->alphabetSize-1);
//             tableFile.write((char*) &kmerAslong,sizeof(kmerAslong));
//             // std::cout<< kmerAslong<<std::endl;
//             sum += kmerAslong;
//             ++count;
//         }
//     }
    
//     Debug(Debug::INFO) << "sum: " << sum <<"\n" << "count: "<< count << "\n" << "last written: " 
//         << index2long(idx.workspace,KMER_SIZE,subMat->alphabetSize-1) << "\n";
//     tableFile.close();
// }




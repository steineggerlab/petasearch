#include "Parameters.h"
#include "Command.h"
#include "Debug.h"
#include "IndexReader.h"
#include "DBReader.h"
#include "NucleotideMatrix.h"
#include "fstream"
#include "MathUtil.h"

#define KMER_SIZE 5
#define SPACED_KMER true
// #define KMERTABLEFILE "temp/kmertableDecoded"

long index2long(size_t index[], size_t kmerSize, size_t alphabetSize);
void writeKmerTableUsingOfStream(char * filename, unsigned int kmerCountTable[], BaseMatrix *subMat);
void writeKmerTableFWrite(char * filename, unsigned char kmerCountTable[], BaseMatrix *subMat );

int createkmertable(int argc, const char **argv, const Command& command){
    Parameters& par = Parameters::getInstance();
    par.kmerSize = KMER_SIZE;
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, 2, false);
    int indexSrcType = IndexReader::SEQUENCES;
    Debug(Debug::INFO) << par.db1 << "\n";
    
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
    for(size_t i = 0; i < reader.sequenceReader->getSize(); i++){
        maxLen = std::max(maxLen, reader.sequenceReader->getSeqLens(i));
    }

    size_t idxSize = MathUtil::ipow<size_t>(subMat->alphabetSize-1, par.kmerSize);
    Debug(Debug::INFO) << "Index Size: " << idxSize << "\n";
    struct timeval startTime;
    struct timeval endTime; 
    Debug(Debug::INFO) << "Starting zeroing memory" << "\n";
    gettimeofday(&startTime, NULL);
    //unsigned char* kmerCountTable=new unsigned char[idxSize];
    //memset((long*)kmerCountTable,  0l, sizeof(unsigned char)*idxSize);
    unsigned char * kmerCountTable = (unsigned char *) calloc(idxSize,sizeof(unsigned char));

    gettimeofday(&endTime, NULL);
    double timediff = (endTime.tv_sec - startTime.tv_sec) + 1e-6 * (endTime.tv_usec - startTime.tv_usec);
    Debug(Debug::INFO) << "memory zerod. Required time: " << timediff << " seconds \n";
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
                size_t kmerIdx = (isNucl) ? Indexer::computeKmerIdx(kmer, par.kmerSize) : idx.int2index(kmer, 0, par.kmerSize);
                //no need for syncronized, in all cases the index get increased at least by 1 which is sufficent
                // kmerCountTable[kmerIdx]++;
        #pragma omp atomic write
		kmerCountTable[kmerIdx] = 1;
                //__sync_fetch_and_add(&kmerCountTable[kmerIdx], 1);
            }
        }
    }

	gettimeofday(&startTime, NULL);
    writeKmerTableFWrite((char*)par.db2.c_str(),kmerCountTable,subMat);
    // writeKmerTableUsingOfStream(KMERTABLEFILE,kmerCountTable,subMat);
    gettimeofday(&endTime, NULL);
    timediff = (endTime.tv_sec - startTime.tv_sec) + 1e-6 * (endTime.tv_usec - startTime.tv_usec);
    Debug(Debug::INFO)<<"Time required to write k-mer table to file (wall clock writing time): " << timediff << " s\n";

    delete [] kmerCountTable;

    return EXIT_SUCCESS;
}

long index2long(size_t index[], size_t kmerSize, size_t alphabetSize){
    long kmerAsLong = 0;
    for(size_t i = 0; i < kmerSize; ++i){
        kmerAsLong += index[i]*MathUtil::ipow<long>(alphabetSize, i);
    }
    return kmerAsLong;

}

void writeKmerTableFWrite(char* filename, unsigned char kmerCountTable[], BaseMatrix *subMat ){
    FILE* handle = fopen(filename, "ab+");
    Indexer idx(subMat->alphabetSize-1, KMER_SIZE);
    size_t idxSize = MathUtil::ipow<size_t>(subMat->alphabetSize-1, KMER_SIZE);
    long sum =0;
    size_t count = 0;
    for(size_t i = 0; i < idxSize; ++i){
        int kmerCount = kmerCountTable[i];
        if(kmerCount){
            idx.index2int(idx.workspace, i, KMER_SIZE);
            long kmerAsLong= index2long(idx.workspace,KMER_SIZE,subMat->alphabetSize-1);
            fwrite(&kmerAsLong, sizeof(kmerAsLong),1,handle);
            sum += kmerAsLong;
            ++count;
        }
    }

    Debug(Debug::INFO) << "sum: " << sum <<"\n" << "count: "<< count << "\n" << "last written: " 
        << index2long(idx.workspace,KMER_SIZE,subMat->alphabetSize-1) << "\n";
    fclose(handle);

}

void writeKmerTableUsingOfStream(char * filename, unsigned int kmerCountTable[], BaseMatrix *subMat ){
    Indexer idx(subMat->alphabetSize-1, KMER_SIZE);
    size_t idxSize = MathUtil::ipow<size_t>(subMat->alphabetSize-1, KMER_SIZE);
    std::ofstream tableFile (filename, std::ios::out | std::ios::binary);
    long sum =0;
    size_t count = 0;
    for(size_t i = 0; i < idxSize; ++i){
        int kmerCount = kmerCountTable[i];
        if(kmerCount){
            idx.index2int(idx.workspace, i, KMER_SIZE);
            long kmerAslong = index2long(idx.workspace,KMER_SIZE,subMat->alphabetSize-1);
            tableFile.write((char*) &kmerAslong,sizeof(kmerAslong));
            // std::cout<< kmerAslong<<std::endl;
            sum += kmerAslong;
            ++count;
        }
    }
    
    Debug(Debug::INFO) << "sum: " << sum <<"\n" << "count: "<< count << "\n" << "last written: " 
        << index2long(idx.workspace,KMER_SIZE,subMat->alphabetSize-1) << "\n";
    tableFile.close();
}




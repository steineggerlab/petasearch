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
void writeKmerTableFWrite(char * filename, unsigned int kmerCountTable[], BaseMatrix *subMat );

int createkmertable(int argc, const char **argv, const Command& command){
    char *KMERTABLEFILE = "tmp/kmerTableDecodedAsLong";    
    Parameters& par = Parameters::getInstance();
    par.kmerSize = KMER_SIZE;
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, 1, false);
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
    unsigned int * kmerCountTable=new unsigned int[idxSize];
    memset(kmerCountTable, 0, sizeof(unsigned int)*idxSize);

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
                kmerCountTable[kmerIdx]++;
            }
        }
    }

    // writeKmerTableFWrite(KMERTABLEFILE,kmerCountTable,subMat);
    writeKmerTableUsingOfStream(KMERTABLEFILE,kmerCountTable,subMat);



    // Indexer idx(subMat->alphabetSize-1, par.kmerSize);
    // size_t count = 0;
    // char * filename = KMERTABLEFILE;
    // FILE* handle = fopen(filename, "ab+");
    // int fd=fileno(handle);
    // struct stat fileStat;
    // fstat(fd, &fileStat);
    // char *buffer = new char[KMER_SIZE];
    // 
    // for(size_t i = 0; i < idxSize; ++i){
    //     int kmerCount = kmerCountTable[i];
    //     if(kmerCount){
    //         ++count;
    //         idx.index2int(idx.workspace, i, par.kmerSize);
    //         long kmerAsLong= index2long(idx.workspace,KMER_SIZE,19);
    //         fwrite(&kmerAsLong, sizeof(kmerAsLong),1,handle);
    //         std::cout<< kmerAsLong<<std::endl;
    //         sum+=kmerAsLong;
    //     }

    // }
    //  idx.index2int(idx.workspace, idxSize-1, par.kmerSize);
    //  std::cout<<index2long(idx.workspace,KMER_SIZE,19)<<std::endl;
    //  std::cout<<"sum: "<< sum<<std::endl;

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

void writeKmerTableFWrite(char * filename, unsigned int kmerCountTable[], BaseMatrix *subMat ){
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




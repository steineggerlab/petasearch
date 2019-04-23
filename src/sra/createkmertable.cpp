#include "Parameters.h"
#include "Command.h"
#include "Debug.h"
#include "IndexReader.h"
#include "DBReader.h"
#include "NucleotideMatrix.h"

#define KMER_SIZE 5
#define SPACED_KMER true
#define KMERTABLEFILE "temp/kmertableDecoded"


int createkmertable(int argc, const char **argv, const Command& command){
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
                //no need for syncronized, in all cases the index get increased at least by 1 which is enough
                kmerCountTable[kmerIdx]++;
                // __sync_fetch_and_add(&(kmerCountTable[kmerIdx]), 1);
            }
        }
    }
    Indexer idx(subMat->alphabetSize-1, par.kmerSize);
    size_t count=0;
    for(size_t i=0; i<idxSize;++i){
        int kmerCount = kmerCountTable[i];
        if(kmerCount){
            ++count;
            idx.index2int(idx.workspace, i, par.kmerSize);
            for(int k = 0; k < par.kmerSize; k++){
                std::cout << subMat->int2aa[idx.workspace[k]];
            }
            std::cout<<"\n"<<"";
        }

    }
    std::cout<<count<<std::endl;
    delete [] kmerCountTable;

    return EXIT_SUCCESS;
}


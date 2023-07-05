#ifndef FIXED_KMERGENERATOR_H 
#define FIXED_KMERGENERATOR_H 
#include <string>
#include <vector>
#include <Indexer.h>
#include <ScoreMatrix.h>

class FixedKmerGenerator {
    public: 
        FixedKmerGenerator(size_t kmerSize, size_t alphabetSize, short threshold, unsigned int maxKmers);
        ~FixedKmerGenerator();

        // calculates the kmer list
        std::pair<size_t *, size_t> generateKmerList(const unsigned char* intSeq, bool addIdentity = false);

        // kmer splitting stragety (3,2)
        // fill up the divide step and calls init_result_list
        void setDivideStrategy(ScoreMatrix* three, ScoreMatrix* two);

        // kmer splitting stragety (1)
        // fill up the divide step and calls init_result_list
        void setDivideStrategy(ScoreMatrix** one);

        void setThreshold(short threshold);
    private:
        // size of kmer
        size_t kmerSize;

        // min score
        short threshold;

        // maximum number of kmers to be generated
        unsigned int maxKmers;

        // partition steps of the kmer size in (2,3)
        size_t divideStepCount;
        // divider of the steps (2,3)
        unsigned int* divideStep;
        size_t* stepMultiplicator;
        Indexer indexer;
        ScoreMatrix** matrixLookup;
        short* outputScoreArray;
        size_t* outputIndexArray;
        std::pair<short *, int>* scoreArrays;
        std::pair<unsigned int*, int>* indexArrays;

        // init the output vectors for the kmer calculation
        void initDataStructure();
};

#endif

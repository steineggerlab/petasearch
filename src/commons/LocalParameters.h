#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

class LocalParameters : public Parameters {
public:
    static void initInstance() {
        new LocalParameters;
    }

    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initInstance();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }

    std::vector<MMseqsParameter*> createkmertable;
    std::vector<MMseqsParameter*> compare2kmertables;
    std::vector<MMseqsParameter*> computeAlignment;
    ///example

    std::vector<MMseqsParameter*> predictexonsworkflow;

    PARAMETER(PARAM_REVERSE_FRAGMENTS)
    int reverseFragments;
    //example end

    PARAMETER(PARAM_CREATE_TARGET_TABLE)
    int createTargetTable;

    PARAMETER(PARAM_REQ_KMER_MATCHES)
    unsigned int requieredKmerMatches;

private:
    LocalParameters() : 
        Parameters(),

        PARAM_REVERSE_FRAGMENTS(PARAM_REVERSE_FRAGMENTS_ID,"--reverse-fragments", "Reverse AA Fragments", "reverse AA fragments to compute under null [0,1]", typeid(int), (void *) &reverseFragments, "^[0-1]{1}$"),
        PARAM_CREATE_TARGET_TABLE(PARAM_CREATE_TARGET_TABLE_ID,"--createTargetTable","creating target table ?","create target table (1) or query table (0) (default: 1) [0,1]",typeid(int), (void *) &createTargetTable,"^[0-1]{1}$"),
        PARAM_REQ_KMER_MATCHES(PARAM_REQ_KMER_MATCHES_ID,"--req-kmer-matches","required k-mer matches per query and target sequence pair","aount of required k-mer matches per query/target pair to increase specifity of matches [0-4]",typeid(int), (void*) &requieredKmerMatches,"^[0-4]{1}$")
        {
            predictexonsworkflow.push_back(&PARAM_REVERSE_FRAGMENTS);

            createkmertable.push_back(&PARAM_K);
            createkmertable.push_back(&PARAM_SPACED_KMER_MODE);
            createkmertable.push_back(&PARAM_THREADS);
            createkmertable.push_back(&PARAM_V);
            createkmertable.push_back(&PARAM_CREATE_TARGET_TABLE);
            createkmertable.push_back(&PARAM_MAX_SEQ_LEN);

            compare2kmertables.push_back(&PARAM_THREADS);
            compare2kmertables.push_back(&PARAM_SPACED_KMER_MODE);
            compare2kmertables.push_back(&PARAM_REQ_KMER_MATCHES);
            compare2kmertables.push_back(&PARAM_COMPRESSED);

            computeAlignment.push_back(&PARAM_SPACED_KMER_MODE);
            computeAlignment.push_back(&PARAM_THREADS);
            computeAlignment.push_back(&PARAM_COMPRESSED);
            computeAlignment.push_back(&PARAM_NO_COMP_BIAS_CORR);
            computeAlignment.push_back(&PARAM_K);
            computeAlignment.push_back(&PARAM_GAP_OPEN);
            computeAlignment.push_back(&PARAM_GAP_EXTEND);
            computeAlignment.push_back(&PARAM_E);
            computeAlignment.push_back(&PARAM_MAX_SEQ_LEN);


            // default value 0 means no reverse of AA fragments
            reverseFragments = 0;
            //default value 1 means to create the target table
            createTargetTable = 1;

            requieredKmerMatches = 2;
        }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);

    

};

#endif

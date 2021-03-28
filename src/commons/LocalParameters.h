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
    std::vector<MMseqsParameter*> computeAlignments;
    std::vector<MMseqsParameter*> petasearchworkflow;
    std::vector<MMseqsParameter*> easypetasearchworkflow;


    PARAMETER(PARAM_REQ_KMER_MATCHES)
    unsigned int requiredKmerMatches;

private:
    LocalParameters() : Parameters(),
        PARAM_REQ_KMER_MATCHES(PARAM_REQ_KMER_MATCHES_ID, "--req-kmer-matches", "required k-mer matches per query and target sequence pair", "amount of required k-mer matches per query/target pair to increase specificity of matches [0-4]", typeid(int), (void*) &requiredKmerMatches, "^[0-4]{1}$")
        {
            createkmertable.push_back(&PARAM_SEED_SUB_MAT);

            createkmertable.push_back(&PARAM_K);
            createkmertable.push_back(&PARAM_SPACED_KMER_MODE);
            createkmertable.push_back(&PARAM_MAX_SEQ_LEN);
            createkmertable.push_back(&PARAM_THREADS);
            createkmertable.push_back(&PARAM_V);

            compare2kmertables.push_back(&PARAM_EXACT_KMER_MATCHING);
            compare2kmertables.push_back(&PARAM_SEED_SUB_MAT);
            compare2kmertables.push_back(&PARAM_K);
            compare2kmertables.push_back(&PARAM_K_SCORE);
            compare2kmertables.push_back(&PARAM_SPACED_KMER_MODE);
            compare2kmertables.push_back(&PARAM_MAX_SEQ_LEN);
            compare2kmertables.push_back(&PARAM_REQ_KMER_MATCHES);
            compare2kmertables.push_back(&PARAM_COMPRESSED);
            compare2kmertables.push_back(&PARAM_THREADS);
            compare2kmertables.push_back(&PARAM_V);

            computeAlignments.push_back(&PARAM_SUB_MAT);
            computeAlignments.push_back(&PARAM_K);
            computeAlignments.push_back(&PARAM_SPACED_KMER_MODE);
            computeAlignments.push_back(&PARAM_NO_COMP_BIAS_CORR);
            computeAlignments.push_back(&PARAM_GAP_OPEN);
            computeAlignments.push_back(&PARAM_GAP_EXTEND);
            computeAlignments.push_back(&PARAM_E);
            computeAlignments.push_back(&PARAM_MAX_SEQ_LEN);
            computeAlignments.push_back(&PARAM_COMPRESSED);
            computeAlignments.push_back(&PARAM_THREADS);
            computeAlignments.push_back(&PARAM_V);

            petasearchworkflow = combineList(createkmertable, compare2kmertables);
            petasearchworkflow = combineList(petasearchworkflow, computeAlignments);
            petasearchworkflow = combineList(petasearchworkflow, swapresult);
            petasearchworkflow = combineList(petasearchworkflow, convertalignments);

            easypetasearchworkflow = combineList(petasearchworkflow,convertalignments);



            requiredKmerMatches = 2;

            kmerSize = 9;
            kmerScore = 225;
        }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};

#endif

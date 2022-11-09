#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

class LocalParameters : public Parameters {
public:
    static const int DBTYPE_SRA_DB = 233;
    static void initInstance() {
        new LocalParameters;
    }

    static LocalParameters &getLocalInstance() {
        if (instance == NULL) {
            initInstance();
        }
        return static_cast<LocalParameters &>(LocalParameters::getInstance());
    }

    std::vector<MMseqsParameter *> createkmertable;
    std::vector<MMseqsParameter *> comparekmertables;
    std::vector<MMseqsParameter *> computeAlignments;
    std::vector<MMseqsParameter *> convert2sradb;
    std::vector<MMseqsParameter *> petasearchworkflow;
    std::vector<MMseqsParameter *> easypetasearchworkflow;
    std::vector<MMseqsParameter *> convertsraalignments;
    std::vector<MMseqsParameter *> readandprint;

    PARAMETER(PARAM_REQ_KMER_MATCHES)
    unsigned int requiredKmerMatches;

    PARAMETER(PARAM_X_DROP)
    int8_t xdrop;

    PARAMETER(PARAM_RANGE_MIN)
    uintptr_t rangeMin;

    PARAMETER(PARAM_RANGE_MAX)
    uintptr_t rangeMax;

private:
    LocalParameters() : Parameters(),
                        PARAM_REQ_KMER_MATCHES(PARAM_REQ_KMER_MATCHES_ID,
                                               "--req-kmer-matches",
                                               "required k-mer matches per query and target sequence pair",
                                               "amount of required k-mer matches per query/target pair to increase specificity of matches [0-4]",
                                               typeid(int),
                                               (void *) &requiredKmerMatches,
                                               "^[0-4]{1}$"),
                        PARAM_X_DROP(PARAM_X_DROP_ID,
                                     "--xdrop",
                                     "the xdrop value for block aligner",
                                     "the xdrop value for block aligner [0-100]",
                                     typeid(int),
                                     (void *) &xdrop,
                                     "^[0-9]+$"),
                        PARAM_RANGE_MIN(PARAM_RANGE_MIN_ID,
                                        "--range-min",
                                        "the minimum size of block-aligner's matrix",
                                        "the minimum size of block-aligner' matrix, must be of power of 2 [32-1024]",
                                        typeid(int),
                                        (void *) &rangeMin,
                                        "^[0-9]+$"),
                        PARAM_RANGE_MAX(PARAM_RANGE_MAX_ID,
                                        "--range-max",
                                        "the maximum size of block-aligner's matrix",
                                        "the maximum size of block-aligner' matrix, must be of power of 2 [32-1024]",
                                        typeid(int),
                                        (void *) &rangeMax,
                                        "^[0-9]+$") {
        createkmertable.push_back(&PARAM_SEED_SUB_MAT);
        createkmertable.push_back(&PARAM_K);
        createkmertable.push_back(&PARAM_SPACED_KMER_MODE);
        createkmertable.push_back(&PARAM_MAX_SEQ_LEN);
        createkmertable.push_back(&PARAM_THREADS);
        createkmertable.push_back(&PARAM_V);

        comparekmertables.push_back(&PARAM_EXACT_KMER_MATCHING);
        comparekmertables.push_back(&PARAM_SEED_SUB_MAT);
        comparekmertables.push_back(&PARAM_K);
        comparekmertables.push_back(&PARAM_K_SCORE);
        comparekmertables.push_back(&PARAM_SPACED_KMER_MODE);
        comparekmertables.push_back(&PARAM_MAX_SEQ_LEN);
        comparekmertables.push_back(&PARAM_REQ_KMER_MATCHES);
        comparekmertables.push_back(&PARAM_COMPRESSED);
        comparekmertables.push_back(&PARAM_THREADS);
        comparekmertables.push_back(&PARAM_V);

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
        computeAlignments.push_back(&PARAM_RANGE_MIN);
        computeAlignments.push_back(&PARAM_RANGE_MAX);
        computeAlignments.push_back(&PARAM_X_DROP);

        convertsraalignments.push_back(&PARAM_FORMAT_MODE);
        convertsraalignments.push_back(&PARAM_FORMAT_OUTPUT);
        convertsraalignments.push_back(&PARAM_THREADS);
        convertsraalignments = combineList(convertsraalignments, Parameters::convertalignments);

        petasearchworkflow = combineList(createkmertable, comparekmertables);
        petasearchworkflow = combineList(petasearchworkflow, computeAlignments);
        petasearchworkflow = combineList(petasearchworkflow, swapresult);
        petasearchworkflow = combineList(petasearchworkflow, convertalignments);
        petasearchworkflow = combineList(petasearchworkflow, convertsraalignments);

        easypetasearchworkflow = combineList(petasearchworkflow, convertalignments);

        requiredKmerMatches = 2;
        xdrop = 10;
        rangeMin = 32;
        rangeMax = 32;

        kmerSize = 9;
        kmerScore = 225;
    }

    LocalParameters(LocalParameters const &);

    ~LocalParameters() {};

    void operator=(LocalParameters const &);
};

#endif

#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include "Parameters.h"

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
    std::vector<MMseqsParameter *> blockalign;
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

    PARAMETER(PARAM_MAX_KMER_PER_POS)
    int maxKmerPerPos;

private:
    LocalParameters() : Parameters(),
        PARAM_REQ_KMER_MATCHES(
            PARAM_REQ_KMER_MATCHES_ID,
            "--req-kmer-matches",
            "required k-mer matches per query and target sequence pair",
            "amount of required k-mer matches per query/target pair to increase specificity of matches [0-4]",
            typeid(int),
            (void *) &requiredKmerMatches,
            "^[0-4]{1}$"),
        PARAM_X_DROP(
            PARAM_X_DROP_ID,
            "--xdrop",
            "block-aligner x-drop",
            "block-aligner x-drop [0-100]",
            typeid(int),
            (void *) &xdrop,
            "^[0-9]+$"),
        PARAM_RANGE_MIN(
            PARAM_RANGE_MIN_ID,
            "--range-min",
            "block-aligner matrix min size",
            "Minimum size of block-aligner' matrix, must be of power of 2 [32-1024]",
            typeid(int),
            (void *) &rangeMin,
            "^[0-9]+$"),
        PARAM_RANGE_MAX(
            PARAM_RANGE_MAX_ID,
            "--range-max",
            "block-aligner matrix max size",
            "Maximum size of block-aligner' matrix, must be of power of 2 [32-1024]",
            typeid(int),
            (void *) &rangeMax,
            "^[0-9]+$"),
        PARAM_MAX_KMER_PER_POS(
            PARAM_MAX_KMER_PER_POS_ID,
            "--max-kmer-per-pos",
            "Maximum k-mers per position",
            "Maximum k-mers per position [>=1]",
            typeid(int),
            (void *) &maxKmerPerPos,
            "^[0-9]+$")
    {
        createkmertable.push_back(&PARAM_SEED_SUB_MAT);
        createkmertable.push_back(&PARAM_K);
        createkmertable.push_back(&PARAM_SPACED_KMER_MODE);
        createkmertable.push_back(&PARAM_SPACED_KMER_PATTERN);
        createkmertable.push_back(&PARAM_MAX_SEQ_LEN);
        createkmertable.push_back(&PARAM_THREADS);
        createkmertable.push_back(&PARAM_V);

        comparekmertables.push_back(&PARAM_EXACT_KMER_MATCHING);
        comparekmertables.push_back(&PARAM_SEED_SUB_MAT);
        comparekmertables.push_back(&PARAM_K);
        comparekmertables.push_back(&PARAM_K_SCORE);
        comparekmertables.push_back(&PARAM_SPACED_KMER_MODE);
        comparekmertables.push_back(&PARAM_SPACED_KMER_PATTERN);
        comparekmertables.push_back(&PARAM_MAX_SEQ_LEN);
        comparekmertables.push_back(&PARAM_REQ_KMER_MATCHES);
        comparekmertables.push_back(&PARAM_MAX_KMER_PER_POS);
        comparekmertables.push_back(&PARAM_NO_COMP_BIAS_CORR);
        comparekmertables.push_back(&PARAM_MASK_RESIDUES);
        comparekmertables.push_back(&PARAM_MASK_PROBABILTY);
        comparekmertables.push_back(&PARAM_COMPRESSED);
        comparekmertables.push_back(&PARAM_THREADS);
        comparekmertables.push_back(&PARAM_V);

        blockalign.push_back(&PARAM_E);
        blockalign.push_back(&PARAM_K);
        blockalign.push_back(&PARAM_SUB_MAT);
        blockalign.push_back(&PARAM_SPACED_KMER_MODE);
        blockalign.push_back(&PARAM_SPACED_KMER_PATTERN);
        blockalign.push_back(&PARAM_NO_COMP_BIAS_CORR);
        blockalign.push_back(&PARAM_GAP_OPEN);
        blockalign.push_back(&PARAM_GAP_EXTEND);
        blockalign.push_back(&PARAM_MAX_SEQ_LEN);
        blockalign.push_back(&PARAM_RANGE_MIN);
        blockalign.push_back(&PARAM_RANGE_MAX);
        blockalign.push_back(&PARAM_X_DROP);
        blockalign.push_back(&PARAM_COMPRESSED);
        blockalign.push_back(&PARAM_THREADS);
        blockalign.push_back(&PARAM_V);

        convertsraalignments.push_back(&PARAM_FORMAT_MODE);
        convertsraalignments.push_back(&PARAM_FORMAT_OUTPUT);
        convertsraalignments.push_back(&PARAM_THREADS);
        convertsraalignments = combineList(convertsraalignments, Parameters::convertalignments);

        petasearchworkflow = combineList(createkmertable, comparekmertables);
        petasearchworkflow = combineList(petasearchworkflow, blockalign);
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

        maxKmerPerPos = 20;

        rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    }

    LocalParameters(LocalParameters const &);

    ~LocalParameters() {};

    void operator=(LocalParameters const &);
};

#endif

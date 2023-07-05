// Created by Match on 8/13/2021.
#ifndef SRASEARCH_BLOCKALIGNER_H
#define SRASEARCH_BLOCKALIGNER_H

#include "Matcher.h"
#include "block_aligner.h"

class BlockAligner {
public:
    BlockAligner(
        size_t maxSequenceLength,
        uintptr_t min, uintptr_t max,
        int8_t gapOpen, int8_t gapExtend
    );

    ~BlockAligner();

    void initTarget(Sequence &target);

    Matcher::result_t align(
        Sequence &query,
        DistanceCalculator::LocalAlignment alignment,
        EvalueComputation *evaluer,
        int xdrop,
        BaseMatrix *subMat = nullptr,
        bool useProfile = false
    );

private:
    PaddedBytes* a;
    PaddedBytes* b;
    AAProfile* aProfile;
    AAProfile* bProfile;
    AAMatrix* matrix;
    BlockHandle blockTrace;
    BlockHandle blockNoTrace;
    Cigar* cigar;

    char* querySeqRev;

    size_t targetLength;
    const char* targetSeq;
    char* targetSeqRev;

    SizeRange range;
    Gaps gaps;
    const char PSSMAlphabet[20] = {'A', 'C', 'D', 'E', 'F',
                                   'G', 'H', 'I', 'K', 'L',
                                   'M', 'N', 'P', 'Q', 'R',
                                   'S', 'T', 'V', 'W', 'Y'};

    void initializeProfile(
        const int8_t *rawProfileMatrix,
        size_t seqStart,
        size_t seqEnd,
        size_t seqLen,
        AAProfile *result,
        bool reverse = false
    );
};

#endif

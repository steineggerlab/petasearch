//
// Created by Match on 8/13/2021.
//

#include <Matcher.h>
#include "block_aligner.h"

#ifndef SRASEARCH_BLOCKALIGNER_H
#define SRASEARCH_BLOCKALIGNER_H


class BlockAligner {
public:
    BlockAligner(size_t maxSequenceLength,
                 uintptr_t min, uintptr_t max,
                 int8_t gapOpen, int8_t gapExtend);

    ~BlockAligner();

    void initQuery(Sequence *query);

    Matcher::result_t align(Sequence *targetSeqObj,
                            DistanceCalculator::LocalAlignment alignment,
                            EvalueComputation *evaluer,
                            int xdrop,
                            BaseMatrix *subMat = nullptr,
                            bool useProfile = false);

private:
    char *targetSeqRev;
    char *querySeq;
    char *targetSeq;
    size_t querySeqLen;
    char *querySeqRev;
    SizeRange range;
    Gaps gaps{};
    const char PSSMAlphabet[20] = {'A', 'C', 'D', 'E', 'F',
                                   'G', 'H', 'I', 'K', 'L',
                                   'M', 'N', 'P', 'Q', 'R',
                                   'S', 'T', 'V', 'W', 'Y'};

    void initializeProfile(const int8_t *rawProfileMatrix,
                           size_t seqStart,
                           size_t seqEnd,
                           size_t seqLen,
                           AAProfile *result,
                           bool reverse = false);
};


#endif //SRASEARCH_BLOCKALIGNER_H

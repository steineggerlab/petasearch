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
                            int xdrop);

private:
    char *targetSeqRev;
    char *querySeq;
    char *targetSeq;
    size_t querySeqLen;
    char *querySeqRev;
    SizeRange range;
    Gaps gaps{};
};


#endif //SRASEARCH_BLOCKALIGNER_H

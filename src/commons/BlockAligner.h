//
// Created by Match on 8/13/2021.
//

#include <Matcher.h>
#include "block_aligner.h"

#ifndef SRASEARCH_BLOCKALIGNER_H
#define SRASEARCH_BLOCKALIGNER_H


class BlockAligner {
public:
    BlockAligner(BaseMatrix *m, int maxSeqLen, int gapOpen, int gapExtend);

    ~BlockAligner();

    void initQuery(Sequence *query);

    Matcher::result_t align(Sequence *targetSeqObj,
                            int diagonal,
                            EvalueComputation *evaluer);

private:
    SubstitutionMatrix::FastMatrix fastMatrix;
    uint8_t *targetSeqRev;
    int targetSeqRevDataLen;
    uint8_t *querySeq;
    uint8_t *querySeqRev;
    int querySeqRevDataLen;
    Sequence *querySeqObj;
    int8_t *mat;
    AAMatrix *subMat;
    SizeRange range;
    Gaps gaps;
};


#endif //SRASEARCH_BLOCKALIGNER_H

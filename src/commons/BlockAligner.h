//
// Created by Match on 8/13/2021.
//

#include <Matcher.h>
#include "block_aligner.h"

#ifndef SRASEARCH_BLOCKALIGNER_H
#define SRASEARCH_BLOCKALIGNER_H


class BlockAligner {
public:
    BlockAligner(BaseMatrix *m, int gapOpen, int gapExtend);

    ~BlockAligner();

    void initQuery(Sequence *query);

    Matcher::result_t align(Sequence *targetSeqObj,
                            int diagonal,
                            EvalueComputation *evaluer);

private:
    SubstitutionMatrix::FastMatrix fastMatrix;
    const char *queryID;
    char *targetSeqRev;
    const char *querySeq;
    int querySeqLen;
    char *querySeqRev;
    int8_t *mat;
    AAMatrix *subMat;
    SizeRange range;
    Gaps gaps;
};


#endif //SRASEARCH_BLOCKALIGNER_H

// Created by Match on 8/13/2021.
#ifndef SRASEARCH_BLOCKALIGNER_H
#define SRASEARCH_BLOCKALIGNER_H

#include "Matcher.h"
#include "DistanceCalculator.h"
#include "block_aligner.h"

class BlockAligner {
public:
    BlockAligner(
        size_t maxSequenceLength,
        uintptr_t min, uintptr_t max,
        int8_t gapOpen, int8_t gapExtend,
        int dbtype = Parameters::DBTYPE_AMINO_ACIDS
    );

    ~BlockAligner();

    void initTarget(Sequence &target);

    Matcher::result_t align(
        Sequence &query,
        DistanceCalculator::LocalAlignment alignment,
        EvalueComputation *evaluer,
        int xdrop,
        BaseMatrix *subMat = NULL
    );

private:
    PaddedBytes* a;
    PaddedBytes* b;
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
    int dbtype;
};

#endif

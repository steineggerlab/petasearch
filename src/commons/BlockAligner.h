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

    Matcher::result_t align(
        const char* targetSeq,
        unsigned int targetLength,
        const char* querySeq,
        unsigned int queryLength,
        DistanceCalculator::LocalAlignment alignment,
        EvalueComputation *evaluer,
        int xdrop,
        BaseMatrix& subMat
    );

private:
    PaddedBytes* a;
    PaddedBytes* b;
    AAProfile* bProfile;
    AAMatrix* matrix;
    BlockHandle blockTrace;
    BlockHandle blockNoTrace;
    Cigar* cigar;

    SizeRange range;
    Gaps gaps;
    int dbtype;
};

#endif

//
// Created by Match on 8/13/2021.
//

#include <DistanceCalculator.h>
#include <cstring>
#include "BlockAligner.h"
#include "Sequence.h"
#include "block_aligner.h"
#include "Util.h"
#include "Parameters.h"
#include "SRAUtil.h"

inline void swap(int &a, int &b) {
    int temp = b;
    b = a;
    a = temp;
}

BlockAligner::BlockAligner(size_t maxSequenceLength,
                           uintptr_t min = 32, uintptr_t max = 32,
                           int8_t gapOpen = -11, int8_t gapExtend = -1) :
        range({min, max}),
        gaps({gapOpen, gapExtend}) {
    targetSeq = static_cast<char *>(calloc(maxSequenceLength + 1, sizeof(char)));
    querySeq = static_cast<char *>(calloc(maxSequenceLength + 1, sizeof(char)));
    targetSeqRev = static_cast<char *>(calloc(maxSequenceLength + 1, sizeof(char)));
    querySeqRev = static_cast<char *>(calloc(maxSequenceLength + 1, sizeof(char)));
}

BlockAligner::~BlockAligner() {
    free(querySeq);
    free(targetSeq);
    free(querySeqRev);
    free(targetSeqRev);
}

void BlockAligner::initQuery(Sequence *query) {
    SRAUtil::stripInvalidChars(query->getSeqData(), querySeq);
    querySeqLen = strlen(querySeq); // query->L;
    SRAUtil::strrev(querySeqRev, querySeq, querySeqLen);
}


Matcher::result_t
BlockAligner::align(Sequence *targetSeqObj,
                    DistanceCalculator::LocalAlignment alignment,
                    EvalueComputation *evaluer,
                    int xdrop) {
    int aaIds = 0;
    std::string backtrace;

    SRAUtil::stripInvalidChars(targetSeqObj->getSeqData(), targetSeq);
    SRAUtil::strrev(targetSeqRev, targetSeq, targetSeqObj->L);

    unsigned int qUngappedEndPos, dbUngappedEndPos;

    unsigned int distanceToDiagonal = alignment.distToDiagonal;

    if (alignment.diagonal < 0) {
        qUngappedEndPos = alignment.endPos + distanceToDiagonal;
        dbUngappedEndPos = alignment.endPos;
    } else {
        qUngappedEndPos = alignment.endPos;
        dbUngappedEndPos = alignment.endPos + distanceToDiagonal;
    }

    // get middle position of ungapped alignment
    long tmp = ((long)querySeqLen - (long)qUngappedEndPos) - 1;
    unsigned int qStartRev = tmp < 0 ? 0 : tmp ; // - 1
    unsigned int qEndRev = querySeqLen;
    char *querySeqRevAlign = SRAUtil::substr(querySeqRev, qStartRev, qEndRev);
    size_t len_querySeqRevAlign = std::strlen(querySeqRevAlign);

    unsigned int tStartRev = (targetSeqObj->L - dbUngappedEndPos) - 1;
    unsigned int tEndRev = targetSeqObj->L;
    char *targetSeqRevAlign = SRAUtil::substr(targetSeqRev, tStartRev, tEndRev);
    size_t len_targetSeqRevAlign = std::strlen(targetSeqRevAlign);

    PaddedBytes *queryRevPadded = block_new_padded_aa(len_querySeqRevAlign, range.max);
    PaddedBytes *targetRevPadded = block_new_padded_aa(len_targetSeqRevAlign, range.max);

    block_set_bytes_padded_aa(queryRevPadded, (const uint8_t *)querySeqRevAlign, len_querySeqRevAlign, range.max);
    block_set_bytes_padded_aa(targetRevPadded, (const uint8_t *)targetSeqRevAlign, len_targetSeqRevAlign, range.max);

    BlockHandle blockRev = block_new_aa_trace_xdrop(len_querySeqRevAlign, len_targetSeqRevAlign, range.max);
    block_align_aa_trace_xdrop(blockRev, queryRevPadded, targetRevPadded, &BLOSUM62, gaps, range, xdrop);

    AlignResult resRev = block_res_aa_trace_xdrop(blockRev);

    unsigned int qStartPos = querySeqLen - (qStartRev + resRev.query_idx);
    unsigned int qEndPosAlign = querySeqLen;
    char *querySeqAlign = SRAUtil::substr(querySeq, qStartPos, qEndPosAlign);
    size_t len_querySeqAlign = std::strlen(querySeqAlign);

    unsigned int tStartPos = targetSeqObj->L - (tStartRev + resRev.reference_idx);
    unsigned int tEndPosAlign = targetSeqObj->L;
    char *targetSeqAlign = SRAUtil::substr(targetSeq, tStartPos, tEndPosAlign);
    size_t len_targetSeqAlign = std::strlen(targetSeqAlign);

//    PaddedBytes *queryPadded = block_make_padded_aa(querySeqAlign, range.max);
//    PaddedBytes *targetPadded = block_make_padded_aa(targetSeqAlign, range.max);
    PaddedBytes *queryPadded = block_new_padded_aa(len_querySeqAlign, range.max);
    PaddedBytes *targetPadded = block_new_padded_aa(len_targetSeqAlign, range.max);

    block_set_bytes_padded_aa(queryPadded, (const uint8_t *)querySeqAlign, len_querySeqAlign, range.max);
    block_set_bytes_padded_aa(targetPadded, (const uint8_t *)targetSeqAlign, len_targetSeqAlign, range.max);

    BlockHandle block = block_new_aa_trace_xdrop(len_querySeqAlign, len_targetSeqAlign, range.max);
    block_align_aa_trace_xdrop(block, queryPadded, targetPadded, &BLOSUM62, gaps, range, xdrop);

    AlignResult res = block_res_aa_trace_xdrop(block);

    bool reverseCigar = false;

    if (resRev.query_idx > res.query_idx && resRev.reference_idx > res.reference_idx) {
        res = block_res_aa_trace_xdrop(blockRev);
        reverseCigar = true;
    }

    Cigar *cigar = block_new_cigar(res.query_idx, res.reference_idx);
    if (reverseCigar) {
         block_cigar_aa_trace_xdrop(blockRev, res.query_idx, res.reference_idx, cigar);
    } else {
        block_cigar_aa_trace_xdrop(block, res.query_idx, res.reference_idx, cigar);
    }

    size_t cigarLen = block_len_cigar(cigar);

    Matcher::result_t realResult;
    //    result.cigar = retCigar;
    int bitScore = static_cast<int>(evaluer->computeBitScore(res.score) + 0.5);
    int qEndPos = qStartPos + res.query_idx - 1;
    int dbEndPos = tStartPos + res.reference_idx - 1;
    float qcov = SmithWaterman::computeCov(qStartPos, qEndPos, querySeqLen);
    float dbcov = SmithWaterman::computeCov(tStartPos, dbEndPos, targetSeqObj->L);
    double evalue = evaluer->computeEvalue(res.score, querySeqLen);

    char ops_char[] = {' ', 'M', 'I', 'D'};

    int alnLength = Matcher::computeAlnLength(qStartPos, qEndPos, tStartPos, dbEndPos);
    if (cigarLen > 0) {
        int32_t targetPos = 0, queryPos = 0;
        for (unsigned long c = 0; c < cigarLen; ++c) {
            OpLen cigar_op = block_get_cigar(cigar, c);
            char letter = ops_char[cigar_op.op];

            uint32_t length = cigar_op.len;

            backtrace.reserve(length);

            for (uint32_t i = 0; i < length; ++i) {
                if (letter == 'M') {
                    if (reverseCigar) {
                        if (targetSeqRevAlign[targetPos] == 'X' or querySeqRevAlign[queryPos] == 'X') {
                            aaIds++;
                        } else if (targetSeqRevAlign[targetPos] == querySeqRevAlign[queryPos]) {
                            aaIds++;
                        }
                    } else {
                        if (targetSeqAlign[targetPos] == 'X' or querySeqAlign[queryPos] == 'X') {
                            aaIds++;
                        } else if (targetSeqAlign[targetPos] == querySeqAlign[queryPos]) {
                            aaIds++;
                        }
                    }
                    ++queryPos;
                    ++targetPos;
                    backtrace.append("M");
                } else {
                    if (letter == 'I') {
                        ++queryPos;
                        backtrace.append("I");
                    } else {
                        ++targetPos;
                        backtrace.append("D");
                    }
                }
            }
        }
        alnLength = backtrace.length();
    }

    int seqIdMode = Parameters::SEQ_ID_ALN_LEN;

    float seqId = Util::computeSeqId(seqIdMode, aaIds, querySeqLen, targetSeqObj->L, alnLength);
    if (reverseCigar) {
        std::reverse(backtrace.begin(), backtrace.end());
    }

    realResult = Matcher::result_t(targetSeqObj->getDbKey(), bitScore, qcov, dbcov, seqId, evalue, alnLength,
                                   qStartPos, qEndPos, querySeqLen, tStartPos, dbEndPos, targetSeqObj->L,
                                   backtrace);

    block_free_padded_aa(queryRevPadded);
    block_free_padded_aa(targetRevPadded);
    block_free_padded_aa(queryPadded);
    block_free_padded_aa(targetPadded);
    block_free_cigar(cigar);
    block_free_aa_trace_xdrop(block);
    block_free_aa_trace_xdrop(blockRev);

    free(querySeqRevAlign);
    free(targetSeqRevAlign);
    free(querySeqAlign);
    free(targetSeqAlign);
    return realResult;
}

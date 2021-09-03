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

inline void swap(int &a, int &b) {
    int temp = b;
    b = a;
    a = temp;
}

void strrev(char *strRev, const char *str, int len) {
    int start = 0;
    while (LIKELY(start <= len)) {
        strRev[start] = str[len];
        strRev[len] = str[start];
        ++start;
        --len;
    }
}

/**
 * @brief Make a slice origStr[start:end], start inclusive, end exclusive
 * */
char *substr(char *origStr, int start, int end) {
    char *subStr = static_cast<char *>(calloc(end - start + 1, sizeof(char)));
    strncpy(subStr, origStr + start, end - start);
    return subStr;
}

BlockAligner::BlockAligner(size_t maxSequenceLength, int8_t gapOpen, int8_t gapExtend) : range({32, 32}) {
    targetSeqRev = static_cast<char *>(calloc(maxSequenceLength + 1, sizeof(char)));
    querySeqRev = static_cast<char *>(calloc(maxSequenceLength + 1, sizeof(char)));
    range.min = 32;
    range.max = 4096;
    gaps = {gapOpen, gapExtend};
}

BlockAligner::~BlockAligner() {
    free(querySeqRev);
    free(targetSeqRev);
}

void BlockAligner::initQuery(Sequence *query) {
    querySeq = query->getSeqData();
    querySeqLen = query->L;
//    querySeqRev = static_cast<char *>(calloc(query->L + 1, sizeof(char)));
    strrev(querySeqRev, querySeq, querySeqLen);
}


Matcher::result_t
BlockAligner::align(Sequence *targetSeqObj, DistanceCalculator::LocalAlignment alignment, EvalueComputation *evaluer,
                    int
                    xdrop) {
    int aaIds = 0;

    // TODO: make ungapped alignment result pass in as parameter
//    int xdrop = 50;

    std::string backtrace;

    const char *targetSeq = targetSeqObj->getSeqData();
//    targetSeqRev = static_cast<char *>(calloc(targetSeqObj->L + 1, sizeof(char)));
    strrev(targetSeqRev, targetSeq, targetSeqObj->L - 1);

    int qUngappedStartPos, qUngappedEndPos, dbUngappedStartPos, dbUngappedEndPos;

//    DistanceCalculator::LocalAlignment alignment;
//    alignment = DistanceCalculator::computeUngappedAlignment(
//            querySeq, querySeqLen, targetSeqObj->getSeqData(), targetSeqObj->L,
//            diagonal, fastMatrix.matrix, Parameters::RESCORE_MODE_ALIGNMENT);


    unsigned int distanceToDiagonal = alignment.distToDiagonal;
//    diagonal = alignment.diagonal;

    if (alignment.diagonal < 0) {
        qUngappedStartPos = alignment.startPos + distanceToDiagonal;
        qUngappedEndPos = alignment.endPos + distanceToDiagonal;
        dbUngappedStartPos = alignment.startPos;
        dbUngappedEndPos = alignment.endPos;
    } else {
        qUngappedStartPos = alignment.startPos;
        qUngappedEndPos = alignment.endPos;
        dbUngappedStartPos = alignment.startPos + distanceToDiagonal;
        dbUngappedEndPos = alignment.endPos + distanceToDiagonal;
    }

    Debug(Debug::INFO) << "qUngappedStartPos = " << qUngappedStartPos << "\n"
                       << "qUngappedEndPos = " << qUngappedEndPos << "\n"
                       << "dbUngappedStartPos = " << dbUngappedStartPos << "\n"
                       << "dbUngappedEndPos = " << dbUngappedEndPos << "\n\n";

    // get middle position of ungapped alignment
    if (qUngappedStartPos > qUngappedEndPos) {
        swap(qUngappedStartPos, qUngappedEndPos);
    }
    int qStartRev = (querySeqLen - qUngappedEndPos) - 1;
    int qEndRev = querySeqLen;
//    if (qEndRev < qStartRev) {
//        swap(qStartRev, qEndRev);
//    }
    Debug(Debug::INFO) << "qStartRev: " << qStartRev << "\n"
                       << "qEndRev：" << qEndRev << "\n";
    // TODO: add comments for substr
    char *querySeqRevAlign = substr(querySeqRev, qStartRev, qEndRev);
//    Debug(Debug::INFO) << "querySeqRevAlign: " << querySeqRevAlign << "\n\n";

    if (dbUngappedStartPos > dbUngappedEndPos) {
        swap(dbUngappedStartPos, dbUngappedEndPos);
    }
    int tStartRev = (targetSeqObj->L - dbUngappedEndPos) - 1;
    int tEndRev = targetSeqObj->L;
    if (tEndRev < tStartRev) {
        swap(tStartRev, tEndRev);
    }
    Debug(Debug::INFO) << "tStartRev: " << tStartRev << "\n"
                       << "tEndRev：" << tEndRev << "\n";
    char *targetSeqRevAlign = substr(targetSeqRev, tStartRev, tEndRev);
//    Debug(Debug::INFO) << "targetSeqRevAlign: " << targetSeqRevAlign << "\n\n";

    PaddedBytes *queryRevPadded = block_make_padded_aa(querySeqRevAlign, range.max);
    PaddedBytes *targetRevPadded = block_make_padded_aa(targetSeqRevAlign, range.max);
    BlockHandle blockRev = block_align_aa_trace_xdrop(queryRevPadded, targetRevPadded, &BLOSUM62, gaps, range, xdrop);
    AlignResult resRev = block_res_aa_trace_xdrop(blockRev);

    int qStartPos = (int) querySeqLen - (qStartRev + resRev.query_idx);
    int qEndPosAlign = querySeqLen;
//    if (qEndPosAlign < qStartPos) {
//        swap(qStartPos, qEndPosAlign);
//    }
    Debug(Debug::INFO) << "qStart: " << qStartPos << "\n"
                       << "qEnd：" << qEndPosAlign << "\n";
    char *querySeqAlign = substr((char *) querySeq, qStartPos, qEndPosAlign);
//    Debug(Debug::INFO) << "querySeqAlign: " << querySeqAlign << "\n\n";


    int tStartPos = targetSeqObj->L - (tStartRev + resRev.reference_idx);
    int tEndPosAlign = targetSeqObj->L;
//    if (tEndPosAlign < tStartPos) {
//        swap(tStartPos, tEndPosAlign);
//    }
    Debug(Debug::INFO) << "tStart: " << tStartPos << "\n"
                       << "tEnd：" << tEndPosAlign << "\n";
    char *targetSeqAlign = substr((char *) targetSeq, tStartPos, tEndPosAlign);
//    Debug(Debug::INFO) << "targetSeqAlign: " << targetSeqAlign << "\n\n";

    PaddedBytes *queryPadded = block_make_padded_aa(querySeqAlign, range.max);
    PaddedBytes *targetPadded = block_make_padded_aa(targetSeqAlign, range.max);
    BlockHandle block = block_align_aa_trace_xdrop(queryPadded, targetPadded, &BLOSUM62, gaps, range, xdrop);
    AlignResult res = block_res_aa_trace_xdrop(block);

    bool reverseCigar = false;

    if (resRev.query_idx > res.query_idx && resRev.reference_idx > res.reference_idx) {
        res = block_res_aa_trace_xdrop(blockRev);
        reverseCigar = true;
    }

    CigarVec cigar;
    if (reverseCigar) {
        cigar = block_cigar_aa_trace_xdrop(blockRev);
    } else {
        cigar = block_cigar_aa_trace_xdrop(block);
    }

    Matcher::result_t realResult;
    //    result.cigar = retCigar;
    int32_t cigarLen = cigar.len;
    int bitScore = static_cast<int>(evaluer->computeBitScore(res.score) + 0.5);
    int qEndPos = qStartPos + res.query_idx - 1;
    int dbEndPos = tStartPos + res.reference_idx - 1;
    float qcov = SmithWaterman::computeCov(qStartPos, qEndPos - 1, querySeqLen);
    float dbcov = SmithWaterman::computeCov(tStartPos, dbEndPos - 1, targetSeqObj->L);
    double evalue = evaluer->computeEvalue(res.score, querySeqLen);

    char ops_char[] = {' ', 'M', 'I', 'D'};

    unsigned int alnLength = Matcher::computeAlnLength(qStartPos, qEndPos, tStartPos, dbEndPos);
    if (cigar.len > 0) {
        int32_t targetPos = 0, queryPos = 0;
        for (int32_t c = 0; c < cigarLen; ++c) {
            char letter = ops_char[cigar.ptr[c].op];
            uint32_t length = cigar.ptr[c].len;
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

    unsigned int seqIdMode = Parameters::SEQ_ID_ALN_LEN;

    float seqId = Util::computeSeqId(seqIdMode, aaIds, querySeqLen, targetSeqObj->L, alnLength);
    Debug(Debug::INFO) << backtrace << "\n";
    if (reverseCigar) {
        std::reverse(backtrace.begin(), backtrace.end());
        Debug(Debug::INFO) << backtrace << "\n";
    }

    Debug(Debug::INFO) << "aaIds: " << aaIds << "\n";

    realResult = Matcher::result_t(targetSeqObj->getDbKey(), bitScore, qcov, dbcov, seqId, evalue, alnLength,
                                   qStartPos, qEndPos, querySeqLen, tStartPos, dbEndPos, targetSeqObj->L, backtrace);

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
    //    free(ezAlign.cigar);
}

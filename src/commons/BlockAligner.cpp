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

/**
 * @brief Reverse a string and put the results in strRev
 * @param strRev the destination to put the reversed string
 * @param str the string to be revesed
 * @param len the length of the string
 */
void strrev(char *strRev, const char *str, int len) {
    int start = 0;
    int end = len - 1;
    while (LIKELY(start <= end)) {
        strRev[start] = str[end];
        strRev[end] = str[start];
        ++start;
        --end;
    }
    strRev[len] = '\0';
}

/**
 * @brief Make a slice origStr[start:end], start inclusive, end exclusive
 * */
char *substr(char *origStr, unsigned int start, unsigned int end) {
    char *subStr = static_cast<char *>(calloc(end - start + 1, sizeof(char)));
    strncpy(subStr, origStr + start, end - start);
    return subStr;
}

/**
 * @brief Strip invalid characters from a string (@, *, newline, tab, etc)
 */
void stripInvalidChars(const char *src, char *dest) {
    size_t j, n = strlen(src);
    for (size_t i = j = 0; i < n; i++) {
        // FIXME: this branch is likely not useful now
        if (src[i] == '\n' || src[i] == '@') {
            continue;
        } else if (src[i] == '*') {
            dest[j++] = 'X';
        } else {
            dest[j++] = src[i];
        }
    }
    dest[j] = '\0';
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
    stripInvalidChars(query->getSeqData(), querySeq);
    querySeqLen = strlen(querySeq); // query->L;
    strrev(querySeqRev, querySeq, querySeqLen);
}


Matcher::result_t
BlockAligner::align(Sequence *targetSeqObj,
                    DistanceCalculator::LocalAlignment alignment,
                    EvalueComputation *evaluer,
                    int xdrop) {
    int aaIds = 0;
    std::string backtrace;

    stripInvalidChars(targetSeqObj->getSeqData(), targetSeq);
    strrev(targetSeqRev, targetSeq, targetSeqObj->L);

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

    Debug(Debug::INFO) << "querySeqRev: " << querySeqRev << "\n";
    Debug(Debug::INFO) << "qUngappedEndPos: " << qUngappedEndPos << " tmp: " << tmp << "\n";
    unsigned int qStartRev = tmp < 0 ? 0 : tmp ; // - 1
    unsigned int qEndRev = querySeqLen;
    char *querySeqRevAlign = substr(querySeqRev, qStartRev, qEndRev);
    Debug(Debug::INFO) << "qStartRev: " << qStartRev << " qEndRev: " << qEndRev << "\n";
    Debug(Debug::INFO) << "querySeqRevAlign: " << querySeqRevAlign << "\n";

    unsigned int tStartRev = (targetSeqObj->L - dbUngappedEndPos) - 1;
    unsigned int tEndRev = targetSeqObj->L;
    char *targetSeqRevAlign = substr(targetSeqRev, tStartRev, tEndRev);

    PaddedBytes *queryRevPadded = block_make_padded_aa(querySeqRevAlign, range.max);
    PaddedBytes *targetRevPadded = block_make_padded_aa(targetSeqRevAlign, range.max);
    BlockHandle blockRev = block_align_aa_trace_xdrop(queryRevPadded, targetRevPadded, &BLOSUM62, gaps, range, xdrop);
    AlignResult resRev = block_res_aa_trace_xdrop(blockRev);

    unsigned int qStartPos = querySeqLen - (qStartRev + resRev.query_idx);
    unsigned int qEndPosAlign = querySeqLen;
    char *querySeqAlign = substr(querySeq, qStartPos, qEndPosAlign);

    unsigned int tStartPos = targetSeqObj->L - (tStartRev + resRev.reference_idx);
    unsigned int tEndPosAlign = targetSeqObj->L;
    char *targetSeqAlign = substr(targetSeq, tStartPos, tEndPosAlign);

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
    unsigned long cigarLen = cigar.len;
    int bitScore = static_cast<int>(evaluer->computeBitScore(res.score) + 0.5);
    int qEndPos = qStartPos + res.query_idx - 1;
    int dbEndPos = tStartPos + res.reference_idx - 1;
    float qcov = SmithWaterman::computeCov(qStartPos, qEndPos, querySeqLen);
    float dbcov = SmithWaterman::computeCov(tStartPos, dbEndPos, targetSeqObj->L);
    double evalue = evaluer->computeEvalue(res.score, querySeqLen);

    char ops_char[] = {' ', 'M', 'I', 'D'};

    int alnLength = Matcher::computeAlnLength(qStartPos, qEndPos, tStartPos, dbEndPos);
    if (cigar.len > 0) {
        int32_t targetPos = 0, queryPos = 0;
        for (unsigned long c = 0; c < cigarLen; ++c) {
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
    //    free(ezAlign.cigar);
}

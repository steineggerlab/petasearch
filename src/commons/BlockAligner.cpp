//
// Created by Match on 8/13/2021.
//

#include <DistanceCalculator.h>
#include <stdio.h>
#include <string.h>
#include "BlockAligner.h"
#include "Sequence.h"
#include "block_aligner.h"
#include "Util.h"

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

char *substr(char *origStr, int start, int end) {
    char *subStr = static_cast<char *>(calloc(end - start + 1, sizeof(char)));
    strncpy(subStr, origStr, end - start);
    return subStr;
}

BlockAligner::BlockAligner(BaseMatrix *subMat, int gapOpen, int gapExtend) :
        fastMatrix(SubstitutionMatrix::createAsciiSubMat(*subMat)) {
    mat = new int8_t[subMat->alphabetSize * subMat->alphabetSize];
    this->subMat = (AAMatrix *) subMat;
    for (int i = 0; i < subMat->alphabetSize; i++) {
        for (int j = 0; j < subMat->alphabetSize; j++) {
            mat[i * subMat->alphabetSize + j] = subMat->subMatrix[i][j];
        }
    }

    range.min = 32;
    range.max = 4096;

    gaps.extend = -gapExtend;
    gaps.open = -gapOpen;
}

BlockAligner::~BlockAligner() {
    delete[] fastMatrix.matrixData;
    delete[] fastMatrix.matrix;
    delete[] mat;
}

void BlockAligner::initQuery(Sequence *query) {
    querySeq = query->getSeqData();
    querySeqLen = query->L;

    querySeqRev = static_cast<char *>(calloc(query->L + 1, sizeof(char)));

    strrev(querySeqRev, querySeq, querySeqLen - 1);
}


Matcher::result_t BlockAligner::align(Sequence *targetSeqObj,
                                      int diagonal,
                                      EvalueComputation *evaluer) {
    int aaIds = 0;

    // TODO: make this pass in as a value
    // TODO: make ungapped alignment result pass in as parameter
    int xdrop = 2048;

    std::string backtrace;

    const char *targetSeq = targetSeqObj->getSeqData();
    targetSeqRev = static_cast<char *>(calloc(targetSeqObj->L + 1, sizeof(char)));
    strrev(targetSeqRev, targetSeq, targetSeqObj->L - 1);

    int qUngappedStartPos, qUngappedEndPos, dbUngappedStartPos, dbUngappedEndPos;

    DistanceCalculator::LocalAlignment alignment;
    int queryLen = querySeqLen;
    int origQueryLen = queryLen;
    alignment = DistanceCalculator::computeUngappedAlignment(
            querySeq, querySeqLen, targetSeqObj->getSeqData(), targetSeqObj->L,
            diagonal, fastMatrix.matrix, Parameters::RESCORE_MODE_ALIGNMENT);


    unsigned int distanceToDiagonal = alignment.distToDiagonal;
    diagonal = alignment.diagonal;

    if (diagonal >= 0) {
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

    int qStartPos = querySeqLen - (qStartRev + resRev.query_idx);
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
    int qEndPos = qStartPos + res.query_idx;
    int dbEndPos = tStartPos + res.reference_idx;
    int qcov = SmithWaterman::computeCov(qStartPos, qEndPos, querySeqLen);
    int dbcov = SmithWaterman::computeCov(tStartPos, dbEndPos, targetSeqObj->L);
    double evalue = evaluer->computeEvalue(res.score, origQueryLen);

    char ops_char[] = {' ', 'M', 'I', 'D'};

    unsigned int alnLength = Matcher::computeAlnLength(qStartPos, qEndPos, tStartPos, dbEndPos);
    if (cigar.len > 0) {
        int32_t targetPos = tStartPos, queryPos = qStartPos;
        for (int32_t c = 0; c < cigarLen; ++c) {
            char letter = ops_char[cigar.ptr[c].op];
            uint32_t length = cigar.ptr[c].len;
            backtrace.reserve(length);

            for (uint32_t i = 0; i < length; ++i) {
                if (letter == 'M') {
                    if (targetSeq[targetPos] == querySeq[queryPos]) {
                        aaIds++;
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

    unsigned int seqIdMode = Matcher::SCORE_COV_SEQID;

    float seqId = Util::computeSeqId(seqIdMode, aaIds, querySeqLen, targetSeqObj->L, alnLength);

    if (reverseCigar == false) {
        std::reverse(backtrace.begin(), backtrace.end());
    }

    realResult = Matcher::result_t(targetSeqObj->getDbKey(), bitScore, qcov, dbcov, seqId, evalue, alnLength,
                                   qStartPos, qEndPos, origQueryLen, tStartPos, dbEndPos, targetSeqObj->L, backtrace);

    block_free_aa_trace_xdrop(block);
    block_free_aa_trace_xdrop(blockRev);
    return realResult;
    //    free(ezAlign.cigar);
}

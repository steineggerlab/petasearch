//
// Created by Match on 8/13/2021.
//

#include <DistanceCalculator.h>
#include "BlockAligner.h"
#include "Sequence.h"

BlockAligner::BlockAligner(BaseMatrix *subMat, int maxSeqLen, int gapOpen, int gapExtend) :
        fastMatrix(SubstitutionMatrix::createAsciiSubMat(*subMat)) {
    targetSeqRevDataLen = maxSeqLen;
    targetSeqRev = static_cast<uint8_t *>(malloc(targetSeqRevDataLen + 1));

    querySeqRevDataLen = maxSeqLen;
    querySeqRev = static_cast<uint8_t *>(malloc(querySeqRevDataLen + 1));

    mat = new int8_t[subMat->alphabetSize * subMat->alphabetSize];
    this->subMat = (AAMatrix *) subMat;
    for (int i = 0; i < subMat->alphabetSize; i++) {
        for (int j = 0; j < subMat->alphabetSize; j++) {
            mat[i * subMat->alphabetSize + j] = subMat->subMatrix[i][j];
        }
    }

    gaps.extend = gapExtend;
    gaps.open = gapOpen;
}

BlockAligner::~BlockAligner() {
    free(targetSeqRev);
    free(querySeqRev);
    delete[] fastMatrix.matrixData;
    delete[] fastMatrix.matrix;
    delete[] mat;
}

void BlockAligner::initQuery(Sequence *query) {
    querySeqObj = query;
    querySeq = query->numSequence;
    if (query->L >= querySeqRevDataLen) {
        querySeqRev = static_cast<uint8_t *>(realloc(querySeqRev, query->L + 1));
        querySeqRevDataLen = query->L;
    }
    SmithWaterman::seq_reverse((int8_t *) querySeqRev, (int8_t *) querySeq, query->L);
}


Matcher::result_t BlockAligner::align(Sequence *targetSeqObj,
                                      int diagonal,
                                      EvalueComputation *evaluer) {
    char *queryCharSeqAlign = (char *) querySeqObj->getSeqData();
    uint8_t *querySeqRevAlign = querySeqRev;
    uint8_t *querySeqAlign = querySeq;

    int aaIds = 0;

    std::string backtrace;

    const unsigned char *targetSeq = targetSeqObj->numSequence;
    if (targetSeqObj->L >= targetSeqRevDataLen) {
        targetSeqRev = static_cast<uint8_t *>(realloc(targetSeqRev, targetSeqObj->L + 1));
        targetSeqRevDataLen = targetSeqObj->L;
    }
    SmithWaterman::seq_reverse((int8_t *) targetSeqRev, (int8_t *) targetSeq, targetSeqObj->L);

    int qUngappedStartPos, qUngappedEndPos, dbUngappedStartPos, dbUngappedEndPos;

    DistanceCalculator::LocalAlignment alignment;
    int queryLen = querySeqObj->L;
    int origQueryLen = queryLen;
    alignment = DistanceCalculator::computeUngappedAlignment(
            queryCharSeqAlign, querySeqObj->L, targetSeqObj->getSeqData(), targetSeqObj->L,
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

//        int32_t cigarLen;
//        int bitScore, qEndPos, dbEndPos, qcov, dbcov;
//        double evalue;
//        unsigned int alnLenth;
//
//    // identical
//    if (qUngappedEndPos - qUngappedStartPos == origQueryLen - 1
//        && dbUngappedStartPos == 0 && dbUngappedEndPos == targetSeqObj->L - 1) {
//        s_align result;
//        uint32_t *retCigar = new uint32_t[1];
//        retCigar[0] = 0;
//        retCigar[0] = origQueryLen << 4;
//        result.cigar = retCigar;
//        result.cigarLen = 1;
//        result.score1 = alignment.score;
//        result.qStartPos1 = qUngappedStartPos;
//        result.qEndPos1 = qUngappedEndPos;
//        result.dbEndPos1 = dbUngappedEndPos;
//        result.dbStartPos1 = dbUngappedStartPos;
//        result.qCov = SmithWaterman::computeCov(result.qStartPos1, result.qEndPos1, querySeqObj->L);
//        result.tCov = SmithWaterman::computeCov(result.dbStartPos1, result.dbEndPos1, targetSeqObj->L);
//        result.evalue = evaluer->computeEvalue(result.score1, origQueryLen);
//        for (int i = qUngappedStartPos; i <= qUngappedEndPos; i++) {
//            aaIds += (querySeqAlign[i] == targetSeq[dbUngappedStartPos + (i - qUngappedStartPos)]) ? 1 : 0;
//        }
//        for (int pos = 0; pos < origQueryLen; pos++) {
//            backtrace.append("M");
//        }
//        return result;
//    }

    // get middle position of ungapped alignment
    int qStartRev = (querySeqObj->L - qUngappedEndPos) - 1;
    int tStartRev = (targetSeqObj->L - dbUngappedEndPos) - 1;

//    ksw_extz_t ez;
//    int flag = 0;
//    flag |= KSW_EZ_SCORE_ONLY;
//    flag |= KSW_EZ_EXTZ_ONLY;

//    int queryRevLenToAlign = querySeqObj->L - qStartRev;

    // void ksw_extz2_sse(
    // void *km, int qlen, const uint8_t *query, int tlen,
    // const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int flag,
    // ksw_extz_t *ez)
//    ksw_extz2_sse(0, queryRevLenToAlign, querySeqRevAlign + qStartRev, targetSeqObj->L - tStartRev,
//                  targetSeqRev + tStartRev, 5, mat, gapo, gape, 64, zdrop, flag, &ez);


    // TODO: fix this reinterpret_cast issue
    PaddedBytes *queryRevPadded = block_make_padded_aa(reinterpret_cast<const char *>(querySeqRevAlign + qStartRev),
                                                       range.max);
    PaddedBytes *targetRevPadded = block_make_padded_aa(reinterpret_cast<const char *>(targetSeqRev + tStartRev),
                                                        range.max);
    BlockHandle blockRev = block_align_aa_trace(queryRevPadded, targetRevPadded, subMat, gaps, range);
    AlignResult resRev = block_res_aa_trace(blockRev);

    int qStartPos = querySeqObj->L - (qStartRev + resRev.query_idx) - 1;
    int tStartPos = targetSeqObj->L - (tStartRev + resRev.reference_idx) - 1;

//    int alignFlag = 0;
//    alignFlag |= KSW_EZ_EXTZ_ONLY;

//    ksw_extz_t Align;
//    memset(&ezAlign, 0, sizeof(ksw_extz_t));

//    int queryLenToAlign = querySeqObj->L - qStartPos;

//    ksw_extz2_sse(0, queryLenToAlign, querySeqAlign + qStartPos, targetSeqObj->L - tStartPos, targetSeq + tStartPos, 5,
//                  mat, gapo, gape, 64, zdrop, alignFlag, &ezAlign);

    PaddedBytes *queryPadded = block_make_padded_aa(reinterpret_cast<const char *>(querySeqAlign + qStartPos),
                                                    range.max);
    PaddedBytes *targetPadded = block_make_padded_aa(reinterpret_cast<const char *>(targetSeq + tStartPos),
                                                     range.max);
    BlockHandle block = block_align_aa_trace(queryPadded, targetPadded, subMat, gaps, range);
    AlignResult res = block_res_aa_trace(block);

    std::string letterCode = "MID";
    // uint32_t *retCigar;

//    if (resRev.query_idx > res.query_idx && resRev.reference_idx > res.reference_idx) {
//        Debug(Debug::ERROR) << "I don't understand what to do here\n";
//        EXIT(EXIT_FAILURE);
//
//        ksw_extz2_sse(0, queryRevLenToAlign, querySeqRevAlign + qStartRev, targetSeqObj->L - tStartRev,
//                      targetSeqRev + tStartRev, 5, mat, gapo, gape, 64, zdrop, alignFlag, &ezAlign);
//
//        retCigar = new uint32_t[ezAlign.n_cigar];
//        for (int i = 0; i < ezAlign.n_cigar; i++) {
//            retCigar[i] = ezAlign.cigar[ezAlign.n_cigar - 1 - i];
//        }
//    } else {
//        retCigar = new uint32_t[ezAlign.n_cigar];
//        for (int i = 0; i < ezAlign.n_cigar; i++) {
//            retCigar[i] = ezAlign.cigar[i];
//        }
//    }

    CigarVec cigar = block_cigar_aa_trace(block);

    Matcher::result_t realResult;
//    result.cigar = retCigar;
    int32_t cigarLen = cigar.len;
    int bitScore = static_cast<int>(evaluer->computeBitScore(res.score) + 0.5);
    int qEndPos = qStartPos + res.query_idx;
    int dbEndPos = tStartPos + res.reference_idx;
    int qcov = SmithWaterman::computeCov(qStartPos, qEndPos, querySeqObj->L);
    int dbcov = SmithWaterman::computeCov(tStartPos, dbEndPos, targetSeqObj->L);
    double evalue = evaluer->computeEvalue(res.score, origQueryLen);

    unsigned int alnLength = Matcher::computeAlnLength(qStartPos, qEndPos, tStartPos, dbEndPos);
    if (cigar.len > 0) {
        int32_t targetPos = tStartPos, queryPos = qStartPos;
        for (int32_t c = 0; c < cigarLen; ++c) {
            char letter = SmithWaterman::cigar_int_to_op(cigar.ptr[c].op);
            uint32_t length = cigar.ptr[c].len;
            backtrace.reserve(length);

            for (uint32_t i = 0; i < length; ++i) {
                if (letter == 'M') {
                    if (targetSeq[targetPos] == querySeqAlign[queryPos]) {
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

    float seqId = Util::computeSeqId(seqIdMode, aaIds, origQueryLen, targetSeqObj->L, alnLength);

    realResult = Matcher::result_t(targetSeqObj->getDbKey(), bitScore, qcov, dbcov, seqId, evalue, alnLength,
                                   qStartPos, qEndPos, origQueryLen, tStartPos, dbEndPos, targetSeqObj->L, backtrace);

    return realResult;
//    free(ezAlign.cigar);
}




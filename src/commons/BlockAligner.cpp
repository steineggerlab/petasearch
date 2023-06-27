// Created by Match on 8/13/2021.
#include "DistanceCalculator.h"
#include "BlockAligner.h"
#include "Sequence.h"
#include "block_aligner.h"
#include "Util.h"
#include "Parameters.h"
#include "SRAUtil.h"

BlockAligner::BlockAligner(size_t maxSequenceLength,
                           uintptr_t min = 32, uintptr_t max = 32,
                           int8_t gapOpen = -11, int8_t gapExtend = -1) :
        range({min, max}),
        gaps({gapOpen, gapExtend}) {
    targetSeq = static_cast<char *>(calloc(maxSequenceLength + 1, sizeof(char)));
    querySeq = static_cast<char *>(calloc(maxSequenceLength + 1, sizeof(char)));
    targetSeqRev = static_cast<char *>(calloc(maxSequenceLength + 1, sizeof(char)));
    querySeqRev = static_cast<char *>(calloc(maxSequenceLength + 1, sizeof(char)));
    block = block_new_aa_trace_xdrop(maxSequenceLength + 1, maxSequenceLength + 1, range.max);
    blockRev = block_new_aa_trace_xdrop(maxSequenceLength + 1, maxSequenceLength + 1, range.max);
}

BlockAligner::~BlockAligner() {
    free(querySeq);
    free(targetSeq);
    free(querySeqRev);
    free(targetSeqRev);
    block_free_aa_trace_xdrop(block);
    block_free_aa_trace_xdrop(blockRev);
}

void BlockAligner::initQuery(Sequence *query) {
    // SRAUtil::stripInvalidChars(query->getSeqData(), querySeq);
    memcpy(querySeq, query->getSeqData(), query->L);
    querySeqLen = query->L;
    SRAUtil::strrev(querySeqRev, querySeq, querySeqLen);
}


void BlockAligner::initializeProfile(const int8_t *rawProfileMatrix,
                                     size_t seqStart,
                                     size_t seqEnd,
                                     size_t seqLen,
                                     AAProfile *result,
                                     bool reverse) {
    if (reverse) {
        for (size_t i = seqStart; i < seqEnd; i++) {
            for (size_t j = 0; j < 20; j++) {
                block_set_aaprofile(result, i - seqStart, PSSMAlphabet[j], rawProfileMatrix[seqLen - 1 - i + j * seqLen]);
            }
        }

    } else {
        for (size_t i = seqStart; i < seqEnd; i++) { // iterate through position
            for (size_t j = 0; j < 20; j++) { // iterate through alphabet
                block_set_aaprofile(result, i - seqStart, PSSMAlphabet[j], rawProfileMatrix[i + j * seqLen]);
            }
        }
    }

}


Matcher::result_t
BlockAligner::align(Sequence *targetSeqObj,
                    DistanceCalculator::LocalAlignment alignment,
                    EvalueComputation *evaluer,
                    int xdrop,
                    BaseMatrix *subMat,
                    bool useProfile) {

    int aaIds = 0;
    std::string backtrace;
    const int8_t *rawProfileMatrix = nullptr;

    if (useProfile && subMat == nullptr) {
        Debug(Debug::ERROR) << "Must provide subMat to BlockAligner if useProfile is set to true!\n";
        EXIT(EXIT_FAILURE);
    }

    if (useProfile) {
        rawProfileMatrix = targetSeqObj->getAlignmentProfile();
        SRAUtil::stripInvalidChars(
                SRAUtil::extractProfileSequence(targetSeqObj->getSeqData(), targetSeqObj->L, subMat).c_str(),
                targetSeq
        );
    } else {
        memcpy(targetSeq, targetSeqObj->getSeqData(), targetSeqObj->L);
        // SRAUtil::stripInvalidChars(targetSeqObj->getSeqData(), targetSeq);
    }
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
    long tmp = ((long) querySeqLen - (long) qUngappedEndPos) - 1;
    unsigned int qStartRev = tmp < 0 ? 0 : tmp; // - 1
    unsigned int qEndRev = querySeqLen;
    char *querySeqRevAlign = querySeqRev + qStartRev;
    size_t len_querySeqRevAlign = qEndRev - qStartRev;
    PaddedBytes *queryRevPadded = block_new_padded_aa(len_querySeqRevAlign, range.max);
    block_set_bytes_padded_aa(queryRevPadded, (const uint8_t *) querySeqRevAlign, len_querySeqRevAlign, range.max);

    unsigned int tStartRev = (targetSeqObj->L - dbUngappedEndPos) - 1;
    unsigned int tEndRev = targetSeqObj->L;
    // This is for sequence alignment
    size_t len_targetSeqRevAlign = tEndRev - tStartRev;
    char *targetSeqRevAlign = targetSeqRev + tStartRev;

    // profile to PSSM with specific range
    AAProfile *targetRevProfile = nullptr;
    PaddedBytes *targetRevPadded = nullptr;
    if (useProfile) {
        targetRevProfile = block_new_aaprofile(len_targetSeqRevAlign, range.max, gaps.extend);
        initializeProfile(rawProfileMatrix, tStartRev, tEndRev, targetSeqObj->L, targetRevProfile, true);
        for (size_t i = 0; i < len_targetSeqRevAlign; i++) {
            block_set_gap_open_C_aaprofile(targetRevProfile, i, gaps.open);
            block_set_gap_close_C_aaprofile(targetRevProfile, i, 0);
            block_set_gap_open_R_aaprofile(targetRevProfile, i, gaps.open);
        }
    } else {
        targetRevPadded = block_new_padded_aa(len_targetSeqRevAlign, range.max);
        block_set_bytes_padded_aa(targetRevPadded, (const uint8_t *) targetSeqRevAlign, len_targetSeqRevAlign, range.max);
    }

    if (useProfile) {
        block_align_profile_aa_trace_xdrop(blockRev, queryRevPadded, targetRevProfile, range, xdrop);
    } else {
        block_align_aa_trace_xdrop(blockRev, queryRevPadded, targetRevPadded, &BLOSUM62, gaps, range, xdrop);
    }

    AlignResult resRev = block_res_aa_trace_xdrop(blockRev);

    unsigned int qStartPos = querySeqLen - (qStartRev + resRev.query_idx);
    unsigned int qEndPosAlign = querySeqLen;
    char *querySeqAlign = querySeq + qStartPos;
    size_t len_querySeqAlign = qEndPosAlign - qStartPos;
    PaddedBytes *queryPadded = block_new_padded_aa(len_querySeqAlign, range.max);
    block_set_bytes_padded_aa(queryPadded, (const uint8_t *) querySeqAlign, len_querySeqAlign, range.max);

    unsigned int tStartPos = targetSeqObj->L - (tStartRev + resRev.reference_idx);
    unsigned int tEndPosAlign = targetSeqObj->L;
    char *targetSeqAlign = targetSeq + tStartPos;
    size_t len_targetSeqAlign = tEndPosAlign - tStartPos;

    AAProfile *targetProfile = nullptr;
    PaddedBytes *targetPadded = nullptr;
    if (useProfile) {
        targetProfile = block_new_aaprofile(len_targetSeqAlign, range.max, gaps.extend);
        initializeProfile(rawProfileMatrix, tStartPos, tEndPosAlign, targetSeqObj->L, targetProfile, false);
        for (size_t i = 0; i < len_targetSeqAlign; i++) {
            block_set_gap_open_C_aaprofile(targetProfile, i, gaps.open);
            block_set_gap_close_C_aaprofile(targetProfile, i, 0);
            block_set_gap_open_R_aaprofile(targetProfile, i, gaps.open);
        }
    } else {
        targetPadded = block_new_padded_aa(len_targetSeqAlign, range.max);
        block_set_bytes_padded_aa(targetPadded, (const uint8_t *) targetSeqAlign, len_targetSeqAlign, range.max);
    }

    if (useProfile) {
        block_align_profile_aa_trace_xdrop(block, queryPadded, targetProfile, range, xdrop);
    } else {
        block_align_aa_trace_xdrop(block, queryPadded, targetPadded, &BLOSUM62, gaps, range, xdrop);
    }

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
                        if (targetSeqRevAlign[targetPos] == 'X' || querySeqRevAlign[queryPos] == 'X') {
                            aaIds++;
                        } else if (targetSeqRevAlign[targetPos] == querySeqRevAlign[queryPos]) {
                            aaIds++;
                        }
                    } else {
                        if (targetSeqAlign[targetPos] == 'X' || querySeqAlign[queryPos] == 'X') {
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

    float seqId = Util::computeSeqId(seqIdMode, aaIds, querySeqLen, targetSeqObj->L, alnLength) * (useProfile ? 10.0 : 1.0);
    if (reverseCigar) {
        std::reverse(backtrace.begin(), backtrace.end());
    }

    Matcher::result_t  realResult = Matcher::result_t(targetSeqObj->getDbKey(), bitScore, qcov, dbcov, seqId, evalue, alnLength,
                                   qStartPos, qEndPos, querySeqLen, tStartPos, dbEndPos, targetSeqObj->L,
                                   backtrace);

    block_free_padded_aa(queryRevPadded);
    block_free_padded_aa(queryPadded);
    if (targetRevPadded != nullptr) {
        block_free_padded_aa(targetRevPadded);
    }
    if (targetPadded != nullptr) {
        block_free_padded_aa(targetPadded);
    }
    if (targetRevProfile != nullptr) {
        block_free_aaprofile(targetRevProfile);
    }
    if (targetProfile != nullptr) {
        block_free_aaprofile(targetProfile);
    }
    block_free_cigar(cigar);

    return realResult;
}

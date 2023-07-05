// Created by Match on 8/13/2021.
#include "DistanceCalculator.h"
#include "BlockAligner.h"
#include "Sequence.h"
#include "block_aligner.h"
#include "Util.h"
#include "Parameters.h"
#include "SubstitutionMatrix.h"
#include "SRAUtil.h"

BlockAligner::BlockAligner(
    size_t maxSequenceLength,
    uintptr_t min = 32, uintptr_t max = 32,
    int8_t gapOpen = -11, int8_t gapExtend = -1 //,
    // SubstitutionMatrix& subMat
) : range({min, max}),
    gaps({gapOpen, gapExtend}) {
    a = block_new_padded_aa(maxSequenceLength, max);
    b = block_new_padded_aa(maxSequenceLength, max);
    // aProfile = block_new_aaprofile(maxSequenceLength, max, gaps.extend);
    // bProfile = block_new_aaprofile(maxSequenceLength, max, gaps.extend);
    blockTrace = block_new_aa_trace_xdrop(maxSequenceLength, maxSequenceLength, max);
    blockNoTrace = block_new_aa_xdrop(maxSequenceLength, maxSequenceLength, max);
    cigar = block_new_cigar(maxSequenceLength, maxSequenceLength);

    // matrix = block_new_simple_aamatrix(1, -1);
    // for (int i = 0; i < subMat.alphabetSize; i++) {
    //     for (int j = 0; j < subMat.alphabetSize; j++) {
    //         block_set_simple_aamatrix(
    //             matrix,
    //             subMat.num2aa[i],
    //             subMat.num2aa[j],
    //             subMat.subMatrix[i][j]
    //         );
    //     }
    // }

    querySeqRev = static_cast<char *>(malloc(maxSequenceLength * sizeof(char)));
    targetSeqRev = static_cast<char *>(malloc(maxSequenceLength * sizeof(char)));
}

BlockAligner::~BlockAligner() {
    block_free_cigar(cigar);
    block_free_aa_trace_xdrop(blockTrace);
    block_free_aa_xdrop(blockNoTrace);
    block_free_padded_aa(a);
    block_free_padded_aa(b);
    // block_free_aaprofile(aProfile);
    // block_free_aaprofile(bProfile);
    // block_free_simple_aamatrix(matrix);

    free(querySeqRev);
    free(targetSeqRev);
}

void BlockAligner::initTarget(Sequence &target) {
    // SRAUtil::stripInvalidChars(query->getSeqData(), querySeq);
    // memcpy(querySeq, query->getSeqData(), query->L);
    // querySeqLen = query->L;
    targetSeq = target.getSeqData();
    targetLength = target.L;
    SRAUtil::strrev(targetSeqRev, targetSeq, targetLength);
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

typedef struct LocalAln {
    size_t a_start;
    size_t b_start;
    size_t a_end;
    size_t b_end;
    int32_t score;
} LocalAln;

// note: traceback cigar string will be reversed, but LocalAln will contain correct start and end positions
LocalAln align_local(BlockHandle block_trace, BlockHandle block_no_trace, size_t a_len, const char* a_str, const char* a_rev, PaddedBytes* a, size_t b_len, const char* b_str, const char* b_rev, PaddedBytes* b, const AAMatrix* matrix, Gaps gaps, size_t a_idx, size_t b_idx, Cigar* cigar, SizeRange range) {
    LocalAln res_aln;
    AlignResult res;

    int32_t x_drop = -(range.min * gaps.extend + gaps.open);

    // forwards alignment starting at (a_idx, b_idx)
    block_set_bytes_padded_aa(a, (uint8_t*)(a_str + a_idx), a_len - a_idx, range.max);
    block_set_bytes_padded_aa(b, (uint8_t*)(b_str + b_idx), b_len - b_idx, range.max);

    block_align_aa_xdrop(block_no_trace, a, b, matrix, gaps, range, x_drop);
    res = block_res_aa_xdrop(block_no_trace);

    res_aln.a_end = a_idx + res.query_idx;
    res_aln.b_end = b_idx + res.reference_idx;

    // reversed alignment starting at the max score location from forwards alignment
    a_idx = a_len - (a_idx + res.query_idx);
    b_idx = b_len - (b_idx + res.reference_idx);

    block_set_bytes_padded_aa(a, (uint8_t*)(a_rev + a_idx), a_len - a_idx, range.max);
    block_set_bytes_padded_aa(b, (uint8_t*)(b_rev + b_idx), b_len - b_idx, range.max);

    // start at a reasonable min_size based on the forwards alignment
    // min_size >>= ITER;

    block_align_aa_trace_xdrop(block_trace, a, b, matrix, gaps, range, x_drop);
    res = block_res_aa_trace_xdrop(block_trace);
    block_cigar_aa_trace_xdrop(block_trace, res.query_idx, res.reference_idx, cigar);

    res_aln.a_start = a_len - (a_idx + res.query_idx);
    res_aln.b_start = b_len - (b_idx + res.reference_idx);
    res_aln.score = res.score;
    return res_aln;
}


Matcher::result_t
BlockAligner::align(Sequence &query,
                    DistanceCalculator::LocalAlignment alignment,
                    EvalueComputation *evaluer,
                    int xdrop,
                    BaseMatrix *subMat,
                    bool useProfile) {

    // const int8_t *rawProfileMatrix = nullptr;
        // rawProfileMatrix = targetSeqObj->getAlignmentProfile();
        // SRAUtil::stripInvalidChars(
        //         SRAUtil::extractProfileSequence(targetSeqObj->getSeqData(), targetSeqObj->L, subMat).c_str(),
        //         targetSeq
        // );


    // AAProfile *targetRevProfile = nullptr;
    // PaddedBytes *targetRevPadded = nullptr;
    // if (useProfile) {
    //     targetRevProfile = block_new_aaprofile(len_targetSeqRevAlign, range.max, gaps.extend);
    //     initializeProfile(rawProfileMatrix, tStartRev, tEndRev, targetSeqObj->L, targetRevProfile, true);
    //     for (size_t i = 0; i < len_targetSeqRevAlign; i++) {
    //         block_set_gap_open_C_aaprofile(targetRevProfile, i, gaps.open);
    //         block_set_gap_close_C_aaprofile(targetRevProfile, i, 0);
    //         block_set_gap_open_R_aaprofile(targetRevProfile, i, gaps.open);
    //     }

    // AAProfile *targetProfile = nullptr;
    // PaddedBytes *targetPadded = nullptr;
    // if (useProfile) {
    //     targetProfile = block_new_aaprofile(len_targetSeqAlign, range.max, gaps.extend);
    //     initializeProfile(rawProfileMatrix, tStartPos, tEndPosAlign, targetSeqObj->L, targetProfile, false);
    //     for (size_t i = 0; i < len_targetSeqAlign; i++) {
    //         block_set_gap_open_C_aaprofile(targetProfile, i, gaps.open);
    //         block_set_gap_close_C_aaprofile(targetProfile, i, 0);
    //         block_set_gap_open_R_aaprofile(targetProfile, i, gaps.open);
    //     }
        // block_align_profile_aa_trace_xdrop(block, queryPadded, targetProfile, range, xdrop);

    unsigned int qKey = query.getDbKey();
    const char* b_str = query.getSeqData();
    size_t b_len = query.L;
    SRAUtil::strrev(querySeqRev, b_str, b_len);
    const char* b_rev = querySeqRev;

    const char* a_str = targetSeq;
    size_t a_len = targetLength;
    const char* a_rev = targetSeqRev;

    // unsigned int qUngappedStartPos = alignment.startPos + ((alignment.diagonal < 0) ? alignment.distToDiagonal : 0);
    unsigned int qUngappedEndPos = alignment.endPos + ((alignment.diagonal < 0) ? alignment.distToDiagonal : 0) - 1;
    // unsigned int dbUngappedStartPos = alignment.startPos + ((alignment.diagonal < 0) ? 0 : alignment.distToDiagonal);
    unsigned int dbUngappedEndPos = alignment.endPos + ((alignment.diagonal < 0) ? 0 : alignment.distToDiagonal) - 1;

    // Debug(Debug::INFO) << "qUnGapStart: " << qUngappedStartPos << " qUnGapEnd: " << qUngappedEndPos << " qlen: " << a_len << "\n";
    // Debug(Debug::INFO) << "dbUnGapStart: " << dbUngappedStartPos << " dbUnGapEnd: " << dbUngappedEndPos << " dblen: " << b_len << "\n";
    // Debug(Debug::INFO) << std::string(a_str + qUngappedStartPos, qUngappedEndPos - qUngappedStartPos) << "\n";
    // Debug(Debug::INFO) << std::string(b_str + dbUngappedStartPos, dbUngappedEndPos - dbUngappedStartPos) << "\n";


    LocalAln local_aln = align_local(blockTrace, blockNoTrace, a_len, a_str, a_rev, a, b_len, b_str, b_rev, b, &BLOSUM62, gaps, qUngappedEndPos, dbUngappedEndPos, cigar, range);
    // printf("a: %s\nb: %s\nscore: %d\nstart idx: (%lu, %lu)\nend idx: (%lu, %lu)\n",
    //     a_str,
    //     b_str,
    //     local_aln.score,
    //     local_aln.a_start,
    //     local_aln.b_start,
    //     local_aln.a_end,
    //     local_aln.b_end
    // );

    float qcov = SmithWaterman::computeCov(local_aln.a_start, local_aln.a_end, a_len);
    float dbcov = SmithWaterman::computeCov(local_aln.b_start, local_aln.b_end, b_len);
    
    int bitScore = static_cast<int>(evaluer->computeBitScore(local_aln.score) + 0.5);
    double evalue = evaluer->computeEvalue(local_aln.score, a_len);

    // Note: 'M' signals either a match or mismatch
    char ops_char[] = {' ', 'M', '=', 'X', 'I', 'D'};

    int alnLength = Matcher::computeAlnLength(local_aln.a_start, local_aln.a_end, local_aln.b_start, local_aln.b_end);

    std::string backtrace;
    size_t cigarLength = block_len_cigar(cigar);
    size_t aaIds = 0;
    if (cigarLength > 0) {
        int32_t targetPos = 0, queryPos = 0;
        for (unsigned long c = 0; c < cigarLength; ++c) {
            OpLen o = block_get_cigar(cigar, cigarLength - 1 - c);
            char letter = ops_char[o.op];
            uint32_t length = o.len;

            switch (letter) {
                case '=':
                    aaIds += length;
                    // FALLTHROUGH
                case 'X':
                    // FALLTHROUGH
                case 'M':
                    queryPos += length;
                    targetPos += length;
                    backtrace.append('M', length);
                    break;
                case 'I':
                    queryPos += length;
                    backtrace.append('I', length);
                    break;
                case 'D':
                    targetPos += length;
                    backtrace.append('D', length);
                    break;
            }
        }
        alnLength = backtrace.length();
    }

    // Debug(Debug::INFO) << "Alignment length: " << alnLength << "\n";
    // Debug(Debug::INFO) << "aaIds: " << aaIds << "\n";
    // Debug(Debug::INFO) << "a_len: " << a_len << "\n";
    // Debug(Debug::INFO) << "b_len: " << b_len << "\n";

    int seqIdMode = Parameters::SEQ_ID_ALN_LEN;
    float seqId = Util::computeSeqId(seqIdMode, aaIds, a_len, b_len, alnLength) * (useProfile ? 10.0 : 1.0);
        // Debug(Debug::INFO) << "seqId: " << seqId << "\n";

        // EXIT(EXIT_FAILURE);

    return Matcher::result_t(
        qKey, bitScore, qcov, dbcov, seqId, evalue, alnLength,
        local_aln.a_start, local_aln.a_end, a_len, local_aln.b_start, local_aln.b_end, b_len,
        backtrace
    );
}

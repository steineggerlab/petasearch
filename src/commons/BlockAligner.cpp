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
    uintptr_t min,
    uintptr_t max,
    int8_t gapOpen,
    int8_t gapExtend,
    int dbtype
) : range({min, max}),
    gaps({gapOpen, gapExtend}),
    dbtype(dbtype) {
    a = block_new_padded_aa(maxSequenceLength, max);
    if (dbtype == Parameters::DBTYPE_AMINO_ACIDS) {
        b = block_new_padded_aa(maxSequenceLength, max);
        matrix = block_new_simple_aamatrix(1, -1);
    } else {
        bProfile = block_new_aaprofile(maxSequenceLength, max, gaps.extend);
    }
    blockTrace = block_new_aa_trace_xdrop(maxSequenceLength, maxSequenceLength, max);
    blockNoTrace = block_new_aa_xdrop(maxSequenceLength, maxSequenceLength, max);
    cigar = block_new_cigar(maxSequenceLength, maxSequenceLength);
}

BlockAligner::~BlockAligner() {
    block_free_cigar(cigar);
    block_free_aa_trace_xdrop(blockTrace);
    block_free_aa_xdrop(blockNoTrace);
    block_free_padded_aa(a);
    if (dbtype == Parameters::DBTYPE_AMINO_ACIDS) {
        block_free_padded_aa(b);
        block_free_aamatrix(matrix);
    } else {
        block_free_aaprofile(bProfile);
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
LocalAln align_local(
    BlockHandle block_trace, BlockHandle block_no_trace,
    const char* a_str, size_t a_len, PaddedBytes* a,
    const char* b_str, size_t b_len, PaddedBytes* b,
    const AAMatrix* matrix, Gaps gaps,
    size_t a_idx, size_t b_idx,
    Cigar* cigar, SizeRange range, int32_t x_drop
) {
    LocalAln res_aln;
    AlignResult res;

    // forwards alignment starting at (a_idx, b_idx)
    block_set_bytes_padded_aa(a, (uint8_t*)(a_str + a_idx), a_len - a_idx, range.max);
    block_set_bytes_padded_aa(b, (uint8_t*)(b_str + b_idx), b_len - b_idx, range.max);

    block_align_aa_xdrop(block_no_trace, a, b, matrix, gaps, range, x_drop);
    res = block_res_aa_xdrop(block_no_trace);

    res_aln.a_end = a_idx + res.query_idx;
    res_aln.b_end = b_idx + res.reference_idx;

    // reversed alignment starting at the max score location from forwards alignment
    block_set_bytes_rev_padded_aa(a, (uint8_t*)a_str, res_aln.a_end, range.max);
    block_set_bytes_rev_padded_aa(b, (uint8_t*)b_str, res_aln.b_end, range.max);

    block_align_aa_trace_xdrop(block_trace, a, b, matrix, gaps, range, x_drop);
    res = block_res_aa_trace_xdrop(block_trace);
    block_cigar_eq_aa_trace_xdrop(block_trace, a, b, res.query_idx, res.reference_idx, cigar);

    res_aln.a_start = res_aln.a_end - res.query_idx;
    res_aln.b_start = res_aln.b_end - res.reference_idx;
    res_aln.score = res.score;
    return res_aln;
}

LocalAln align_local_profile(
    BlockHandle block_trace, BlockHandle block_no_trace,
    const char* a_str, const size_t a_len, PaddedBytes* a,
    const char* b_str, const size_t b_len, AAProfile* bProfile,
    Gaps gaps, BaseMatrix& subMat,
    size_t a_idx, size_t b_idx,
    Cigar* cigar, SizeRange range, int32_t x_drop
) {
    LocalAln res_aln;
    AlignResult res;

    // forwards alignment starting at (a_idx, b_idx)
    block_set_bytes_padded_aa(a, (uint8_t*)(a_str + a_idx), a_len - a_idx, range.max);

    int aa = Sequence::PROFILE_AA_SIZE;
    block_clear_aaprofile(bProfile, b_len / aa - b_idx, range.max);
    block_set_all_aaprofile(bProfile, (uint8_t*)subMat.num2aa, aa, (int8_t*)(b_str + b_idx * aa), b_len - b_idx * aa);
    block_set_all_gap_open_C_aaprofile(bProfile, gaps.open);
    block_set_all_gap_close_C_aaprofile(bProfile, 0);
    block_set_all_gap_open_R_aaprofile(bProfile, gaps.open);

    block_align_profile_aa_xdrop(block_no_trace, a, bProfile, range, x_drop);
    res = block_res_aa_xdrop(block_no_trace);

    res_aln.a_end = a_idx + res.query_idx;
    res_aln.b_end = b_idx + res.reference_idx;

    // reversed alignment starting at the max score location from forwards alignment
    block_set_bytes_rev_padded_aa(a, (uint8_t*)a_str, res_aln.a_end, range.max);

    block_clear_aaprofile(bProfile, res_aln.b_end, range.max);
    block_set_all_rev_aaprofile(bProfile, (uint8_t*)subMat.num2aa, aa, (int8_t*)b_str, res_aln.b_end * aa);
    block_set_all_gap_open_C_aaprofile(bProfile, gaps.open);
    block_set_all_gap_close_C_aaprofile(bProfile, 0);
    block_set_all_gap_open_R_aaprofile(bProfile, gaps.open);

    block_align_profile_aa_trace_xdrop(block_trace, a, bProfile, range, x_drop);
    res = block_res_aa_trace_xdrop(block_trace);
    block_cigar_aa_trace_xdrop(block_trace, res.query_idx, res.reference_idx, cigar);

    res_aln.a_start = res_aln.a_end - res.query_idx;
    res_aln.b_start = res_aln.b_end - res.reference_idx;
    res_aln.score = res.score;
    return res_aln;
}


Matcher::result_t
BlockAligner::align(
    const char* targetSeq,
    unsigned int targetLength,
    const char* querySeq,
    unsigned int queryLength,
    DistanceCalculator::LocalAlignment alignment,
    EvalueComputation *evaluer,
    int xdrop,
    BaseMatrix& subMat
) {
    if (dbtype == Parameters::DBTYPE_AMINO_ACIDS) {
        for (int i = 0; i < subMat.alphabetSize; i++) {
            for (int j = 0; j < subMat.alphabetSize; j++) {
                block_set_aamatrix(
                    matrix,
                    subMat.num2aa[i],
                    subMat.num2aa[j],
                    subMat.subMatrix[i][j]
                );
            }
        }
    }

    unsigned int qUngappedStartPos = alignment.startPos + ((alignment.diagonal < 0) ? alignment.distToDiagonal : 0);
    unsigned int qUngappedEndPos = alignment.endPos + ((alignment.diagonal < 0) ? alignment.distToDiagonal : 0);
    unsigned int dbUngappedStartPos = alignment.startPos + ((alignment.diagonal < 0) ? 0 : alignment.distToDiagonal);
    unsigned int dbUngappedEndPos = alignment.endPos + ((alignment.diagonal < 0) ? 0 : alignment.distToDiagonal);
    qUngappedEndPos = std::max(1u, qUngappedEndPos) - 1;
    dbUngappedEndPos = std::max(1u, dbUngappedEndPos) - 1;
    qUngappedStartPos = (qUngappedStartPos > qUngappedEndPos) ? qUngappedEndPos : qUngappedStartPos;
    dbUngappedStartPos = (dbUngappedStartPos > dbUngappedEndPos) ? dbUngappedEndPos : dbUngappedStartPos;
    qUngappedEndPos = (qUngappedEndPos < qUngappedStartPos) ? qUngappedStartPos : qUngappedEndPos;
    dbUngappedEndPos = (dbUngappedEndPos < dbUngappedStartPos) ? dbUngappedStartPos : dbUngappedEndPos;

    LocalAln local_aln;
    if (dbtype == Parameters::DBTYPE_AMINO_ACIDS) {
        local_aln = align_local(blockTrace, blockNoTrace, targetSeq, targetLength, a, querySeq, queryLength, b, matrix, gaps, qUngappedEndPos, dbUngappedEndPos, cigar, range, xdrop);
    } else {
        local_aln = align_local_profile(blockTrace, blockNoTrace, targetSeq, targetLength, a, querySeq, queryLength, bProfile, gaps, subMat, qUngappedEndPos, dbUngappedEndPos, cigar, range, xdrop);
    }

    float qcov = SmithWaterman::computeCov(local_aln.a_start, local_aln.a_end, targetLength);
    float dbcov = SmithWaterman::computeCov(local_aln.b_start, local_aln.b_end, queryLength);
    
    int bitScore = static_cast<int>(evaluer->computeBitScore(local_aln.score) + 0.5);
    double evalue = evaluer->computeEvalue(local_aln.score, targetLength);

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
                    backtrace.append(length, 'M');
                    break;
                case 'I':
                    queryPos += length;
                    backtrace.append(length, 'I');
                    break;
                case 'D':
                    targetPos += length;
                    backtrace.append(length, 'D');
                    break;
            }
        }
        alnLength = backtrace.length();
    }

    int seqIdMode = Parameters::SEQ_ID_ALN_LEN;
    float seqId = Util::computeSeqId(seqIdMode, aaIds, targetLength, queryLength, alnLength); // * (useProfile ? 10.0 : 1.0);

    return Matcher::result_t(
        0, bitScore, qcov, dbcov, seqId, evalue, alnLength,
        local_aln.a_start, local_aln.a_end, targetLength, local_aln.b_start, local_aln.b_end, queryLength,
        backtrace
    );
}

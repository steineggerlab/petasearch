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
    // SubstitutionMatrix& subMat
) : range({min, max}),
    gaps({gapOpen, gapExtend}),
    dbtype(dbtype) {
    a = block_new_padded_aa(maxSequenceLength, max);
    if (dbtype == Parameters::DBTYPE_AMINO_ACIDS) {
        b = block_new_padded_aa(maxSequenceLength, max);
    } else {
        // bProfile = block_new_aaprofile(maxSequenceLength, max, gaps.extend);
        // volatile int val = gaps.open;
        // for (size_t i = 0; i < maxSequenceLength; i++) {
        //     block_set_gap_open_C_aaprofile(bProfile, i, val);
        //     block_set_gap_close_C_aaprofile(bProfile, i, 0);
        //     block_set_gap_open_R_aaprofile(bProfile, i, val);
        // }
    }
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
    if (dbtype == Parameters::DBTYPE_AMINO_ACIDS) {
        block_free_padded_aa(b);
    // } else {
        // block_free_aaprofile(bProfile);
    }
    // block_free_simple_aamatrix(matrix);

    free(querySeqRev);
    free(targetSeqRev);
}

void BlockAligner::initTarget(Sequence &target) {
    targetSeq = target.getSeqData();
    targetLength = target.L;
    SRAUtil::strrev(targetSeqRev, targetSeq, targetLength);
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
    size_t a_len, const char* a_str, const char* a_rev, PaddedBytes* a,
    size_t b_len, const char* b_str, const char* b_rev, PaddedBytes* b,
    const AAMatrix* matrix, Gaps gaps,
    size_t a_idx, size_t b_idx,
    Cigar* cigar, SizeRange range, int32_t x_drop
) {
    LocalAln res_aln;
    AlignResult res;

    // int32_t x_drop = -(range.min * gaps.extend + gaps.open);

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

LocalAln align_local_profile(
    BlockHandle block_trace, BlockHandle block_no_trace,
    const size_t a_len, const char* a_str, const char* a_rev, PaddedBytes* a,
    const size_t b_len, const char* b_str, AAProfile* b,
    Gaps gaps, BaseMatrix& subMat,
    size_t a_idx, size_t b_idx,
    Cigar* cigar, SizeRange range, int32_t x_drop
) {
    LocalAln res_aln;
    AlignResult res;

    // int32_t x_drop = 10000; // -(range.min * gaps.extend + gaps.open);
    // forwards alignment starting at (a_idx, b_idx)
    block_set_bytes_padded_aa(a, (uint8_t*)(a_str + a_idx), a_len - a_idx, range.max);

    // int curr_max = INT_MAX;
    // block_set_len_aaprofile(b, b_len * 10);
    b = block_new_aaprofile(b_len - b_idx, range.max, gaps.extend);
    for (size_t i = 0; i < (b_len - b_idx); ++i) {
        // printf("%d ", b_idx + i);
        for (uint8_t j = 0; j < Sequence::PROFILE_AA_SIZE; ++j) {
            volatile int val = static_cast<int8_t>(static_cast<short>(b_str[(b_idx + i) * Sequence::PROFILE_READIN_SIZE + j]) / 4);
            // if (val < curr_max) {
            //     curr_max = val;
            // }
            block_set_aaprofile(b, i + 1, subMat.num2aa[j], val);
            // block_set_aaprofile(b, i + 1, subMat.num2aa[j], 2);
            // printf("%d ", b_profile_matrix[(j * b_len) + (b_idx + i)]);
        }
        // printf("\n");
    }
    volatile int val = gaps.open;
    for (size_t i = 0; i < b_len - b_idx; i++) {
        block_set_gap_open_C_aaprofile(b, i, val);
        block_set_gap_close_C_aaprofile(b, i, 0);
        block_set_gap_open_R_aaprofile(b, i, val);
    }
    // block_clear_aaprofile(b, b_len);
    // for (size_t i = 0; i < (b_len - b_idx); i++) {
    //     block_set_gap_open_C_aaprofile(b, i, gaps.open);
    //     block_set_gap_close_C_aaprofile(b, i, 0);
    //     block_set_gap_open_R_aaprofile(b, i, gaps.open);
    // }

    // std::cout << "curr_max: " << curr_max << "\n";

    block_align_profile_aa_xdrop(block_no_trace, a, b, range, x_drop);
    res = block_res_aa_xdrop(block_no_trace);

    // std::cout << res.score << std::endl;
    // EXIT(EXIT_FAILURE);

    res_aln.a_end = a_idx + res.query_idx;
    res_aln.b_end = b_idx + res.reference_idx;

    a_idx = a_len - (a_idx + res.query_idx);
    b_idx = b_len - (b_idx + res.reference_idx);

    // std::cout << "res.query_idx: " << res.query_idx << std::endl;
    // std::cout << "res.reference_idx: " << res.reference_idx << std::endl;
    // std::cout << "a_idx: " << a_idx << std::endl;
    // std::cout << "b_idx: " << b_idx << std::endl;

    block_set_bytes_padded_aa(a, (uint8_t*)(a_rev + a_idx), a_len - a_idx, range.max);
    // block_set_len_aaprofile(b, b_len * 10);
    block_free_aaprofile(b);
    b = block_new_aaprofile(b_len - b_idx, range.max, gaps.extend);
    for (size_t i = 0; i < (b_len - b_idx); ++i) {
        // printf("%d ", b_len - b_idx - i);
        for (uint8_t j = 0; j < Sequence::PROFILE_AA_SIZE; ++j) {
            // volatile int val = b_profile_matrix[(j * b_len) + (b_len - b_idx - i)];
            volatile int val = static_cast<int8_t>(static_cast<short>(b_str[(b_len - b_idx - i) * Sequence::PROFILE_READIN_SIZE + j]) / 4);
            // if (val < curr_max) {
            //     curr_max = val;
            // }
            block_set_aaprofile(b, i + 1, subMat.num2aa[j], val);
            // block_set_aaprofile(b, i, subMat.num2aa[j], 2);
            // printf("%d ", b_profile_matrix[(j * b_len) + (b_len - b_idx - i)]);
        }
        // printf("\n");
    }
    for (size_t i = 0; i < b_len - b_idx; i++) {
        block_set_gap_open_C_aaprofile(b, i, val);
        block_set_gap_close_C_aaprofile(b, i, 0);
        block_set_gap_open_R_aaprofile(b, i, val);
    }
        // std::cout << "curr_max: " << curr_max << "\n";


    // start at a reasonable min_size based on the forwards alignment
    // min_size >>= ITER;

    block_align_profile_aa_trace_xdrop(block_trace, a, b, range, x_drop);
    res = block_res_aa_trace_xdrop(block_trace);
    block_cigar_aa_trace_xdrop(block_trace, res.query_idx, res.reference_idx, cigar);

    // std::cout << "score: " << res.score << std::endl;

    res_aln.a_start = a_len - (a_idx + res.query_idx);
    res_aln.b_start = b_len - (b_idx + res.reference_idx);
    res_aln.score = res.score;

    block_free_aaprofile(b);

    return res_aln;
}


Matcher::result_t
BlockAligner::align(
    const char* querySeq,
    unsigned int queryLength,
    DistanceCalculator::LocalAlignment alignment,
    EvalueComputation *evaluer,
    int xdrop,
    BaseMatrix* subMat
) {
    size_t b_len = queryLength;
    const char* b_str = querySeq;
    if (dbtype == Parameters::DBTYPE_AMINO_ACIDS) {
        SRAUtil::strrev(querySeqRev, b_str, b_len);
    }
    const char* b_rev = querySeqRev;

    size_t a_len = targetLength;
    const char* a_str = targetSeq;
    const char* a_rev = targetSeqRev;

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

    // Debug(Debug::INFO) << "startPos: " << alignment.startPos << " endPos: " << alignment.endPos << " diagonal: " << alignment.diagonal << " distToDiagonal: " << alignment.distToDiagonal << "\n";
    // Debug(Debug::INFO) << "qUnGapStart: " << qUngappedStartPos << " qUnGapEnd: " << qUngappedEndPos << " qlen: " << a_len << "\n";
    // Debug(Debug::INFO) << "dbUnGapStart: " << dbUngappedStartPos << " dbUnGapEnd: " << dbUngappedEndPos << " dblen: " << b_len << "\n";
    // Debug(Debug::INFO) << std::string(a_str + qUngappedStartPos, qUngappedEndPos - qUngappedStartPos) << "\n";
    // if (useProfile) {
    //     std::string realSeq = SRAUtil::extractProfileSequence(b_str, b_len - 1, subMat);
    //     Debug(Debug::INFO) << std::string(realSeq.c_str() + dbUngappedStartPos, dbUngappedEndPos - dbUngappedStartPos) << "\n";
    // } else {
    //     Debug(Debug::INFO) << std::string(b_str + dbUngappedStartPos, dbUngappedEndPos - dbUngappedStartPos) << "\n";
    // }

    LocalAln local_aln;
    if (dbtype == Parameters::DBTYPE_AMINO_ACIDS) {
        local_aln = align_local(blockTrace, blockNoTrace, a_len, a_str, a_rev, a, b_len, b_str, b_rev, b, &BLOSUM62, gaps, qUngappedEndPos, dbUngappedEndPos, cigar, range, xdrop);
    } else {
        local_aln = align_local_profile(blockTrace, blockNoTrace, a_len, a_str, a_rev, a, b_len, b_str, bProfile, gaps, *subMat, qUngappedEndPos, dbUngappedEndPos, cigar, range, xdrop);
    }
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

    // Debug(Debug::INFO) << "Alignment length: " << alnLength << "\n";
    // Debug(Debug::INFO) << "aaIds: " << aaIds << "\n";
    // Debug(Debug::INFO) << "a_len: " << a_len << "\n";
    // Debug(Debug::INFO) << "b_len: " << b_len << "\n";

    int seqIdMode = Parameters::SEQ_ID_ALN_LEN;
    float seqId = Util::computeSeqId(seqIdMode, aaIds, a_len, b_len, alnLength); // * (useProfile ? 10.0 : 1.0);
        // Debug(Debug::INFO) << "seqId: " << seqId << "\n";

        // EXIT(EXIT_FAILURE);

    return Matcher::result_t(
        0, bitScore, qcov, dbcov, seqId, evalue, alnLength,
        local_aln.a_start, local_aln.a_end, a_len, local_aln.b_start, local_aln.b_end, b_len,
        backtrace
    );
}

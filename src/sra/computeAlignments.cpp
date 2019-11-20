#include "LocalParameters.h"
#include "Command.h"
#include "Debug.h"
#include "QueryTableEntry.h"
#include "IndexReader.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "NucleotideMatrix.h"
#include "EvalueComputation.h"
#include "DistanceCalculator.h"
#include "Matcher.h"

#include "omptl/omptl_algorithm"

#ifdef OPENMP
#include <omp.h>
#endif

int blockByDiagSort(const QueryTableEntry &first, const QueryTableEntry &second);

int computeAlignments(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.spacedKmer = false;
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    par.evalThr = 10000;
    par.parseParameters(argc, argv, command, true, 1, 0);

    // query target prev_result new_result is the order

    DBReader<unsigned int> querySequenceReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    querySequenceReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> targetSequenceReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    targetSequenceReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_WRITABLE);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();

    BaseMatrix *subMat;
    int seqType = targetSequenceReader.getDbtype();
    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.nucleotides, 1.0, 0.0);
    }
    else {
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.aminoacids, 2.0, 0.0);
    }
    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(*subMat);
    EvalueComputation evaluer(targetSequenceReader.getAminoAcidDBSize(), subMat);

    Debug::Progress progress(resultReader.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Sequence querySeq(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false);
        Sequence targetSeq(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false);
        // Sequence nextSequence(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false);

        Indexer idx(subMat->alphabetSize - 1, par.kmerSize);

        Matcher matcher(seqType, par.maxSeqLen, subMat, &evaluer, par.compBiasCorrection, par.gapOpen, par.gapExtend);

        char buffer[1024];
        std::vector<Matcher::result_t> results;
        results.reserve(300);

        std::string data;
        data.reserve(1000);

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < resultReader.getSize(); ++i) {
            progress.updateProgress();

            size_t targetId = resultReader.getDbKey(i);
            unsigned int targetKey = targetSequenceReader.getDbKey(targetId);

            //determine block size of current target sequence
//            size_t targetId = targetSequenceReader.getId(targetKey);
            const char *targetSeqData = targetSequenceReader.getData(targetId, thread_idx);
            unsigned int targetSeqLen = targetSequenceReader.getSeqLen(targetId);
            targetSeq.mapSequence(targetId, targetKey, targetSeqData, targetSeqLen);

            // TODO: this might be wasted if no single hit hit the target
            matcher.initQuery(&targetSeq);

            // prefetch next sequence
            // nextTargetID = (++currentPos)->targetSequenceID;
            // char *nextData = targetSequenceReader.getData(nextTargetID, 0);
            // TODO: prefetch instruction
            // nextSequence.mapSequence(nextTargetID,0,nextData);

            QueryTableEntry *targetBlock = (QueryTableEntry *)resultReader.getData(i, thread_idx);
            size_t targetBlockCount = resultReader.getEntryLen(i) / sizeof(QueryTableEntry);
            //calculate diags for each query block in the target block
            for (size_t j = 0; j < targetBlockCount;) {
                QueryTableEntry *queryBlockStart = targetBlock + j;
                QueryTableEntry *queryBlockEnd = queryBlockStart;
                unsigned int lastQueryId;
                do {
                    lastQueryId = queryBlockEnd->querySequenceId;
//                    idx.printKmer(queryBlockEnd->Query.kmer, par.kmerSize, subMat->int2aa);
//                    Debug(Debug::ERROR) << "\n";
                    while (targetSeq.hasNextKmer()) {
                        const int *kmer = targetSeq.nextKmer();
//                        idx.printKmer(idx.int2index(kmer, 0, par.kmerSize), par.kmerSize, subMat->int2aa);
//                        Debug(Debug::ERROR) << "\n";
                        if (queryBlockEnd->Query.kmer ==
                            idx.int2index(kmer, 0, par.kmerSize)) { //TODO test performance (cache invalidations)
//                            Debug(Debug::ERROR) << "!!!!!";
                            queryBlockEnd->Result.diag =
                                    queryBlockEnd->Query.kmerPosInQuery - targetSeq.getCurrentPosition();
                            break;
                        }
                    }
                    targetSeq.resetCurrPos();
                    ++queryBlockEnd;
                    ++j;
                } while (j < targetBlockCount && lastQueryId == queryBlockEnd->querySequenceId);

                std::sort(queryBlockStart, queryBlockEnd - 1, blockByDiagSort);

                unsigned int shortestDiagDistance = UINT_MAX;
                for (QueryTableEntry *k = queryBlockStart + 1; k < queryBlockEnd && shortestDiagDistance > 4; ++k) {
                    shortestDiagDistance = std::min(shortestDiagDistance, k->Result.diag - (k - 1)->Result.diag);
                }

                // SUPER IMPORTANT DOCUMENT THIS
                //only compute the alignment if we found at least 2 matches which are close to each other in the sequence --> increases sensitivity
                if (shortestDiagDistance > 4) {
                    continue;
                }

                int maxScore = INT_MIN;
                const char *querySeqData = querySequenceReader.getData(queryBlockStart->querySequenceId, thread_idx);
                unsigned int querySeqLen = querySequenceReader.getSeqLen(queryBlockStart->querySequenceId);
                for (QueryTableEntry *k = queryBlockStart; k < queryBlockEnd; ++k) {
                    // TODO: ??
                    if (k > queryBlockStart && k->Result.diag == (k - 1)->Result.diag) {
                        k->Result.score = (k - 1)->Result.score;
                        continue;
                    }

                    DistanceCalculator::LocalAlignment aln =
                            DistanceCalculator::computeUngappedAlignment(
                                    // TODO check order of Q and T
                                    querySeqData, querySeqLen,
                                    targetSeqData, targetSeqLen,
                                    k->Result.diag, fastMatrix.matrix, par.rescoreMode);
                    k->Result.score = aln.score;


//                    Debug(Debug::ERROR) << querySeqLen << "\t" << std::string(querySeqData, querySeqLen) << "\n";
//                    Debug(Debug::ERROR) << targetSeqLen << "\t" << std::string(k->Result.diag, ' ') << std::string(targetSeqData, targetSeqLen) << "\n";


                    // if (protein) {
                    // different than wiki, explain swap afterwards
                    double eval = evaluer.computeEvalue(aln.score, querySeqLen);
//                    Debug(Debug::ERROR) << k->Result.diag  << "\t" << aln.score <<  "\t" << eval << "\n";

                    // this is bad for nucleotide petasearch, we need to know the best diagonal
                    if (eval <= par.evalThr) {
                        maxScore = aln.score;
                        break;
                    }
                    // } else {
                    // actually calculate maxscore
                    // }
                }

                if (maxScore == INT_MIN) {
                    writer.writeData("", 0, targetKey, thread_idx);
                    continue;
                }

                unsigned int queryKey = querySequenceReader.getDbKey(queryBlockStart->querySequenceId);
                querySeq.mapSequence(queryBlockStart->querySequenceId, queryKey, querySeqData, querySeqLen);

                // we have to swap coverage mode etc
                // replace 0 with actual best diagonal for nucleotide alignments
                Matcher::result_t res = matcher.getSWResult(&querySeq, 0, false, par.covMode, par.covThr, par.evalThr,
                                                            Matcher::SCORE_ONLY, par.seqIdMode, false);
                results.emplace_back(res);
            }
            std::sort(results.begin(), results.end(), Matcher::compareHits);
            for (size_t j = 0; j < results.size(); ++j) {
                size_t len = Matcher::resultToBuffer(buffer, results[j], false, false);
                data.append(buffer, len);
            }
            // we want to write this sorted
            writer.writeData(data.c_str(), data.size(), targetKey, thread_idx);
            results.clear();
            data.clear();
            // after this module in workflow call "swapresults"
            // TODO what about convertalis if we need pairwise alignments or if we need profiles??
        }
    }
    writer.close();

    resultReader.close();
    targetSequenceReader.close();

    return EXIT_SUCCESS;
}

int blockByDiagSort(const QueryTableEntry &first, const QueryTableEntry &second) {
    if (first.Result.diag < second.Result.diag){
        return true;
    }
    if (second.Result.diag < first.Result.diag){
        return false;
    }
    return false;
}
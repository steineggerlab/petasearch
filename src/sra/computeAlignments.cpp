
#include "LocalParameters.h"
#include "Command.h"
#include "Debug.h"
#include "IndexReader.h"
#include "DBReader.h"
#include "SRADBReader.h"
#include "DBWriter.h"
#include "NucleotideMatrix.h"
#include "EvalueComputation.h"
#include "DistanceCalculator.h"
#include "Matcher.h"
#include "QueryTableEntry.h"
#include "BlockAligner.h"

#include "omptl/omptl_algorithm"
#include "SRAUtil.h"

#ifdef OPENMP

#include <omp.h>

#endif

/**
 * Print alignment from bactrace. Debug only
 * */
std::string printAlnFromBt(const char *seq, unsigned int offset, const std::string &bt, bool reverse) {
    std::string out;
    unsigned int seqPos = 0;
    for (uint32_t i = 0; i < bt.size(); ++i) {
        char seqChar = seq[offset + seqPos];
        switch (bt[i]) {
            case 'M':
                out.append(1, seqChar);
                seqPos += 1;
                break;
            case 'I':
                if (reverse) {
                    out.append(1, '-');
                } else {
                    out.append(1, seqChar);
                    seqPos += 1;
                }
                break;
            case 'D':
                if (reverse) {
                    out.append(1, seqChar);
                    seqPos += 1;
                } else {
                    out.append(1, '-');
                }
                break;
        }
    }
    return out;
}

int blockByDiagSort(const QueryTableEntry &first, const QueryTableEntry &second);

bool matcherResultsSort(const Matcher::result_t &first, const Matcher::result_t &second);

const unsigned int INVALID_DIAG = (unsigned int) -1;

bool isWithinNDiagonals(const std::vector<QueryTableEntry> &queries, unsigned int N) {
    unsigned int shortestDiagDistance = UINT_MAX;
    for (size_t i = 1; i < queries.size() && shortestDiagDistance > N; ++i) {
        shortestDiagDistance = std::min(shortestDiagDistance, queries[i].Result.diag - queries[i - 1].Result.diag);
    }
    // SUPER IMPORTANT DOCUMENT THIS
    //only compute the alignment if we found at least 2 matches which are close to each other in the sequence --> increases sensitivity
    return shortestDiagDistance <= N;
}

DistanceCalculator::LocalAlignment ungappedDiagFilter(
        std::vector<QueryTableEntry> &queries,
        const char *querySeqData, size_t querySeqLen,
        const char *targetSeqData, size_t targetSeqLen,
        const char **matrix, EvalueComputation &evaluer, int rescoreMode, double evalThr) {
    int maxScore = INT_MIN;
    DistanceCalculator::LocalAlignment alignmentResult = DistanceCalculator::LocalAlignment();
    unsigned int lastDiagonal = INVALID_DIAG;
    for (size_t i = 1; i < queries.size(); ++i) {
        if (queries[i].Result.diag == lastDiagonal) {
            lastDiagonal = queries[i].Result.diag;
            continue;
        }
        lastDiagonal = queries[i].Result.diag;

        alignmentResult = DistanceCalculator::computeUngappedAlignment(querySeqData, querySeqLen, targetSeqData,
                                                                       targetSeqLen, queries[i].Result.diag,
                                                                       matrix, rescoreMode);
        queries[i].Result.score = alignmentResult.score;

        // if (protein) {
        // different than wiki, explain swap afterwards
        double eval = evaluer.computeEvalue(alignmentResult.score, querySeqLen);
        //                    Debug(Debug::ERROR) << k->Result.diag  << "\t" << aln.score <<  "\t" << eval << "\n";

        // this is bad for nucleotide petasearch, we need to know the best diagonal
        if (eval <= evalThr) {
            maxScore = alignmentResult.score;
            break;
        }
        // } else {
        // actually calculate maxscore
        // }
    }

    if (maxScore == INT_MIN) {
        alignmentResult.diagonal = INVALID_DIAG;
    }
    return alignmentResult;
}

struct Kmer {
    unsigned long kmer = -1;
    int kmerPos = -1;

    Kmer() = default;

    Kmer(unsigned long kmer, int kmerPos) : kmer(kmer), kmerPos(kmerPos) {}
};

bool kmerComparator(const Kmer &kmer1, const Kmer &kmer2) {
    if (kmer1.kmer != kmer2.kmer) {
        return kmer1.kmer < kmer2.kmer;
    }
    return kmer1.kmerPos < kmer2.kmerPos;
}

struct BlockIterator {
    void reset(char *data) {
        buffer = data;
        lastId = (unsigned int) -1;
    }

    bool getNext(std::vector<QueryTableEntry> &block) {
        block.clear();
        if (*buffer == '\0') {
            return false;
        }
        do {
            QueryTableEntry query = QueryTableEntry::parseQueryEntry(buffer);
            if (lastId != query.querySequenceId) {
                lastId = query.querySequenceId;
                if (block.empty() == false) {
                    return true;
                }
            }
            block.emplace_back(query);
            buffer = Util::skipLine(buffer);
            lastId = query.querySequenceId;
        } while (*buffer != '\0');
        return true;
    }

    unsigned int lastId;
    char *buffer;
};

int computeAlignments(int argc, const char **argv, const Command &command) {
    Timer timer;
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.spacedKmer = false;
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    par.evalThr = 1000000;
    par.parseParameters(argc, argv, command, true, 0, 0);

    int inputSeqType = FileUtil::parseDbType(par.db1.c_str());
    bool useProfileSearch = Parameters::isEqualDbtype(inputSeqType, Parameters::DBTYPE_HMM_PROFILE);

    // query target prev_result new_result is the order
    DBReader<unsigned int> querySequenceReader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                               DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    querySequenceReader.open(DBReader<unsigned int>::NOSORT);

    SRADBReader targetSequenceReader(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                     DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    targetSequenceReader.open(DBReader<unsigned int>::NOSORT | DBReader<unsigned int>::LINEAR_ACCCESS);
    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                        DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX |
                                        DBReader<unsigned int>::USE_WRITABLE);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed,
                    Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();

    Debug(Debug::INFO) << "Size of prefilter DB: " << resultReader.getSize() << "\n";

    BaseMatrix *subMat;
    int seqType = targetSequenceReader.getDbtype();
    bool isNucDB = Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES);
    if (isNucDB) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    } else {
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
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
        Sequence querySeq(par.maxSeqLen, querySequenceReader.getDbtype(), subMat, useProfileSearch ? 0 : par.kmerSize, par.spacedKmer,
                          false, !useProfileSearch);
        Sequence targetSeq(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false);

        Indexer idx(subMat->alphabetSize - 1, par.kmerSize);

        int xdrop = par.xdrop;

        BlockAligner blockAligner(par.maxSeqLen, par.rangeMin, par.rangeMax,
                                  isNucDB ? -par.gapOpen.values.nucleotide() : -par.gapOpen.values.aminoacid(),
                                  isNucDB ? -par.gapExtend.values.nucleotide() : -par.gapExtend.values.aminoacid());

        char buffer[1024];

        std::vector<Matcher::result_t> results;
        results.reserve(300);

        std::string result;
        result.reserve(1000);

        std::vector<QueryTableEntry> queries;
        queries.reserve(300);

        BlockIterator it;

        unsigned long correct_count = 0;

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < resultReader.getSize(); ++i) {
            size_t targetKey = resultReader.getDbKey(i);
            unsigned int targetId = targetKey;
            const char *targetSeqData = targetSequenceReader.getData(targetId, thread_idx);
            const unsigned int targetSeqLen = targetSequenceReader.getSeqLen(targetId);
            targetSeq.mapSequence(targetId, targetKey, targetSeqData, targetSeqLen);

            unsigned int queryKey = -1;

            // TODO: this might be wasted if no single hit hit the target
            blockAligner.initQuery(&targetSeq);

            std::vector<Kmer> targetKmers;
            targetKmers.reserve(targetSeqLen - par.kmerSize);
            while (targetSeq.hasNextKmer()) {
                const unsigned char *kmer = targetSeq.nextKmer();
                targetKmers.emplace_back(idx.int2index(kmer, 0, par.kmerSize),
                                         targetSeq.getCurrentPosition());
            }
            std::sort(targetKmers.begin(), targetKmers.end(), kmerComparator);

            // TODO: prefetch next sequence

            char *data = resultReader.getData(i, thread_idx);
            it.reset(data);
            while (it.getNext(queries)) {
                for (size_t j = 0; j < queries.size(); ++j) {
                    QueryTableEntry &query = queries[j];
                    bool kmerFound = false;
//                    while (targetSeq.hasNextKmer()) {
//                        const unsigned char *kmer = targetSeq.nextKmer();
                    //                        idx.printKmer(idx.int2index(kmer, 0, par.kmerSize), par.kmerSize, subMat->num2aa);
                    //                        Debug(Debug::INFO) << "\n";
//                        if (query.Query.kmer == idx.int2index(kmer, 0, par.kmerSize)) {
//                            query.Result.diag = query.Query.kmerPosInQuery - targetSeq.getCurrentPosition();
//                            kmerFound = true;
//                            break;
//                        }

                    const auto kmerCandidates = std::equal_range(targetKmers.begin(), targetKmers.end(),
                                                                 Kmer(query.Query.kmer, query.Query.kmerPosInQuery),
                                                                 [](const Kmer &kmer1, const Kmer &kmer2) {
                                                                     return kmer1.kmer < kmer2.kmer;
                                                                 });
//                    }
                    for (auto i = kmerCandidates.first; i != kmerCandidates.second; ++i) {
                        if (query.Query.kmer == i->kmer) {
                            query.Result.diag = query.Query.kmerPosInQuery - i->kmerPos;
                            kmerFound = true;
                            break;
                        }
                    }
//                    kmerFound = kmer != targetKmers.end() && query.Query.kmer == kmer->kmer;
                    if (kmerFound == false) {
//                        query.Result.diag = query.Query.kmerPosInQuery - kmer->kmerPos;
//                    } else {
                        Debug(Debug::ERROR) << "Found no matching k-mers between query and target sequence.\n";
                        EXIT(EXIT_FAILURE);
                    }
//                    targetSeq.resetCurrPos();
                }

                std::sort(queries.begin(), queries.end(), blockByDiagSort);
                if (isWithinNDiagonals(queries, 4) == false) {
                    continue;
                }

                queryKey = queries[0].querySequenceId;
                unsigned int queryId = querySequenceReader.getId(queryKey);
                const char *querySeqData = querySequenceReader.getData(queryId, thread_idx);
                const unsigned int querySeqLen = querySequenceReader.getSeqLen(queryId);
                std::string realSeq = SRAUtil::extractProfileSequence(querySeqData, querySeqLen, subMat);
                querySeq.mapSequence(queryId, queryKey, querySeqData, querySeqLen);

                if (useProfileSearch && realSeq.length() != querySeqLen) {
                    Debug(Debug::ERROR) << "Qeury seq len is wrong!\nCorrect count: " << correct_count << "\n";
                    Debug(Debug::ERROR) << "Retrieved sequence length: " << querySeqLen << "\n";
                    Debug(Debug::ERROR) << "Newly measured sequence length: " << realSeq.length() << "\n";
                    EXIT(EXIT_FAILURE);
                }

                correct_count++;
                DistanceCalculator::LocalAlignment aln = ungappedDiagFilter(queries,
                                                                            useProfileSearch ? realSeq.c_str()
                                                                                             : querySeqData,
                                                                            querySeqLen,
                                                                            targetSeqData,
                                                                            targetSeqLen,
                                                                            fastMatrix.matrix,
                                                                            evaluer,
                                                                            par.rescoreMode,
                                                                            par.evalThr);
                if (aln.diagonal == (int) INVALID_DIAG) {
                    continue;
                }

                Matcher::result_t res = blockAligner.align(&querySeq, aln, &evaluer, xdrop, subMat, useProfileSearch);
                res.dbKey = targetKey;
                res.queryOrfStartPos = queryKey;

                if (res.eval <= par.evalThr) {
                    results.emplace_back(res);
                }
//                Debug(Debug::INFO) << "Backtrace: " << res.backtrace << "\n";
//                Debug(Debug::INFO) << printAlnFromBt(targetSeqData, res.qStartPos, res.backtrace, false) << "\t"
//                                   << targetKey
//                                   << "\t" << res.qStartPos << "\t" << targetSeqLen << "\n";
//                Debug(Debug::INFO) << printAlnFromBt(querySeqData, res.dbStartPos, res.backtrace, true) << "\t"
//                                   << queryKey
//                                   << "\t" << res.dbStartPos << "\t" << querySeqLen << "\n" << res.eval << "\t"
//                                   << res.alnLength
//                                   << "\n\n";
            }

            progress.updateProgress();
//             std::sort(results.begin(), results.end(), Matcher::compareHits);
//            for (size_t j = 0; j < results.size(); ++j) {
//                size_t len = Matcher::resultToBuffer(buffer, results[j], false, false);
//                result.append(buffer, len);
//            }

            // writer.writeData(result.c_str(), result.size(), targetKey, thread_idx);
            //results.clear();
            //result.clear();
        }

        SORT_PARALLEL(results.begin(), results.end(), matcherResultsSort);

        for (size_t i = 0; i < results.size(); ++i) {
            results[i].dbOrfStartPos = (int) results[i].dbKey;
            results[i].dbKey = (unsigned int) results[i].queryOrfStartPos;
            Matcher::result_t::swapResult(results[i], evaluer, true);
            size_t len = Matcher::resultToBuffer(buffer, results[i], false, false, true);
            writer.writeData(buffer, len, results[i].dbKey, thread_idx);
        }
    }



//#pragma omp parallel
//    {
//        unsigned int thread_idx = 0;
//#ifdef OPENMP
//        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
//#endif
//        char buffer[1024];
//#pragma omp for schedule(dynamic, 1)
//        for (size_t i = 0; i < results.size(); ++i) {
////            results[i].dbOrfStartPos = (int) results[i].dbKey;
////            results[i].dbKey = (unsigned int) results[i].queryOrfStartPos;
//            Matcher::result_t::swapResult(results[i], evaluer, true);
//            size_t len = Matcher::resultToBuffer(buffer, results[i], false, false);
//            writer.writeData(buffer, len, results[i].dbKey, thread_idx);
//        }
//    }

    Debug(Debug::INFO) << "Compute Alignment finished, time spent: " << timer.lap() << "\n";

    writer.close(true);

    resultReader.close();
    querySequenceReader.close();
    targetSequenceReader.close();

    delete[] fastMatrix.matrixData;
    delete[] fastMatrix.matrix;
    delete subMat;
    subMat = nullptr;
    return EXIT_SUCCESS;
}

int blockByDiagSort(const QueryTableEntry &first, const QueryTableEntry &second) {
    if (first.Result.diag < second.Result.diag) {
        return true;
    }
    if (second.Result.diag < first.Result.diag) {
        return false;
    }

    if (first.Result.score < second.Result.score) {
        return true;
    }
    if (second.Result.score < first.Result.score) {
        return false;
    }

    if (first.Result.eval > second.Result.eval) {
        return true;
    }
    if (second.Result.eval > first.Result.eval) {
        return false;
    }

    if (first.querySequenceId < second.querySequenceId) {
        return true;
    }
    if (second.querySequenceId < first.querySequenceId) {
        return false;
    }

    if (first.targetSequenceID < second.targetSequenceID) {
        return true;
    }
    if (second.targetSequenceID < first.targetSequenceID) {
        return false;
    }

    return false;
}

bool matcherResultsSort(const Matcher::result_t &first, const Matcher::result_t &second) {
    unsigned int firstQueryKey = (unsigned int) first.queryOrfStartPos;
    unsigned int secondQueryKey = (unsigned int) second.queryOrfStartPos;
    if (firstQueryKey != secondQueryKey) {
        return firstQueryKey < secondQueryKey;
    }
    if (first.eval != second.eval) {
        return first.eval < second.eval;
    }
    if (first.score != second.score) {
        return first.score > second.score;
    }
    if (first.dbLen != second.dbLen) {
        return first.dbLen < second.dbLen;
    }
    return first.dbKey < second.dbKey;
}

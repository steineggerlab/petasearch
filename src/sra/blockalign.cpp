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
    DistanceCalculator::LocalAlignment alignmentResult;
    unsigned int lastDiagonal = INVALID_DIAG;
    for (size_t i = 1; i < queries.size(); ++i) {
        if (queries[i].Result.diag == lastDiagonal) {
            lastDiagonal = queries[i].Result.diag;
            continue;
        }
        lastDiagonal = queries[i].Result.diag;

        alignmentResult = DistanceCalculator::computeUngappedAlignment(
            querySeqData, querySeqLen, targetSeqData, targetSeqLen,
            queries[i].Result.diag, matrix, rescoreMode
        );

        if (alignmentResult.startPos < 0 || alignmentResult.endPos < 0) {
            continue;
        }

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
    unsigned long long kmer = -1;
    int kmerPos = -1;

    Kmer() = default;
    Kmer(unsigned long long kmer, int kmerPos) : kmer(kmer), kmerPos(kmerPos) {}
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

int blockalign(int argc, const char **argv, const Command &command) {
    Timer timer;
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.spacedKmer = false;
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    par.evalThr = 1000000;
    par.parseParameters(argc, argv, command, true, 0, 0);

    const int inputSeqType = FileUtil::parseDbType(par.db1.c_str());
    const bool useProfileSearch = Parameters::isEqualDbtype(inputSeqType, Parameters::DBTYPE_HMM_PROFILE);

    // query target prev_result new_result is the order
    DBReader<unsigned int> querySequenceReader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                               DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    querySequenceReader.open(DBReader<unsigned int>::NOSORT);

    SRADBReader targetSequenceReader(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                     DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    targetSequenceReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                        DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX |
                                        DBReader<unsigned int>::USE_WRITABLE);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed,
                    Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();

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

    const int xdrop = par.xdrop;

    size_t kmerMatch = 0;
    size_t ungappedNum = 0;
    size_t alignmentsNum = 0;
    size_t totalPassedNum = 0;
    Debug::Progress progress(resultReader.getSize());
#pragma omp parallel reduction(+:kmerMatch, ungappedNum, alignmentsNum, totalPassedNum)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Sequence targetSeq(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false, false, par.spacedKmerPattern);

        Indexer idx(subMat->alphabetSize - 1, par.kmerSize);

        BlockAligner blockAligner(
            par.maxSeqLen, par.rangeMin, par.rangeMax,
            isNucDB ? -par.gapOpen.values.nucleotide() : -par.gapOpen.values.aminoacid(),
            isNucDB ? -par.gapExtend.values.nucleotide() : -par.gapExtend.values.aminoacid(),
            querySequenceReader.getDbtype()
        );

        // Sequence querySeq(par.maxSeqLen, querySequenceReader.getDbtype(), subMat, 0, false, false, false);
        // Matcher matcher(querySeq.getSeqType(), targetSeq.getSeqType(), par.maxSeqLen, subMat, &evaluer, false, 1.0, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(), 1.0, 0);

        char buffer[1024];

        std::vector<Matcher::result_t> results;
        results.reserve(300);

        std::string result;
        result.reserve(1000);

        std::string realSeq;
        realSeq.reserve(1000);

        std::vector<QueryTableEntry> queries;
        queries.reserve(300);

        BlockIterator it;

        unsigned long correct_count = 0;
#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < resultReader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int targetKey = resultReader.getDbKey(i);
            const unsigned int targetSeqLen = targetSequenceReader.getSeqLen(targetKey);

            if (targetSeqLen < (unsigned int)par.kmerSize) {
                continue;
            }

            const char *targetSeqData = targetSequenceReader.getData(targetKey, thread_idx);
            targetSeq.mapSequence(targetKey, targetKey, targetSeqData, targetSeqLen);

            unsigned int queryKey = -1;

            bool isBlockAlignerInit = false;

            std::vector<Kmer> targetKmers;
            targetKmers.reserve(targetSeqLen - par.kmerSize);
            while (targetSeq.hasNextKmer()) {
                const unsigned char *kmer = targetSeq.nextKmer();
                targetKmers.emplace_back(idx.int2index(kmer, 0, par.kmerSize), targetSeq.getCurrentPosition());
            }
            SORT_SERIAL(targetKmers.begin(), targetKmers.end(), kmerComparator);

            // TODO: prefetch next sequence
            char *data = resultReader.getData(i, thread_idx);
            it.reset(data);
            while (it.getNext(queries)) {
                for (size_t j = 0; j < queries.size(); ++j) {
                    QueryTableEntry &query = queries[j];
                    const auto kmer = std::lower_bound(targetKmers.begin(), targetKmers.end(),
                                                       Kmer(query.Query.kmer, query.Query.kmerPosInQuery),
                                                       [](const Kmer &kmer1, const Kmer &kmer2) {
                                                           return kmer1.kmer < kmer2.kmer;
                                                       });
                    bool kmerFound = kmer != targetKmers.end() && query.Query.kmer == kmer->kmer;
                    if (kmerFound) {
                        query.Result.diag = query.Query.kmerPosInQuery - kmer->kmerPos;
                    } else {
                        Debug(Debug::ERROR) << "Found no matching k-mers between query and target sequence\n";
                        EXIT(EXIT_FAILURE);
                    }
                }
                kmerMatch++;

                SORT_SERIAL(queries.begin(), queries.end(), blockByDiagSort);
                if (isWithinNDiagonals(queries, 4) == false) {
                    continue;
                }

                queryKey = queries[0].querySequenceId;
                unsigned int queryId = querySequenceReader.getId(queryKey);
                const char *querySeqData = querySequenceReader.getData(queryId, thread_idx);
                const unsigned int querySeqLen = querySequenceReader.getSeqLen(queryId);
                const unsigned int queryEntryLen = querySequenceReader.getEntryLen(queryId);

                if (useProfileSearch) {
                    realSeq.clear();
                    Sequence::extractProfileSequence(querySeqData, queryEntryLen - 1, *subMat, realSeq);
                    if (realSeq.length() != querySeqLen) {
                        Debug(Debug::ERROR) << "Query sequence length is wrong!\n" 
                                            << "Correct count: " << correct_count << "\n"
                                            << "Retrieved sequence length: " << querySeqLen << "\n"
                                            << "Newly measured sequence length: " << realSeq.length() << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                }

                correct_count++;
                DistanceCalculator::LocalAlignment aln = ungappedDiagFilter(
                    queries,
                    useProfileSearch ? realSeq.c_str() : querySeqData,
                    querySeqLen,
                    targetSeqData,
                    targetSeqLen,
                    fastMatrix.matrix,
                    evaluer,
                    par.rescoreMode,
                    par.evalThr
                );
                ungappedNum++;

                // for (size_t i = 0; i < queries.size(); ++i) {
                //     Debug(Debug::WARNING) << queries[i].querySequenceId << "\t" << queries[i].targetSequenceID << "\n";
                //     Debug(Debug::WARNING) << queries[i].Query.kmerPosInQuery << "\t" << queries[i].Query.kmer << "\n";
                // }
                // Debug(Debug::WARNING) << "querySeqData: " << querySeqData << "\n";
                // Debug(Debug::WARNING) << "targetSeqData: " << targetSeqData << "\n";
                // Debug(Debug::WARNING) << "diagonal: " << aln.diagonal << "\n";
                // Debug(Debug::WARNING) << "score: " << aln.score << "\n";
                // Debug(Debug::WARNING) << "distToDiagonal: " << aln.distToDiagonal << "\n";
                // Debug(Debug::WARNING) << "startPos: " << aln.startPos << "\n";
                // Debug(Debug::WARNING) << "endPos: " << aln.endPos << "\n";
                // unsigned int qUngappedStartPos = aln.startPos + ((aln.diagonal >= 0) ? aln.distToDiagonal : 0);
                // unsigned int qUngappedEndPos = aln.endPos + ((aln.diagonal >= 0) ? aln.distToDiagonal : 0);
                // unsigned int dbUngappedStartPos = aln.startPos + ((aln.diagonal >= 0) ? 0 : aln.distToDiagonal);
                // unsigned int dbUngappedEndPos = aln.endPos + ((aln.diagonal >= 0) ? 0 : aln.distToDiagonal);
                // Debug(Debug::INFO) << "qUnGapStart: " << qUngappedStartPos << "\n";
                // Debug(Debug::INFO) << "qUnGapEnd: " << qUngappedEndPos << "\n";
                // Debug(Debug::INFO) << "dbUnGapStart: " << dbUngappedStartPos << "\n";
                // Debug(Debug::INFO) << "dbUnGapEnd: " << dbUngappedEndPos << "\n";
                // Debug(Debug::INFO) << std::string(querySeqData + qUngappedStartPos, qUngappedEndPos - qUngappedStartPos) << "\n";
                // Debug(Debug::INFO) << std::string(targetSeqData + dbUngappedStartPos, dbUngappedEndPos - dbUngappedStartPos) << "\n";
                // EXIT(EXIT_FAILURE);

                if (aln.diagonal == (int) INVALID_DIAG || aln.startPos < 0 || aln.endPos < 0) {
                    continue;
                }

                if (isBlockAlignerInit == false) {
                    isBlockAlignerInit = true;
                }

                // querySeq.mapSequence(queryId, queryKey, querySeqData, querySeqLen);
                // matcher.initQuery(&querySeq);
                // Matcher::result_t res = matcher.getSWResult(&targetSeq, INT_MAX, false, 0, 0.0, par.evalThr, Matcher::SCORE_COV_SEQID, 0, false);
                Matcher::result_t res = blockAligner.align(targetSeq.getSeqData(), targetSeq.L, querySeqData, querySeqLen, aln, &evaluer, xdrop, *subMat);
                res.dbKey = targetKey;
                res.queryOrfStartPos = queryKey;
                alignmentsNum++;

                if (res.eval <= par.evalThr) {
                    results.emplace_back(res);
                    totalPassedNum++;
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
//            std::sort(results.begin(), results.end(), Matcher::compareHits);
//            for (size_t j = 0; j < results.size(); ++j) {
//                size_t len = Matcher::resultToBuffer(buffer, results[j], false, false);
//                result.append(buffer, len);
//            }
//            writer.writeData(result.c_str(), result.size(), targetKey, thread_idx);
//            results.clear();
//            result.clear();
        }

        SORT_SERIAL(results.begin(), results.end(), matcherResultsSort);

        for (size_t i = 0; i < results.size(); ++i) {
            results[i].dbOrfStartPos = (int) results[i].dbKey;
            results[i].dbKey = (unsigned int) results[i].queryOrfStartPos;
            Matcher::result_t::swapResult(results[i], evaluer, true);
            size_t len = Matcher::resultToBuffer(buffer, results[i], par.addBacktrace, true, true);
            writer.writeData(buffer, len, results[i].dbKey, thread_idx);
        }
        results.clear();
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
    writer.close(true);

    Debug(Debug::INFO) << kmerMatch << " before diagonal filter\n";
    Debug(Debug::INFO) << ungappedNum << " ungapped alignments calculated\n";
    Debug(Debug::INFO) << alignmentsNum << " alignments calculated\n";
    Debug(Debug::INFO) << totalPassedNum << " sequence pairs passed the thresholds";
    if (alignmentsNum > 0) {
        Debug(Debug::INFO) << " (" << ((float) totalPassedNum / (float) alignmentsNum) << " of overall calculated)";
    }
    Debug(Debug::INFO) << "\n";
    size_t dbSize = querySequenceReader.getSize();
    if (dbSize > 0) {
        size_t hits = totalPassedNum / dbSize;
        size_t hits_rest = totalPassedNum % dbSize;
        float hits_f = ((float) hits) + ((float) hits_rest) / (float) dbSize;
        Debug(Debug::INFO) << hits_f << " hits per query sequence\n";
    }

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

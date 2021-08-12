#include <QueryTableEntry.h>
#include "LocalParameters.h"
#include "Command.h"
#include "Debug.h"
#include "IndexReader.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "NucleotideMatrix.h"
#include "EvalueComputation.h"
#include "DistanceCalculator.h"
#include "Matcher.h"
#include "QueryTableEntry.h"
#include "BlockAligner.h"

#include "omptl/omptl_algorithm"

#ifdef OPENMP

#include <omp.h>
#include <QueryTableEntry.h>

#endif


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

unsigned int ungappedDiagFilter(std::vector<QueryTableEntry> &queries,
                                const char *querySeqData,
                                size_t querySeqLen,
                                const char *targetSeqData,
                                size_t targetSeqLen,
                                const char **matrix,
                                EvalueComputation &evaluer,
                                int rescoreMode,
                                double evalThr) {

    int maxScore = INT_MIN;
    unsigned int lastDiagonal = INVALID_DIAG;
    for (size_t i = 1; i < queries.size(); ++i) {
        if (queries[i].Result.diag == lastDiagonal) {
            lastDiagonal = queries[i].Result.diag;
            continue;
        }
        lastDiagonal = queries[i].Result.diag;

        DistanceCalculator::LocalAlignment aln =
                DistanceCalculator::computeUngappedAlignment(
                        querySeqData, querySeqLen, targetSeqData, targetSeqLen,
                        queries[i].Result.diag, matrix, rescoreMode);
        queries[i].Result.score = aln.score;

//        Debug(Debug::ERROR) << querySeqLen << "\t" << std::string(querySeqData, querySeqLen) << "\n";
//        Debug(Debug::ERROR) << targetSeqLen << "\t" << std::string(k->Result.diag, ' ') << std::string(targetSeqData, targetSeqLen) << "\n";

        // if (protein) {
        // different than wiki, explain swap afterwards
        double eval = evaluer.computeEvalue(aln.score, querySeqLen);
//                    Debug(Debug::ERROR) << k->Result.diag  << "\t" << aln.score <<  "\t" << eval << "\n";

        // this is bad for nucleotide petasearch, we need to know the best diagonal
        if (eval <= evalThr) {
            maxScore = aln.score;
            break;
        }
        // } else {
        // actually calculate maxscore
        // }
    }

    if (maxScore == INT_MIN) {
        return INVALID_DIAG;
    }
    return lastDiagonal;
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
    DBReader<unsigned int> targetSequenceReader(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    targetSequenceReader.open(DBReader<unsigned int>::NOSORT);
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
        subMat = new NucleotideMatrix(par.scoringMatrixFile.nucleotides, 1.0, 0.0);
    } else {
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
        Sequence querySeq(par.maxSeqLen, querySequenceReader.getDbtype(), subMat, par.kmerSize, par.spacedKmer,
                          false, useProfileSearch ? false : true);
        Sequence targetSeq(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false);

        Indexer idx(subMat->alphabetSize - 1, par.kmerSize);

        // TODO: replace this matcher with block-aligner
//        Matcher matcher(seqType, par.maxSeqLen, subMat, &evaluer, (bool) par.compBiasCorrection,
//                        isNucDB ? par.gapOpen.nucleotides : par.gapOpen.aminoacids,
//                        isNucDB ? par.gapExtend.nucleotides : par.gapExtend.aminoacids);
        BlockAligner blockAligner(subMat, par.maxSeqLen,
                                  isNucDB ? par.gapOpen.nucleotides : par.gapOpen.aminoacids,
                                  isNucDB ? par.gapExtend.nucleotides : par.gapExtend.aminoacids);

        char buffer[1024];
        std::vector<Matcher::result_t> results;
        results.reserve(300);

        std::string result;
        result.reserve(1000);

        std::vector<QueryTableEntry> queries;
        queries.reserve(300);

        BlockIterator it;

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < resultReader.getSize(); ++i) {
            progress.updateProgress();

            size_t targetKey = resultReader.getDbKey(i);
            unsigned int targetId = targetSequenceReader.getId(targetKey);
            const char *targetSeqData = targetSequenceReader.getData(targetId, thread_idx);
            const unsigned int targetSeqLen = targetSequenceReader.getSeqLen(targetId);
            targetSeq.mapSequence(targetId, targetKey, targetSeqData, targetSeqLen);

            // TODO: this might be wasted if no single hit hit the target
            blockAligner.initQuery(&targetSeq);

            // TODO: prefetch next sequence

            char *data = resultReader.getData(i, thread_idx);
            it.reset(data);
            while (it.getNext(queries)) {
                for (size_t j = 0; j < queries.size(); ++j) {
                    QueryTableEntry &query = queries[j];
//                    Debug(Debug::INFO) << SSTR(queries.size()) << "\n";
//                    Debug(Debug::INFO) << SSTR(targetKey) << "\n";
//                    Debug(Debug::INFO) << SSTR(query.querySequenceId) << "\n";
//                    Debug(Debug::INFO) << SSTR(query.Query.kmer) << "\n";
//                    Debug(Debug::INFO) << SSTR(i) << "\n";
//                    Debug(Debug::INFO) << SSTR(j) << "\n";
//                    idx.printKmer(query.Query.kmer, par.kmerSize, subMat->int2aa);
//                    Debug(Debug::INFO) << "\n";
                    bool kmerFound = false;
                    while (targetSeq.hasNextKmer()) {
                        const unsigned char *kmer = targetSeq.nextKmer();
//                        idx.printKmer(idx.int2index(kmer, 0, par.kmerSize), par.kmerSize, subMat->int2aa);
//                        Debug(Debug::INFO) << "\n";
                        if (query.Query.kmer == idx.int2index(kmer, 0, par.kmerSize)) {
                            query.Result.diag = query.Query.kmerPosInQuery - targetSeq.getCurrentPosition();
                            kmerFound = true;
                            break;
                        }
                    }
                    if (kmerFound == false) {
                        Debug(Debug::ERROR) << "Found no matching k-mers between query and target sequence.\n";
                        EXIT(EXIT_FAILURE);
                    }
                    targetSeq.resetCurrPos();
                }

                std::sort(queries.begin(), queries.end(), blockByDiagSort);
                if (isWithinNDiagonals(queries, 4) == false) {
                    continue;
                }

                unsigned int queryKey = queries[0].querySequenceId;
                unsigned int queryId = querySequenceReader.getId(queryKey);
                const char *querySeqData = querySequenceReader.getData(queryId, thread_idx);
                const unsigned int querySeqLen = querySequenceReader.getSeqLen(queryId);
                querySeq.mapSequence(queryId, queryKey, querySeqData, querySeqLen);
                std::string realSeq;
                querySeq.extractProfileSequence(querySeqData, *subMat, realSeq);

                unsigned int diag = ungappedDiagFilter(queries,
                                                       realSeq.c_str(),
                                                       querySeqLen,
                                                       targetSeqData,
                                                       targetSeqLen,
                                                       fastMatrix.matrix,
                                                       evaluer,
                                                       par.rescoreMode,
                                                       par.evalThr);
                if (diag == INVALID_DIAG) {
                    continue;
                }

                // TODO we have to swap coverage mode either here or already in workflow etc
                querySeq.mapSequence(queryId, queryKey, querySeqData, querySeqLen);
                Matcher::result_t res = blockAligner.align(&querySeq, diag, &evaluer);
                results.emplace_back(res);

                if (thread_idx == 0 && i < 640) {
                    Debug(Debug::INFO) << res.backtrace << "\n";
                    Debug(Debug::INFO) << printAlnFromBt(targetSeqData, res.qStartPos, res.backtrace, false) << "\t" << targetKey
                    << "\t" << res.qStartPos << "\t" << targetSeqLen << "\n";
                    Debug(Debug::INFO) << printAlnFromBt(querySeqData, res.dbStartPos, res.backtrace, true) << "\t" << queryKey
                    << "\t" << res.dbStartPos << "\t" << querySeqLen << "\n" << res.eval << "\t" << res.alnLength
                    << "\n\n";
                }
            }

            std::sort(results.begin(), results.end(), Matcher::compareHits);
            for (size_t j = 0; j < results.size(); ++j) {
                size_t len = Matcher::resultToBuffer(buffer, results[j], false, false);
                result.append(buffer, len);
            }

            writer.writeData(result.c_str(), result.size(), targetKey, thread_idx);
            results.clear();
            result.clear();
        }
    }

    Debug(Debug::INFO) << "Compute Alignment finished, time spent: " << timer.lap() << "\n";

    writer.close();

    resultReader.close();
    targetSequenceReader.close();

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

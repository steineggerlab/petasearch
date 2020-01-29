#ifndef QUERY_TABLE_ENTRY_H
#define QUERY_TABLE_ENTRY_H

#include "Debug.h"
#include "Util.h"
#include "itoa.h"

struct __attribute__((__packed__)) QueryTableEntry
{
    unsigned int querySequenceId;
    unsigned int targetSequenceID;
    union {
        struct __attribute__((__packed__)) {
            unsigned int kmerPosInQuery;
            unsigned long long kmer;
        } Query;
        struct __attribute__((__packed__)) {
            unsigned int diag;
            unsigned int score;
            unsigned int eval;
        } Result;
    };

    static size_t queryEntryToBuffer(char *buff1, QueryTableEntry &h) {
        char * basePos = buff1;
        char * tmpBuff = Itoa::u32toa_sse2((uint32_t) h.querySequenceId, buff1);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::u32toa_sse2((uint32_t) h.Query.kmerPosInQuery, tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::u64toa_sse2((uint64_t) h.Query.kmer, tmpBuff);
        *(tmpBuff-1) = '\n';
        *(tmpBuff) = '\0';
        return tmpBuff - basePos;
    }

    static QueryTableEntry parseQueryEntry(const char* data) {
        QueryTableEntry result;
        const char *wordCnt[3];
        size_t cols = Util::getWordsOfLine(data, wordCnt, 3);
        if (cols == 3) {
            result.querySequenceId = Util::fast_atoi<unsigned int>(wordCnt[0]);
            result.Query.kmerPosInQuery = Util::fast_atoi<unsigned int>(wordCnt[1]);
            result.Query.kmer = Util::fast_atoi<unsigned long long>(wordCnt[2]);
        } else {
            Debug(Debug::INFO) << "Invalid query entry cols = " << cols << " wordCnt[0]: " << wordCnt[0] << "\n" ;
            EXIT(EXIT_FAILURE);
        }
        return result;
    }
};

#endif
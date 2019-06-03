#ifndef QUERY_TABLE_ENTRY_H
#define QUERY_TABLE_ENTRY_H
struct __attribute__((__packed__)) QueryTableEntry
{
    unsigned int querySequenceId;
    unsigned int targetSequenceID;
    union {
        struct __attribute__((__packed__)) {
            unsigned int kmerPosInQuerry;
            unsigned long kmer;
        } Query;
        struct __attribute__((__packed__)) {
            unsigned int diag;
            unsigned int score;
            unsigned int eval;
        } Result;
    };
};

#endif
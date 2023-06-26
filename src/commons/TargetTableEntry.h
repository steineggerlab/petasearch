#ifndef TARGET_TABLE_ENTRY_H
#define TARGET_TABLE_ENTRY_H
struct __attribute__((__packed__)) TargetTableEntry
{   
    unsigned long long kmerAsLong;
    unsigned int sequenceID;
    unsigned int sequenceLength;
};
#endif
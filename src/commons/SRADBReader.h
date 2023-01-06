//
// Created by matchy233 on 1/8/21.
//
#ifndef SRASEARCH_SRADBREADER_H
#define SRASEARCH_SRADBREADER_H

#include "DBReader.h"
#include "MemoryTracker.h"
#include "Sequence.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"

#include <cstddef>
#include <utility>
#include <vector>
#include <string>

class SRADBReader : public MemoryTracker {
public:
    SRADBReader(const char* dataFileName, const char* indexFileName, int threads, int mode);
    ~SRADBReader();
    void open(int accessType);
    void close();
    int getDbtype();
    size_t getSize();
    size_t getSeqLen(size_t id);
    unsigned int getDbKey(size_t id);
    char *getData(size_t id, int thread_idx);
    size_t getAminoAcidDBSize();
    void readIndex(char *data, size_t indexDataSize, unsigned long *index);
private:
    int closed;
    int threads;
    int dataMode;
    bool dataMapped;
    int accessType;
    int dbType;
    bool isHeader;
    void checkClosed();

    size_t size;

    char *dataFileName;
    std::vector<std::string> dataFileNames;
    size_t totalDataSize;
    size_t dataFileCnt;

    char *indexFileName;
    unsigned long *index;

    char **dataFiles;
    size_t *dataSizeOffset;

    char **seqBuffer;
    size_t maxSeqLen;

    char *mmapData(FILE *file, size_t *dataSize);
    void unmapData();
    void setSequentialAdvice();

    char *getDataUncompressed(size_t id);

    char *getDataByOffset(size_t offset);
};

#endif

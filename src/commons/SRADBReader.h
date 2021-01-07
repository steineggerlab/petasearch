//
// Created by matchy233 on 1/8/21.
//

#ifndef SRASEARCH_SRADBREADER_H
#define SRASEARCH_SRADBREADER_H

#include "DBReader.h"


class SRADBReader {
public:
    SRADBReader(const char* dataFileName, const char* indexFileName, int threads, int mode);
    ~SRADBReader();
    bool open(int accessType);
    void close();
    int getDbtype();
    size_t getSize();
    size_t getSeqLen(size_t id);
    unsigned int getDbKey(size_t id);
    char *getData(size_t id, int thread_idx);
    size_t getAminoAcideDBSize();
private:
    int closed;
    DBReader<unsigned int> *__dbreader = NULL;
    void checkClosed();
    char* dataFileName;
    char* indexFileName;
    char** dataFiles;
    size_t * dataSizeOffset;
};


#endif //SRASEARCH_SRADBREADER_H

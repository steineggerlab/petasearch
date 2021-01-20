//
// Created by matchy233 on 1/8/21.
//

#include "SRADBReader.h"
#include "BitManipulateMacros.h"

SRADBReader::SRADBReader(const char *dataFileName, const char *indexFileName, int threads, int mode) {
    _dbreader = new DBReader<unsigned int>(dataFileName, indexFileName, threads, mode);
}

bool SRADBReader::open(int accessType) {
    return _dbreader->open(accessType);
}

void SRADBReader::close() {
    _dbreader->close();
}

int SRADBReader::getDbtype() {
    return _dbreader->getDbtype();
}

char *SRADBReader::getData(size_t id, int thread_idx) {
    char *rawString = _dbreader->getData(id, thread_idx);
    unsigned short *packedArray = reinterpret_cast<unsigned short *>(rawString);
    int packedLen = sizeof(rawString) / 2;
    char *resString = (char *)malloc(packedLen * 3);
    int i = 0;
    int j = 0;
    for (; i < packedLen - 1; i++) {
        resString[j] = GET_HIGH_CHAR(packedArray[i]);
        resString[j+1] = GET_MID_CHAR(packedArray[i]);
        resString[j+2] = GET_LOW_CHAR(packedArray[i]);
        j += 3;
    }
    char p1 = GET_HIGH_CHAR(packedArray[i]);
    char p2 = GET_MID_CHAR(packedArray[i]);
    char p3 = GET_LOW_CHAR(packedArray[i]);
    int count = (p1 == 0) + (p2 == 0);
    resString[j+count] = p3;
    if (p2 != 0) {
        resString[j+count-1] = p2;
    }
    if (p3 != 0) {
        resString[j+count-2] = p1;
    }
    return resString;
}

size_t SRADBReader::getSize() {
    return _dbreader->getSize();
}

size_t SRADBReader::getSeqLen(size_t id) {
    return id;
}

unsigned int SRADBReader::getDbKey(size_t id) {
    return _dbreader->getDbKey(id);
}

size_t SRADBReader::getAminoAcidDBSize() {
    return _dbreader->getAminoAcidDBSize() * 3 / 2;
}

void SRADBReader::checkClosed() {
    if (closed == 1) {
        Debug(Debug::ERROR) << "Trying to read a closed database.\n";
        EXIT(EXIT_FAILURE);
    }
}

SRADBReader::~SRADBReader() {
    if (dataFileName != NULL) {
        free(dataFileName);
    }

    if (indexFileName != NULL) {
        free(indexFileName);
    }

    if (dataSizeOffset != NULL) {
        delete[] dataSizeOffset;
    }

    if (dataFiles != NULL) {
        delete[] dataFiles;
    }
    delete _dbreader;
}

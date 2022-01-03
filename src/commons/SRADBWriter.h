
#ifndef SRASEARCH_SRADBWRITER_H
#define SRASEARCH_SRADBWRITER_H
// For parallel write access, one each thread creates its own DB
// After the parallel calculation are done, all DBs are merged into single DB

#include <string>
#include <vector>

#include "DBReader.h"
#include "MemoryTracker.h"

template <typename T> class DBReader;

class SRADBWriter : public MemoryTracker  {
public:
    SRADBWriter(const char* dataFileName, const char* indexFileName, unsigned int threads, size_t mode, int dbtype);

    ~SRADBWriter();

    void open(size_t bufferSize = SIZE_MAX);

    void close(bool merge = false);

    char* getDataFileName() { return dataFileName; }

    char* getIndexFileName() { return indexFileName; }

    void writeStart(unsigned int thrIdx = 0);
    size_t writeAdd(const char* data, size_t dataSize, unsigned int thrIdx = 0);
    void writeEnd(unsigned int key, unsigned int thrIdx = 0, bool addNullByte = true, bool addIndexEntry = true);

    void writeData(const char *data, size_t dataSize, unsigned int key, unsigned int threadIdx = 0, bool addNullByte = true, bool addIndexEntry = true);

    static size_t indexToBuffer(char *buff1, size_t offsetStart);

    static void mergeResults(const std::string &outFileName, const std::string &outFileNameIndex,
                             const std::vector<std::pair<std::string, std::string>> &files);

    void writeIndexEntry(unsigned int key, size_t offset, size_t length, unsigned int thrIdx);

    static void writeDbtypeFile(const char* path, int dbtype, bool isCompressed);

    size_t getStart(unsigned int threadIdx){
        return starts[threadIdx];
    }

    size_t getOffset(unsigned int threadIdx){
        return offsets[threadIdx];
    }
//
//    template <typename T>
//    static void writeIndex(FILE *outFile, size_t indexSize, T *index);
//
//    template <typename T>
//    static void writeIndexEntryToFile(FILE *outFile, char *buff1, T &index);

    bool isClosed(){
        return closed;
    }
private:
    void checkClosed();

    char *makeResultFilename(const char* name, size_t split);

    static void mergeResults(const char *outFileName, const char *outFileNameIndex,
                             const char **dataFileNames, const char **indexFileNames,
                             unsigned long fileCount, bool mergeDatafiles);

    static void mergeIndex(const char** indexFilenames, unsigned int fileCount, const std::vector<size_t> &dataSizes);

    char* dataFileName;
    char* indexFileName;

    FILE** dataFiles;
    char** dataFilesBuffer;
    size_t bufferSize;
    FILE** indexFiles;

    char** dataFileNames;
    char** indexFileNames;

    size_t* starts;
    size_t* offsets;

    const unsigned int threads;
    const size_t mode;
    int dbtype;

    bool closed;

    std::string datafileMode;
};

#endif //SRASEARCH_SRADBWRITER_H

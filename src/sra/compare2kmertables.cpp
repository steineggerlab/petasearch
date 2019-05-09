#include "Parameters.h"
#include "Command.h"
#include "Debug.h"
#include <sys/stat.h>
#include <sys/mman.h>
#include "MathUtil.h"


size_t kmer2index(long kmer);

int compare2kmertables(int argc, const char **argv, const Command& command){
    Parameters& par = Parameters::getInstance();
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, 1, false);

    FILE* handleQuerryKmerTable = fopen(argv[0],"rb");
    int fdQuerryTable = fileno(handleQuerryKmerTable);
    struct stat fileStatsQuerryTable;
    fstat(fdQuerryTable, &fileStatsQuerryTable);
    size_t fileSizeQuerryTable = fileStatsQuerryTable.st_size;

    FILE* handleTargetKmerTable = fopen(argv[1],"rb");
    int fdTargetTable = fileno(handleTargetKmerTable);
    struct stat fileStatsTargetTable;
    fstat(fdTargetTable,&fileStatsTargetTable);
    size_t fileSizeTargetTable = fileStatsTargetTable.st_size;

    long* startPosQuerryTable = (long*) mmap(NULL, fileSizeQuerryTable, PROT_READ,MAP_PRIVATE, fdQuerryTable, 0);
    long* startPosTargetTable = (long*) mmap(NULL, fileSizeTargetTable, PROT_READ,MAP_PRIVATE, fdTargetTable, 0);
    if (posix_madvise (startPosQuerryTable, fileSizeQuerryTable, POSIX_MADV_SEQUENTIAL|POSIX_MADV_WILLNEED) != 0){
        Debug(Debug::ERROR) << "posix_madvise returned an error for k-mer table 2 " << "\n";
    } 
    if (posix_madvise (startPosTargetTable, fileSizeTargetTable, POSIX_MADV_SEQUENTIAL|POSIX_MADV_WILLNEED) != 0){
        Debug(Debug::ERROR)  << "posix_madvise returned an error for k-mer table 2 " << "\n";
    }

    long* currentQuerryPos = startPosQuerryTable;
    long* currentTargetPos = startPosTargetTable;
    long* endTargetPos = startPosTargetTable + fileSizeTargetTable/sizeof(long);
    size_t equalKmers = 0;

    while(currentTargetPos <= endTargetPos){
        if(*currentQuerryPos == *currentTargetPos){
            //Match found
            ++equalKmers;
            ++currentTargetPos;
        }
        while (*currentQuerryPos < *currentTargetPos)
        {
            ++currentQuerryPos;
        }
        while (*currentTargetPos < *currentQuerryPos)
        {
            ++currentTargetPos;
        }
    }
    munmap(startPosQuerryTable, fileSizeQuerryTable);
    munmap(startPosTargetTable, fileSizeTargetTable);
    fclose(handleQuerryKmerTable);
    fclose(handleTargetKmerTable);

    Debug(Debug::INFO)<<"number of equal Kmers: "<<equalKmers<<"\n";
    return equalKmers;
}


size_t kmer2index(long kmer){
    //hornerschema verwenden, bei der division drauf achten,dass wir nicht größer als die alphabet größer werden.
}

// long readNextQuerryKmer(){
//     if(currentQuerryPos<file)
// }
// long readNextTargetKmer(){

// }
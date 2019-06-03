#include "Parameters.h"
#include "Command.h"
#include "Debug.h"
#include <sys/stat.h>
#include <sys/mman.h>
#include "MathUtil.h"

int compare2kmertables(int argc, const char **argv, const Command& command){
    Parameters& par = Parameters::getInstance();
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, 2, false);

    FILE* handleQuerryKmerTable = fopen(par.db1.c_str(),"rb");
    int fdQuerryTable = fileno(handleQuerryKmerTable);
    struct stat fileStatsQuerryTable;
    fstat(fdQuerryTable, &fileStatsQuerryTable);
    size_t fileSizeQuerryTable = fileStatsQuerryTable.st_size;

    FILE* handleTargetKmerTable = fopen(par.db2.c_str(),"rb");
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

    struct timeval startTime;
    struct timeval endTime; 
    gettimeofday(&startTime, NULL);

    //maybe faster:
    // Type your code here, or load an example.
// int square(long * currentTargetPos, long * endTargetPos, long * currentQuerryPos) {
//      long * targetBasePos = currentTargetPos;
//       while(__builtin_expect(currentTargetPos <= endTargetPos, 1)){
//         currentTargetPos += (*currentQuerryPos == *currentTargetPos);
//         //kmer->tidx =  currentTargetPos - targetBasePos;
//         //kmer += (*currentQuerryPos == *currentTargetPos);
//         currentQuerryPos += (*currentQuerryPos < *currentTargetPos);
//         currentTargetPos += (*currentTargetPos < *currentQuerryPos);
//     }
// }
    while(currentTargetPos <= endTargetPos){
        if(*currentQuerryPos == *currentTargetPos){
            //Match found
            ++equalKmers;
            ++currentTargetPos;
        }
        while (*currentQuerryPos < *currentTargetPos){
            ++currentQuerryPos;
        }
        while (*currentTargetPos < *currentQuerryPos){
            ++currentTargetPos;
        }
    }
    gettimeofday(&endTime, NULL);
    double timediff = (endTime.tv_sec - startTime.tv_sec) + 1e-6 * (endTime.tv_usec - startTime.tv_usec);
    munmap(startPosQuerryTable, fileSizeQuerryTable);
    munmap(startPosTargetTable, fileSizeTargetTable);
    fclose(handleQuerryKmerTable);
    fclose(handleTargetKmerTable);
    Debug(Debug::INFO) << timediff<<" s; Rate "<<((fileSizeTargetTable+fileSizeQuerryTable)/1e+9)/timediff
            <<" GB/s \n";
    Debug(Debug::INFO)<<"number of equal Kmers: "<<equalKmers<<"\n";
    return equalKmers;
}
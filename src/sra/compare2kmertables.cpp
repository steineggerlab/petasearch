#include "Parameters.h"
#include "Command.h"
#include "Debug.h"
#include <sys/stat.h>
#include <sys/mman.h>
#include "MathUtil.h"


#define KMER_TABLE_1 "tmp/kmerTableDecodedAsLong"
#define KMER_TABLE_2 "tmp/kmerTableDecodedAsLong2"

size_t kmer2index(long kmer);

int compare2kmertables(int argc, const char **argv, const Command& command){
    Parameters& par = Parameters::getInstance();
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, 1, false);
    


    FILE* handleKmerTable1 = fopen(KMER_TABLE_1,"rb");
    int fdTable1 = fileno(handleKmerTable1);
    struct stat fileStatsTable1;
    fstat(fdTable1, &fileStatsTable1);
    size_t fileSizeTable1 = fileStatsTable1.st_size;

    FILE* handleKmerTable2 = fopen(KMER_TABLE_2,"rb");
    int fdTable2 = fileno(handleKmerTable2);
    struct stat fileStatsTable2;
    fstat(fdTable2,&fileStatsTable2);
    size_t fileSizeTable2 = fileStatsTable2.st_size;

    long* posKmerTable1 = (long*) mmap(NULL, fileSizeTable1, PROT_READ,MAP_PRIVATE, fdTable1, 0);
    long* posKmerTable2 = (long*) mmap(NULL, fileSizeTable2, PROT_READ,MAP_PRIVATE, fdTable2, 0);
    if (posix_madvise (posKmerTable1, fileSizeTable1, POSIX_MADV_SEQUENTIAL|POSIX_MADV_WILLNEED) != 0){
        Debug(Debug::ERROR) << "posix_madvise returned an error for k-mer table 2 " << "\n";
    } 
    if (posix_madvise (posKmerTable2, fileSizeTable2, POSIX_MADV_SEQUENTIAL|POSIX_MADV_WILLNEED) != 0){
        Debug(Debug::ERROR)  << "posix_madvise returned an error for k-mer table 2 " << "\n";
    }

    long* currentQuerryPos = posKmerTable1;
    long* currentTargetPos = posKmerTable2;
    long* endTargetPos = posKmerTable2 + fileSizeTable2/sizeof(long);
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
    Debug(Debug::INFO)<<"number of equal Kmers: "<<equalKmers<<"\n";

    munmap(posKmerTable1, fileSizeTable1);
    munmap(posKmerTable2, fileSizeTable2);
    fclose(handleKmerTable1);
    fclose(handleKmerTable2);

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
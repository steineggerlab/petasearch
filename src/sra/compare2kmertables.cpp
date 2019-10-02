#include "LocalParameters.h"
#include "Command.h"
#include "Debug.h"
#include <sys/stat.h>
#include <sys/mman.h>
#include "MathUtil.h"
#include "QueryTableEntry.h"
#include "omptl/omptl_algorithm"


unsigned short maxshort = 65534;
QueryTableEntry *removeNotHittedSequences(QueryTableEntry * startPos, QueryTableEntry *endPos);
int resultTableSort(const QueryTableEntry &first, const QueryTableEntry &second);

int compare2kmertables(int argc, const char **argv, const Command& command){
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, 3, false);
    Timer timer;

    FILE* handleQueryKmerTable = fopen(par.db1.c_str(),"rb+");
    int fdQueryTable = fileno(handleQueryKmerTable);
    struct stat fileStatsQueryTable;
    fstat(fdQueryTable, &fileStatsQueryTable);
    size_t fileSizeQueryTable = fileStatsQueryTable.st_size;

    FILE* handleTargetKmerTable = fopen(par.db2.c_str(),"rb");
    int fdTargetTable = fileno(handleTargetKmerTable);
    struct stat fileStatsTargetTable;
    fstat(fdTargetTable,&fileStatsTargetTable);
    size_t fileSizeTargetTable = fileStatsTargetTable.st_size;

    FILE* handleTargetIDTable = fopen(par.db3.c_str(),"rb");
    int fdTargetIDTable = fileno(handleTargetIDTable);
    struct stat fileStatsTargetIDTable;
    fstat(fdTargetIDTable,&fileStatsTargetIDTable);
    size_t fileSizeTargetIDTable = fileStatsTargetIDTable.st_size;

    QueryTableEntry* startPosQueryTable = (QueryTableEntry*) mmap(NULL, fileSizeQueryTable, PROT_READ | PROT_WRITE,MAP_SHARED, fdQueryTable, 0);
    unsigned short* startPosTargetTable = (unsigned short*) mmap(NULL, fileSizeTargetTable, PROT_READ,MAP_PRIVATE, fdTargetTable, 0);
    unsigned int* startPosIDTable = ( unsigned int *) mmap(NULL,fileSizeTargetIDTable,PROT_READ,MAP_PRIVATE,fdTargetIDTable,0); 
    if (posix_madvise (startPosQueryTable, fileSizeQueryTable, POSIX_MADV_SEQUENTIAL|POSIX_MADV_WILLNEED) != 0){
        Debug(Debug::ERROR) << "posix_madvise returned an error for the query k-mer table\n";
    } 
    if (posix_madvise (startPosTargetTable, fileSizeTargetTable, POSIX_MADV_SEQUENTIAL|POSIX_MADV_WILLNEED) != 0){
        Debug(Debug::ERROR)  << "posix_madvise returned an error for the target k-mer table\n";
    }
    if (posix_madvise (startPosIDTable, fileSizeTargetIDTable, POSIX_MADV_SEQUENTIAL|POSIX_MADV_WILLNEED) != 0){
        Debug(Debug::ERROR)  << "posix_madvise returned an error for the target id table\n";
    }

    QueryTableEntry* currentQueryPos = startPosQueryTable;
    QueryTableEntry* endQueryPos = startPosQueryTable + fileSizeQueryTable/sizeof(QueryTableEntry);
    unsigned short* currentTargetPos = startPosTargetTable;
    unsigned short* endTargetPos = startPosTargetTable + fileSizeTargetTable/sizeof(short);
    unsigned int* currentIDPos = startPosIDTable;
    size_t equalKmers = 0;
    size_t currentKmer =0;
    struct timeval startTime;
    struct timeval endTime; 
    gettimeofday(&startTime, NULL);
    
   

    // cover the rare case that the first (real) target entry is larger than maxshort

    while(currentTargetPos <= endTargetPos && *currentTargetPos == maxshort){
        currentKmer+=maxshort;
        ++currentTargetPos;
        ++currentIDPos;
    }
    currentKmer+= *currentTargetPos; 

    while(__builtin_expect(currentTargetPos < endTargetPos,1)){
        if(currentKmer == currentQueryPos->Query.kmer){
            ++equalKmers;
            currentQueryPos->targetSequenceID = *currentIDPos;
            ++currentQueryPos;
            while(currentQueryPos->Query.kmer == currentKmer && __builtin_expect(currentQueryPos<endQueryPos,1)){
                currentQueryPos->targetSequenceID = *currentIDPos;
                ++currentQueryPos;
            }
            ++currentTargetPos;
            ++currentIDPos;
            while(__builtin_expect(*currentTargetPos == maxshort,0)){
                currentKmer+= maxshort;
                ++currentTargetPos;
                ++currentIDPos;
            }
            currentKmer+=*currentTargetPos;
        }

         while ( currentQueryPos->Query.kmer < currentKmer && __builtin_expect(currentQueryPos<endQueryPos,1) ){
            ++currentQueryPos;
         }

         while(currentKmer < currentQueryPos->Query.kmer && currentTargetPos < endTargetPos){
            ++currentTargetPos;
            ++currentIDPos;
            while(__builtin_expect(*currentTargetPos == maxshort,0)){
                currentKmer+= maxshort;
                ++currentTargetPos;
                ++currentIDPos;
            }
            currentKmer+= *currentTargetPos;
         }      
    }


    gettimeofday(&endTime, NULL);
    double timediff = (endTime.tv_sec - startTime.tv_sec) + 1e-6 * (endTime.tv_usec - startTime.tv_usec);
    Debug(Debug::INFO) << timediff<<" s; Rate "<<((fileSizeTargetTable+fileSizeQueryTable+fileSizeTargetIDTable)/1e+9)/timediff
            <<" GB/s \n";
    Debug(Debug::INFO)<<"number of equal Kmers: "<<equalKmers<<"\n";

    Debug(Debug::INFO)<<"Sorting Resulttable "<<"\n";

    omptl::sort(startPosQueryTable,endQueryPos, resultTableSort);

    Debug(Debug::INFO) << "Removing Sequences with less than two hits \n";
    QueryTableEntry *truncatedResultEndPos = removeNotHittedSequences(startPosQueryTable,endQueryPos);
    Debug(Debug::INFO) << "Truncating resuly table: Removing " << endQueryPos-truncatedResultEndPos 
        << " entries with less than 2 hits. \n";
    if(ftruncate(fdQueryTable,(truncatedResultEndPos-startPosQueryTable) * sizeof(QueryTableEntry))){
        Debug(Debug::ERROR) << "An error occurred while truncating the file. It should be truncated to the size of "
            << (truncatedResultEndPos-startPosQueryTable) * sizeof(QueryTableEntry) << "bytes.\n";
    }


    munmap(startPosQueryTable, fileSizeQueryTable);
    munmap(startPosTargetTable, fileSizeTargetTable);
    munmap(startPosIDTable, fileSizeTargetIDTable);
    fclose(handleQueryKmerTable);
    fclose(handleTargetKmerTable);
    fclose(handleTargetIDTable);

    return 0;
}

QueryTableEntry *removeNotHittedSequences(QueryTableEntry *startPos, QueryTableEntry *endPos){
    QueryTableEntry *currentReadPos = startPos;
    QueryTableEntry *currentWritePos = startPos;
    while(currentReadPos < endPos){
        size_t count = 0;
        while(currentReadPos < endPos && (currentReadPos + 1)->targetSequenceID == currentReadPos->targetSequenceID){
            count++;
            ++currentReadPos;
        }
        if(count > 1){
            memcpy(currentWritePos,currentReadPos - count,sizeof(QueryTableEntry) * count );
            currentWritePos+=count;
            //memcopy the sequences
        }
        ++currentReadPos;
    }
    return currentWritePos;
}

int resultTableSort(const QueryTableEntry &first, const QueryTableEntry &second){   
    if(first.targetSequenceID < second.targetSequenceID)
        return true;
    if(second.targetSequenceID < first.targetSequenceID)
        return false;
    if (first.querySequenceId > second.querySequenceId)
        return true;
    if (second.querySequenceId > first.querySequenceId)
        return false;
    if (first.Query.kmerPosInQuery < second.Query.kmerPosInQuery)
        return true;
    if (second.Query.kmerPosInQuery < first.Query.kmerPosInQuery)
        return false;
    if (first.Query.kmer < second.Query.kmer)
        return true;
    if (second.Query.kmer < first.Query.kmer)
        return false;
    return false;
}
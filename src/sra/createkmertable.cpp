#include "Parameters.h"
#include "Command.h"
#include "Debug.h"

#define KMER_SIZE 5
#define SPACED_KMER true


int createkmertable(int argc, const char **argv, const Command& command){
    Parameters& par = Parameters::getInstance();
    par.kmerSize = KMER_SIZE;
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, 1, false);

    Debug(Debug::INFO) << par.db1 << "\n";
    
    return EXIT_SUCCESS;
}
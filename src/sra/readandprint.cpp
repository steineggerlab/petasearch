//
// Created by matchy on 1/24/22.
//
#include "LocalParameters.h"
#include "SRADBReader.h"

int readandprint(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, true, 0, LocalParameters::PARSE_VARIADIC);

    SRADBReader reader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                       DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    Debug(Debug::INFO) << "DB size: " << reader.getSize() << "\n";

    for (size_t i = 0; i < reader.getSize(); ++i) {
        char *data = reader.getData(i, 0);
        unsigned int seqLen = reader.getSeqLen(i);
        Debug(Debug::INFO) << "sequence " << i << ": " << data << "\n";
        Debug(Debug::INFO) << "length: " << seqLen << "\n";
    }

    return EXIT_SUCCESS;
}

#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const char* binary_name = "srasearch";
const char* tool_name = "srasearch";
const char* tool_introduction = "SRA Search.";
const char* main_author = "Jonas Huegel";
const char* show_extended_help = "1";
const char* show_bash_info = "1";
bool hide_base_commands = true;
void updateValidation();
void (*validatorUpdate)(void) = updateValidation;

LocalParameters& localPar = LocalParameters::getLocalInstance();

std::vector<struct Command> commands = {
        {
            "petasearch", petasearch, &localPar.petasearchworkflow, COMMAND_MAIN,
             "",
             NULL,
             "Jonas Hügel <jonas.huegel@mpibpc.mpg.de> ",
             "<i:queryDB> <i:targetDBsFile> <i:compResultDB> <o:alignmentResults> <o:tmp>",
             CITATION_MMSEQS2,
             {{"queryDB",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
              {"targetDBsFile",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
              {"compResultDBsFile",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
              {"alignmentResults",DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile},
              {"tmp",DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}
        },
        {
            "easy-petasearch", easypetasearch, &localPar.easypetasearchworkflow, COMMAND_MAIN,
             "",
             NULL,
             "Jonas Hügel <jonas.huegel@mpibpc.mpg.de> ",
             "<i:queryDB> <i:targetDB> <i:resultDB> <o:tmp>",
             CITATION_MMSEQS2,
             {{"queryDb",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
              {"targetDb",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
              {"resultDb",DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile},
              {"tmp",DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}
        },
        {
            "createkmertable", createkmertable, &localPar.createkmertable, COMMAND_EXPERT,
            "Extracts a table containing all unique k-mers",
            "Extracts a unique k-mer table.",
            "Jonas Hügel <jonas.huegel@mpibpc.mpg.de> ",
            "<i:sequenceDB> <o:kmerTable>",
            CITATION_MMSEQS2,
            {{"sequenceDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
             {"kmerTable",DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}}
        },
        {
            "compare2kmertables", compare2kmertables, &localPar.compare2kmertables, COMMAND_EXPERT,
            "",
            NULL,
            "Jonas Hügel <jonas.huegel@mpibpc.mpg.de> ",
            "<i:queryKmerTable> <i:targetKmerTable> <o:resultTable>",
            CITATION_MMSEQS2,
            {{"queryDb",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
             {"targetKmerTable",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
             {"resultTable",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile}}
        },
        {
            "computeAlignments", computeAlignments, &localPar.computeAlignments, COMMAND_EXPERT,
            "",
            NULL,
            "Jonas Hügel <jonas.huegel@mpibpc.mpg.de> ",
            "<i:querySequenceDB> <i:targetSequenceDB> <i:prevResultTable> <o:alignmentFile>",
            CITATION_MMSEQS2,
            {{"querySequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
             {"targetSequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
             {"prevResultTable", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb},
             {"alignmentFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}}
        }
};

void updateValidation() {}

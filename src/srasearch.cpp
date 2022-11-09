#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"
#include "DownloadDatabase.h"
#include "Prefiltering.h"

const char* binary_name = "srasearch";
const char* tool_name = "srasearch";
const char* tool_introduction = "SRA Search.";
const char* main_author = "Jonas Huegel";
const char* show_extended_help = "1";
const char* show_bash_info = "1";
const char* index_version_compatible = "16";

bool hide_base_commands = true;
bool hide_base_downloads = false;

void updateValidation();
void (*validatorUpdate)(void) = updateValidation;

std::vector<DatabaseDownload> externalDownloads = {};
std::vector<KmerThreshold> externalThreshold = {};


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
            "comparekmertables", comparekmertables, &localPar.comparekmertables, COMMAND_EXPERT,
            "",
            NULL,
            "Jonas Hügel <jonas.huegel@mpibpc.mpg.de> ",
            "<i:queryKmerTable> <i:targetKmerTable> <o:resultTable>",
            CITATION_MMSEQS2,
            {{"queryDb",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
             {"targetKmerTable",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
             {"resultTable",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile}}
        },
        {   // bc compatibility
            "compare2kmertables", comparekmertables, &localPar.comparekmertables, COMMAND_HIDDEN, "", NULL, "", "", 0, {}
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
        },
        {
                "convert2sradb", convert2sradb, &localPar.convert2sradb, COMMAND_EXPERT,
                "",
                NULL,
                "Minghang Li <matchy@snu.ac.kr> ",
                "<i:sequenceDB>|<i:fastaFiles> <o:srasearchDB>",
                CITATION_MMSEQS2,
                {{"dbToConvert", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                 {"outputDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}}
        },
        {
            "convertsraalis", convertsraalignments, &localPar.convertsraalignments, COMMAND_EXPERT,
            "",
            NULL,
            "Minghang Li <matchy@snu.ac.kr> ",
            "<i:queryDB> <i:targetSRADB> <i:sraAlignmentFile> <o:sraAlignmentDB>",
            CITATION_MMSEQS2,
            {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
             {"targetSRADB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
             {"sraAlignmentFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
             {"sraAlignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}}
        },
        {
            "readandprint", readandprint, &localPar.readandprint, COMMAND_EXPERT,
            "",
            NULL,
            "Minghang Li <matchy@snu.ac.kr> ",
            "<i:sraDB>",
            CITATION_MMSEQS2,
            {{"inputdb", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile }}
        }
};

void updateValidation() {}

#include "LocalParameters.h"
#include "FileUtil.h"
#include "Util.h"
#include "CommandCaller.h"

#include "petasearch.sh.h"

#include <cassert>

void setPetaSearchWorkflowDefaults(Parameters *p) {
    p->kmerSize = 9;
    p->spacedKmer = false;
    p->evalThr = 1000;
}

int petasearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setPetaSearchWorkflowDefaults(&par);

    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string tmpDir = par.db5;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;



    cmd.addVariable("CREATE_TTABLE_PAR", par.createParameterString(par.createkmertable).c_str());
    cmd.addVariable("COMP_KMER_TABLES_PAR", par.createParameterString(par.compare2kmertables).c_str());
    cmd.addVariable("COMP_ALI_PAR", par.createParameterString(par.computeAlignments).c_str());
    par.evalThr = 100000;
    cmd.addVariable("SWAP_PAR", par.createParameterString(par.swapresult).c_str());
    cmd.addVariable("CONVERTALIS_PAR",par.createParameterString(par.convertalignments).c_str());

    std::string program = tmpDir + "/petasearch.sh";
    FileUtil::writeFile(program, petasearch_sh, petasearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
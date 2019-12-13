#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "easypetasearch.sh.h"

#include <cassert>

void setEasyPetaSearchWorkflowDefaults(Parameters *p) {
    p->kmerSize = 9;
    p->spacedKmer = false;
    p->evalThr = 1000;
}

int easypetasearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setEasyPetaSearchWorkflowDefaults(&par);

    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(par.filenames, par.easypetasearchworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
//    cmd.addVariable("RUNNER", par.runner.c_str());

    par.createTargetTable = 0;
    cmd.addVariable("CREATE_QTABLE_PAR", par.createParameterString(par.createkmertable).c_str());
    par.createTargetTable = 1;
    cmd.addVariable("CREATE_TTABLE_PAR", par.createParameterString(par.createkmertable).c_str());
    cmd.addVariable("COMP_KMER_TABLES_PAR", par.createParameterString(par.compare2kmertables).c_str());
    cmd.addVariable("COMP_ALI_PAR", par.createParameterString(par.computeAlignments).c_str());
    par.evalThr = 100000;
    cmd.addVariable("SWAP_PAR", par.createParameterString(par.swapresult).c_str());
    cmd.addVariable("CONV_PAR",par.createParameterString(par.convertalignments).c_str());

    std::string program = tmpDir + "/easypetasearch.sh";
    FileUtil::writeFile(program, easypetasearch_sh, easypetasearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
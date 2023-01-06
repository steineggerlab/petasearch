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
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("PETASEARCH_PAR", par.createParameterString(par.petasearchworkflow).c_str());
    cmd.addVariable("CONVERTALIS_PAR",par.createParameterString(par.convertalignments).c_str());

    std::string program = tmpDir + "/easypetasearch.sh";
    FileUtil::writeFile(program, easypetasearch_sh, easypetasearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}

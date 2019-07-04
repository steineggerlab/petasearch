#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

class LocalParameters : public Parameters {
public:
    static void initInstance() {
        new LocalParameters;
    }

    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initInstance();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }

    std::vector<MMseqsParameter*> createkmertable;
    std::vector<MMseqsParameter*> compare2kmertables;
    ///example

    std::vector<MMseqsParameter*> predictexonsworkflow;

    PARAMETER(PARAM_REVERSE_FRAGMENTS)
    int reverseFragments;
    //example end

    PARAMETER(PARAM_CREATE_TARGET_TABLE)
    int createTargetTable;

private:
    LocalParameters() : 
        Parameters(),

      PARAM_REVERSE_FRAGMENTS(PARAM_REVERSE_FRAGMENTS_ID,"--reverse-fragments", "Reverse AA Fragments", "reverse AA fragments to compute under null [0,1]", typeid(int), (void *) &reverseFragments, "^[0-1]{1}$"),
      PARAM_CREATE_TARGET_TABLE(PARAM_CREATE_TARGET_TABLE_ID,"--createTargetTable","creating target table ?","create target table (1) or query table (0) (default: 1) [0,1]",typeid(int), (void *) &createTargetTable,"^[0-1]{1}$")
    {
        predictexonsworkflow.push_back(&PARAM_REVERSE_FRAGMENTS);
        createkmertable.push_back(&PARAM_K);
        createkmertable.push_back(&PARAM_THREADS);
        createkmertable.push_back(&PARAM_V);
        createkmertable.push_back(&PARAM_CREATE_TARGET_TABLE);

        // default value 0 means no reverse of AA fragments
        reverseFragments = 0;
        createTargetTable = 0;
    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);

    

};

#endif

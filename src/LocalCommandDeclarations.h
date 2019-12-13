#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int createkmertable(int argc, const char **argv, const Command& command);
extern int compare2kmertables(int argc, const char **argv, const Command& command);
extern int computeAlignments(int argc, const char **argv, const Command& command);
extern int petasearch(int argc, const char **argv, const Command& command);
extern int easypetasearch(int argc, const char **argv, const Command &command);
#endif

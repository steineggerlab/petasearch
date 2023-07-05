#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int createkmertable(int argc, const char **argv, const Command& command);
extern int comparekmertables(int argc, const char **argv, const Command& command);
extern int blockalign(int argc, const char **argv, const Command& command);
extern int convert2sradb(int argc, const char **argv, const Command& command);
extern int petasearch(int argc, const char **argv, const Command& command);
extern int easypetasearch(int argc, const char **argv, const Command &command);
extern int convertsraalignments(int argc, const char **argv, const Command &command);
extern int readandprint(int argc, const char **argv, const Command &command);
extern int playground(int argc, const char **argv, const Command &command);

#endif

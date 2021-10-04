
#ifndef MAIN_H_
#define MAIN_H_

#include <fstream>

#ifndef _MAXBUFFSIZE
#define _MAXBUFFSIZE 256
#endif

#ifndef INPUTPARSER_HPP_
class InputFile;
#endif

#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef STRUCTURE_HPP_
class Structure;
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

/*
 * If program is not run properly or arg -h was called this
 * function prints a minimal documentation on shell
 */
void outputHelp();

void shutUp(BasicInformation &);

void setupMolecule(Molecule &, Flags &flags);

int
ReadInput(Molecule &, InputFile &, BasicInformation &, Flags &, std::fstream &);

int initializeProgram(BasicInformation &, Flags &, int, char *[]);

/* Analyzes the arguments types at programm start */
int checkArguments(const char *);

/* Calls CheckParameters and saves the respective arguments */
int initializeArguments(int, char *[], BasicInformation &, Flags &);
int grumble(BasicInformation &, Flags &);
void handleMemory(BasicInformation &);
void freeMemory(BasicInformation &);
int handleOutput(BasicInformation &, Flags &);
void killTITANIA(BasicInformation &, Flags &);

void handle_stop(Structure *, BasicInformation &, Flags &);
/* Initializes the default setting of the programm */
void baseSettings(BasicInformation &, Flags &);

#endif

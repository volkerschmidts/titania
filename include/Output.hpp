
#ifndef OUTPUT_HPP_
#define OUTPUT_HPP_

#include <cstdio>
#include <eigen3/Eigen/Core>
#include <fstream>
#include <string>

#ifndef OUTPUT_OPTIONS_TITANIA_
#include <OutputOptions.hpp>
#endif

#ifndef ATOM_HPP_
class Atom;
#endif

#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef RDCSET_HPP_
class RDCset;
#endif

#ifndef STRUCTURE_HPP_
class Structure;
#endif

#ifndef INPUTPARSER_HPP_
class InputFile;
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

#ifndef ALIGNMENT_OUTPUT_TITANIA_
struct AlignmentOutput;
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
enum struct StructureOptions : unsigned int;
#endif

#ifndef REDUNDANTINTERNALS_HPP_
class RedundantInternals;
#endif

#define OUTPUT_INFOBOX_STRING_SIZE 56

const char *to_cstring(double, int, int, bool exponential = false);
const char *to_cstring(int, int);
const char *to_cstring(unsigned, int);
const char *to_cstring(long int, int);
const char *to_cstring(long unsigned, int);

//#define ENABLE_OLD_OUTPUT_

void shutup(BasicInformation &);

void makeHeader(FILE *, BasicInformation &, Flags &, InputFile &);
void makeReferences(FILE *, BasicInformation &);
void makeInputInformation(FILE *, BasicInformation &, Flags &, Molecule &);
void makeTITANIAresults(FILE *,
                        BasicInformation &,
                        Flags &,
                        Molecule &,
                        Structure *);
void makeRuntimeInformation(FILE *, BasicInformation &, Flags &, Molecule &);

#ifdef ENABLE_OLD_OUTPUT_
#include <fstream>
/* Echos the input file to output file */
void echoInput(InputFile &, std::fstream &);
void outputCommentary(InputFile &, std::fstream &);

void generateOutput(Molecule &, BasicInformation &, Flags &);
void output_Angles(Molecule &, Structure *, BasicInformation &, Flags &);
void output_chiralVolumes(Molecule &, BasicInformation &, Flags &);
#endif
void generateAlignmentOutput(Molecule &,
                             Structure *,
                             RDCset *,
                             AlignmentOutput *,
                             BasicInformation &,
                             Flags &);
#ifdef ENABLE_OLD_OUTPUT_
void output_SECONDA(Molecule &, BasicInformation &);
void output_phobos(Molecule &, BasicInformation &);
void output_Q_factor(Molecule &, BasicInformation &, Flags &);
void output_runtimeInformation(BasicInformation &, Flags &, Molecule &);
void output_rdc_motion(Structure *, BasicInformation &, AlignmentOutput *);
#endif
void output_big_data(Molecule &, BasicInformation &, Flags &);

void outputXYZ(Atom *,
               BasicInformation &,
               StructureOptions,
               std::string,
               bool = false);
void
titania2hotFCHT(Molecule &, Structure *, RDCset *, BasicInformation &, Flags &);
void titania2titania(Structure *, BasicInformation &);
// void titania2titania ( Molecule&, Structure*, InputFile&, BasicInformation&
// );

void outputMatrix(Eigen::MatrixXd,
                  unsigned int,
                  unsigned int,
                  std::fstream &,
                  OutputOptions opt = OutputOptions::undefined);
void outputVector(Eigen::VectorXd,
                  unsigned int,
                  unsigned int,
                  std::fstream &,
                  OutputOptions opt = OutputOptions::undefined,
                  bool trans        = false);
std::string centerString(std::string, int);
const char *centerCstring(const char *, int);
std::string ValueError2String(double, double, int, int);

void
prepare_trj_file(Molecule &, RedundantInternals &, BasicInformation &, Flags &);
void outputStructure2trj(BasicInformation &,
                         Flags &,
                         Molecule &,
                         Structure *,
                         RedundantInternals &,
                         int);

#endif

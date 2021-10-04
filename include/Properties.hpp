
#ifndef PROPERTIES_HPP_
#define PROPERTIES_HPP_

#include <string>

#ifndef ATOM_HPP_
class Atom;
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

/* Nucleus information of all isotopes imported from hotFCHT */
struct NucleusList {
  int z;                   /* Atomic number */
  int A;                   /* Isotope */
  double mass;             /* Atomic mass of respective isotope */
  double abundance;        /* Abundance of respective isotope */
  double nuclearspin;      /* Nuclear spin of respective isotope */
  double magneticmoment;   /* Magnetic moment of respective isotope */
  double quadrupolemoment; /* Quadrupolar moment of respective isotope */
};

/* List of all Elements with all isotopes and nmr default isotope */
struct ElementList {
  int z;                 /* Atomic number */
  std::string element;   /* String representation of element */
  int isotopes[255];     /* List of all isotopes */
  double mass;           /* Abundance averaged mass */
  int rdcdefaultnucleus; /* Default isotope for rdcs... */
};

/* Checks if String represents a proper element */
bool isElement(std::string);

/* Checks for String representation in ElementList */
int findNucleusByName(Atom &);

/* Checks for atomic number in ElementList */
int findNucleusByZ(Atom &);

/* Saves respective NucleusList in Atom */
int loadNucleus(Atom &, const int, const int);

/* Calculates Dmax for user defined normalization factor */
int loadRDCNormFact(Atom &, Atom &, BasicInformation &, Flags &);

unsigned int getStandardBonds(Atom *);
void loadHybridisation(Atom *);

#endif

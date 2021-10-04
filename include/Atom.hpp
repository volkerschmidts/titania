
#ifndef ATOM_HPP_
#define ATOM_HPP_

#include <eigen3/Eigen/Core>
#include <string>

#define CHIRAL_CENTER_ATOMS_ 3

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

#ifndef STRUCTURE_HPP_
class Structure;
#endif

#ifndef SPHERICALHARMONICS_HPP_
class SphericalHarmonics;
#endif

#ifndef BOND_HPP_
class Bond;
#endif

#ifndef STRUCTURE_SIMULATOR_HPP_
class StructureSimulator;
#endif

#ifndef POTENTIAL_HPP_
class Potential;
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
enum struct StructureOptions : unsigned int;
#endif

#ifndef HYBRIDISATION_TITANIA_
enum struct Hybridisation : unsigned int;
#endif

#ifndef COORDINATES_TITANIA_
struct Coordinates;
#endif

#ifndef PROPERTIES_HPP_
struct NucleusList;
#endif

class Atom {
 private:
  std::string element; /* Element symbol */
  std::string label;   /* Unique label */
  std::string chiralVolumeAtoms;
  int z;        /* Atomic number */
  int A;        /* Isotope */
  int AtomType; /* Atom type based on mm3 */
  int distance_violations;
  int chiral_index[CHIRAL_CENTER_ATOMS_];
  unsigned int inputIndex; /* Position in input file */
  unsigned int rdcIndex;   /* Number of RDCs adapted numbering scheme */
  unsigned int
      NumberOfBonds; /* Number of bonds (double/triple bonds = 1 bond) */
  unsigned int NumberOfBondRDCs;  /* Number of 1D-RDCs */
  unsigned int NumberOfHarmonics; /* Number of nD-RDCs */
  double gamma;                   /* Gyromagnetic ratio */
  double mass;                /* Abundance averaged mass (for faster access) */
  double partialCharge;       /* Partial charge based on mmff94 bci */
  double Energy;              /* Current atom energy */
  double initialChiralVolume; /* V(chiral) = r(a)^t * ( r(b)xr(b) ) */
  double optimizedChiralVolume;
  double vdW_Potential;
  bool fixed; /* Used for Structure::toxyz to check if optCoord are known */
  bool invertable;
  bool increased_violations;
  NucleusList *nucleus;     /* Information of isotope */
  Coordinates *coordinates; /* Obvious? */
  Coordinates
      *optimized; /* Coordinates obtained after structure optimization */
  Hybridisation hybr;
  Structure *parent; /* Structure model in which atom is defined */
  Atom *nextAtom;    /* Linked */
  Atom *prevAtom;    /* List */
  SphericalHarmonics *
      *harmonics; /* Lists all spherical harmonics describing rdcs */
  Bond **bonds;

  // Add some additional function for faster access within this class
  void setOC(double, double, double);
  void setIC(double, double, double);
  Eigen::Vector3d getOC();
  Eigen::Vector3d getIC();

 public:
  /* Default construcor. Used to define tmp Atoms or HeadAtom */
  Atom();
  /* Copy (like) constructor. Used to append atoms in a structure */
  Atom(Atom *a);
  /* Default destructor. */
  ~Atom();

  /* Method to initialize all atoms to new Structure using Atom(Atom*) */
  Atom *initializeOptimizedAtoms(Structure *, Structure *);

  /* Setting next AtomNode */
  Atom *setNext();

  /* Accessing next/previous AtomNodes */
  Atom *getNext() const
  {
    return nextAtom;
  };
  Atom *getPrev() const
  {
    return prevAtom;
  };

  /* InputIndex depending on input file 1 - n */
  void setIndex(unsigned int index)
  {
    inputIndex = index;
  }
  unsigned int getIndex() const
  {
    return inputIndex;
  }
  void setrdcIndex(unsigned int index)
  {
    rdcIndex = index;
  }
  unsigned int getrdcIndex() const
  {
    return rdcIndex;
  }

  /* Element defined as string according to periodic table */
  void setElement(std::string nEl)
  {
    element = nEl;
  }
  std::string getElement() const
  {
    return element;
  }

  /* Labels defined in input file by user */
  void setLabel(std::string nLabel)
  {
    label = nLabel;
  }
  std::string getLabel() const
  {
    return label;
  }
  std::string getIdentifier() const
  {
    return (element + label);
  }

  /* int representation of atomic number (mainly handled by Properties.cpp)*/
  void setZ(int nZ)
  {
    z = nZ;
  }
  int getZ() const
  {
    return z;
  }

  /* int representation of isotop (if standard: Properties.cpp, else user
   * defined) */
  void setA(int nA)
  {
    A = nA;
  }
  int getA() const
  {
    return A;
  }

  /* Includes all information about nucleus (handled by Properties.cpp ) */
  void setNucleus(NucleusList *nNuc)
  {
    nucleus = nNuc;
  }
  void setHybridisation(Hybridisation h)
  {
    hybr = h;
  }
  Hybridisation &getHybridisation()
  {
    return hybr;
  }
  void loadAtomType(Flags &flags);
  int getAtomType() const
  {
    return AtomType;
  }
  void loadBondorder();

  /* Coordinates of input (previous Structure) or calculated via toxyz (after
   * optimization) */
  void setCoordinates(double, double, double, StructureOptions);
  void setCoordinates(Eigen::Vector3d, StructureOptions);
  void setCoordinates(double *, StructureOptions);

  Coordinates *getCoordinates(StructureOptions) const; // { return coordinates;}
  Eigen::Vector3d Coordinates2Eigen(StructureOptions);
  void shiftCoordinates(Eigen::Vector3d &, StructureOptions);
  void rotateCoordinates(Eigen::Matrix3d &, StructureOptions);
  void retainCoordinates();
  double getDistance(Atom *, StructureOptions);
  double getAngle(Atom *, Atom *, StructureOptions);

  /* Parent defines the class, which saves the current class (here Atom saved in
   * Structure) */
  void setParent(Structure *nParent)
  {
    parent = nParent;
  }
  Structure *getParent() const
  {
    return parent;
  };

  /* setGamma() loads gamma from Properties.cpp */
  void setGamma();
  double getGamma() const
  {
    return gamma;
  }

  /* Adds one bondpartner based on ConnectivityMap (see ReadInput.cpp) */
  void initializeBonds();
  void addBondpartner(Atom *);
  void addBond(Bond *);
  void setBond(Bond *b, unsigned int i)
  {
    bonds[i] = b;
  }
  Bond *getBond(unsigned int i) const
  {
    if (i < NumberOfBonds)
      return bonds[i];
    else
      return NULL;
  }
  Bond *getBond(int i) const
  {
    if (i >= 0 && (unsigned int) i < NumberOfBonds)
      return bonds[i];
    else
      return NULL;
  }
  Bond *getBond(Atom *) const;
  Atom *getBondpartner(unsigned int) const;
  unsigned int getNumberOfBonds() const
  {
    return NumberOfBonds;
  }
  void increaseNumberOfBonds()
  {
    ++NumberOfBonds;
  }
  void setNumberOfBonds(unsigned int b)
  {
    NumberOfBonds = b;
  }
  void set_chiral_volume_indizes();
  bool hasChiralVolume() const
  {
    if (NumberOfBonds >= 3 && !invertable)
      return true;
    else
      return false;
  };
  bool is_invertable() const
  {
    return invertable;
  }
  double getChiralVolume(StructureOptions, bool);
  std::string printChiralVolumeAtoms();
  const char *printChiralVolumeAtoms(int);
  double check_Chiral_Volume_validity(StructureOptions);

  /* Adds spherical harmonics based on defined rdcs */
  void initializeHarmonics(); // { bondHarmonics = (SphericalHarmonics**) malloc
                              // ( NumberOfHarmonics * sizeof (
                              // SphericalHarmonics* ) ); }
  void clearHarmonics();
  void addBondharmonic(SphericalHarmonics *);
  SphericalHarmonics *getBondharmonic(unsigned int i) const;
  SphericalHarmonics *getHarmonic(unsigned int i) const;
  SphericalHarmonics *getHarmonic(Atom *) const;
  unsigned int getNumberOfHarmonics() const
  {
    return NumberOfHarmonics;
  }
  void setNumberOfHarmonics(unsigned int i)
  {
    NumberOfHarmonics = i;
  }
  void increaseNumberOfHarmonics()
  {
    ++NumberOfHarmonics;
  }

  /* Accessing internal booleans used in Structure::toxyz() */
  bool isFixed() const
  {
    return fixed;
  };
  void setFixed()
  {
    fixed = true;
  };
  void setUnfixed()
  {
    fixed = false;
  };

  /* Faster access of atomic mass */
  double getMass() const
  {
    return mass;
  };
  void setMass(double m)
  {
    mass = m;
  };

  void setPartialCharge(double p)
  {
    partialCharge = p;
  }
  double getPartialCharge() const
  {
    return partialCharge;
  }
  void setEnergy(double E)
  {
    Energy = E;
  }
  double getEnergy() const
  {
    return Energy;
  }

  double calculateVdW_Potential(StructureSimulator &, StructureOptions);
  double getVdW_Potential();
  int get_number_of_distance_violations() const
  {
    return distance_violations;
  }
  void increase_distance_violations()
  {
    if (!increased_violations)
    {
      ++distance_violations;
      increased_violations = true;
    }
  }
  void reset_distance_violations()
  {
    distance_violations = 0;
  }
};

void atomByIdentifier(std::string, Atom **);

void atomByrdcIndex(Atom **, unsigned int);

#endif

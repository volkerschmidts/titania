
#ifndef MOLECULE_HPP_
#define MOLECULE_HPP_

#include <eigen3/Eigen/Core>
#include <string>

#ifndef ATOM_HPP_
class Atom;
#endif

#ifndef STRUCTURE_HPP_
class Structure;
#endif

#ifndef RDCDATA_HPP_
class RDCdata;
#endif

#ifndef RDCSET_HPP_
class RDCset;
#endif

#ifndef SPHERICALHARMONICSMATRIX_HPP_
class SphericalHarmonicsMatrix;
#endif

#ifndef INPUTPARSER_HPP_
class InputFile;
class StructureInput;
class RDCinput;
class KeywordInput;
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

class Molecule {
 private:
  std::string label;    /* Saves a label to identify Molecule */
  Structure *HeadStruc; /* First structure (input) of this molecule */
  Structure
      *TailStruc;  /* Last structure of this molecule (MC or opt Structure) */
  RDCset *HeadSet; /* First rdc set obtained from input */
  RDCset *TailSet; /* Last rdc set obtained from input */
  unsigned int NumberOfAtoms; /* Number of atoms defined by structure */
  unsigned int NumberOfRDCs;  /* Number of rdcs per set defined in input */
  unsigned int NumberOfRDCatoms;
  unsigned int NumberOfRDCsets; /* Number of rdc sets defined in input */
  unsigned int NumberOfStructs; /* Saves the number of structures. */
  unsigned int NumberOfLongRangeRDCs;
  unsigned int DBE; /* Saved double bond equivalents of the input structure */
  unsigned int NumberOfDoubleBonds;
  double
      MolecularMass; /* Molecular mass of molecule (same for all Structures) */
  Eigen::MatrixXd rdcMatrix; /* Saves the input rdc matrix */
  Eigen::MatrixXd inputB;    /* Saves the current B matrix */
  Eigen::MatrixXd kappa;
  Eigen::MatrixXd weights;
  Eigen::MatrixXd input_Transformation;
  Eigen::Matrix3d input_InertiaTensor;
  Eigen::Vector3d input_InertiaTensorEigVal;
  Eigen::Matrix3d input_InertiaTensorEigVec;

 public:
  /* Default constructor */
  Molecule(std::string newLabel);

  /* Default destructor */
  ~Molecule();

  /* Label defined in input file by user */
  void setLabel(const std::string newLabel)
  {
    label = newLabel;
  };
  std::string getLabel() const
  {
    return label;
  };

  int loadInput(InputFile &, BasicInformation &, Flags &);
  int loadStructure(InputFile &, StructureInput *, BasicInformation &);
  int loadRDCs(InputFile &, RDCinput *, BasicInformation &);
  int loadKeywords(KeywordInput *, BasicInformation &, Flags &);


  /* Number of Atoms defined in input file */
  unsigned int getNOA() const
  {
    return NumberOfAtoms;
  };
  void setNORA(const unsigned int nora)
  {
    NumberOfRDCatoms = nora;
  }
  unsigned int getNORA() const
  {
    return NumberOfRDCatoms;
  }

  /* Number of rdcs defined per set */
  void raiseNOR()
  {
    ++NumberOfRDCs;
  };
  unsigned int getNOR() const
  {
    return NumberOfRDCs;
  };

  /* Number of rdc sets. setNOS should not be used, better use raise... after
   * every new set! */
  unsigned int getNORsets() const
  {
    return NumberOfRDCsets;
  };
  unsigned int raiseNumberOfRDCsets()
  {
    return ++NumberOfRDCsets;
  };

  /* Access the rdc sets */
  RDCset *addHeadSet();
  RDCset *getHeadSet() const
  {
    return HeadSet;
  };
  RDCset *getTailSet() const
  {
    return TailSet;
  };
  void setTailSet(RDCset *t)
  {
    TailSet = t;
  };
  void setHeadSet(RDCset *t)
  {
    HeadSet = t;
  };
  int linkSets();

  /* Number of Structures */
  unsigned int raiseNumberOfStrucs()
  {
    return ++NumberOfStructs;
  };
  unsigned int getNumberOfStructures() const
  {
    return NumberOfStructs;
  };

  /* Access Structures (including optimized ones) */
  Structure *addHeadStruc();
  Structure *getHeadStruc() const
  {
    return HeadStruc;
  };
  Structure *getTailStruc() const
  {
    return (TailStruc ? TailStruc : HeadStruc);
  };
  void setTailStruc(Structure *s)
  {
    TailStruc = s;
  };
  void removeStructures(int);
  void checkPlanarity(BasicInformation &);

  /* Saves the rdc matrix (just input matrix) */
  void setRDCMatrix(Eigen::MatrixXd newM)
  {
    rdcMatrix = newM;
  };
  Eigen::MatrixXd getRDCmatrix() const
  {
    return rdcMatrix;
  };
  Eigen::MatrixXd getKappaMatrix();
  void initializeWmatrix();
  Eigen::MatrixXd getWmatrix();
  void setInputInertiaTensor(Eigen::Matrix3d);
  void setInputInertiaTensorEigVal(Eigen::Vector3d);
  void setInputInertiaTensorEigVec(Eigen::Matrix3d);
  void setInputTransformation(Eigen::Matrix3d);
  void setInputTransformation(Eigen::Vector3d);
  Eigen::MatrixXd getInputTransformation() const
  {
    return input_Transformation;
  }
  Eigen::Matrix3d getInputInertiaTensor() const
  {
    return input_InertiaTensor;
  }
  Eigen::Vector3d getInputInertiaTensorEval() const
  {
    return input_InertiaTensorEigVal;
  }
  Eigen::Matrix3d getInputInertiaTensorEvec() const
  {
    return input_InertiaTensorEigVec;
  }

  /* Access B (scaled Cosine) matrix */
  void setInputBMatrix(Structure *);
  Eigen::MatrixXd getInputBMatrix() const
  {
    return inputB;
  };

  /* Molecular mass handling */
  double determineMolecularMass();
  double getMolecularMass() const
  {
    return MolecularMass;
  };

  /* Hanling of DBEs & rings */
  void countLongRangeRDCs();
  unsigned int getNumberOfLongRangeRDCs() const
  {
    return NumberOfLongRangeRDCs;
  }
  void raiseNumberOfDB()
  {
    ++NumberOfDoubleBonds;
  }
  unsigned int getNumberOfDBs() const
  {
    return NumberOfDoubleBonds;
  }
};

int atomLabelParser(std::string, Atom *);
int initializeChiralVolumes(Molecule &CurrMol,
                            BasicInformation &baseInformation);

#endif

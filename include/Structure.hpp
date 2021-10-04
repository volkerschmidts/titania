
#ifndef STRUCTURE_HPP_
#define STRUCTURE_HPP_

#include <eigen3/Eigen/Core>
#include <string>

#ifndef ATOM_HPP_
class Atom;
#endif

#ifndef BOND_HPP_
class Bond;
#endif

#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef RDCDATA_HPP_
class RDCdata;
#endif

#ifndef SPHERICALHARMONICSMATRIX_HPP_
class SphericalHarmonicsMatrix;
#endif

#ifndef SPHERICALHARMONICS_HPP_
class SphericalHarmonics;
#endif

#ifndef STRUCTURE_SIMULATOR_HPP_
class StructureSimulator;
#endif

#ifndef LINALG_HPP_
class Tensor;
#endif

#ifndef MONTECARLOSTATISTIC_HPP_
class MonteCarloStatistic;
#endif

#ifndef MC_OUTPUT_TITANIA_
#include <MC_Output.hpp>
#endif

#ifndef RDCMATRIXOPTIONS_TITANIA_
enum struct rdcMatrixOptions : unsigned int;
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

#ifndef MC_STOP_TITANIA_
enum struct MC_stop : unsigned int;
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
#include <StructureOptions.hpp>
#endif

#ifndef REDUNDANTINTERNALS_HPP_
class RedundantInternals;
#endif

class Structure {
 private:
  unsigned int inputIndex; /* (Input) index of current structure -> faster
                              identification/access */
  unsigned int NumberOfBonds;
  unsigned int performedRedundantSteps;
  bool optimized; /* Internal tool to check if this structure results from SCRM
                     optimization */
  bool hasUndefinedRDCs;
  double radiusOfGyration;
  double rdc_rmsd;
  std::string label;    /* Sting to identify current structure. Automatically
                           generated for opt structures */
  Molecule *parent;     /* Molecule which defines this structure */
  Structure *nextStruc; /* Next structure node */
  Structure *prevStruc; /* Previous structure node */
  Atom *HeadAtom;       /* First atom in this structure */
  Atom *TailAtom;       /* Last atom in this structure */
  Bond **ListOfBonds;
  SphericalHarmonicsMatrix *HeadYmatrix; /* First spherical harmonics matrix
                                            based on current coordinates */
  SphericalHarmonicsMatrix *TailYmatrix; /* Last spherical harmonics matrix
                                            based on current coordinates */
  Tensor *
      initialInertiaTensor; /* Inertia tensor with all additional information */
  Tensor *optimizedInertiaTensor; /* Inertia tensor with all additional
                                     information */
  MonteCarloStatistic *monte_carlo_bootstrap;
  MC_Output monte_carlo_output;
  Eigen::MatrixXd optCosineMatrix;
  Eigen::MatrixXd iniCosineMatrix;
  Eigen::MatrixXd rdcMatrix; /* Updated rdc matrix */
  Eigen::MatrixXd
      rdcMatrixScaled; /* Matrix scaled on currend Dmax's from respective
                          coordinates of this structure */
  Eigen::MatrixXd rdcMatrixNorm; /* Matrix normalized on current Dmax's from
                                    respective coordinates of this structure */
  Eigen::Matrix2d rdcVectorSampling;
  Eigen::Vector3d
      centerOfMass; /* Center of Mass. If != (0,0,0) recalculate it! */

 public:
  /* Standard Constructor */
  Structure();
  /* Copy constructor for faster appending of opt structures */
  Structure(Structure *);
  /* Standard destructor */
  ~Structure();

  /* Access Molecule this Structure belongs to */
  Molecule *getParent()
  {
    return parent;
  };
  void setParent(Molecule *nParent)
  {
    parent = nParent;
  };

  /* Access Label of this structure */
  void setLabel(std::string newLabel)
  {
    label = newLabel;
  };
  std::string getLabel()
  {
    return label;
  };

  /* Access index of this structure */
  void setIndex(unsigned int nIndex)
  {
    inputIndex = nIndex;
  };
  unsigned int getIndex()
  {
    return inputIndex;
  };

  /* Access cosine matrix */
  void rescaleCosine(int i, int j)
  {
    iniCosineMatrix.resize(i, j);
    optCosineMatrix.resize(i, j);
  };
  void determineCosineMatrix(StructureOptions);
  Eigen::MatrixXd getCosineMatrix(StructureOptions);

  /* Access rdc matrizes using rdcMatrixOptions (see Molecule.hpp) */
  void setRDCmatrix(rdcMatrixOptions);
  void updateRDCmatrix(BasicInformation &, Flags &, StructureOptions);
  void recalculateRDCs(StructureOptions, Flags &);
  void setUndefined()
  {
    hasUndefinedRDCs = true;
  }
  Eigen::MatrixXd getRDCmatrix(rdcMatrixOptions);
  //      Eigen::MatrixXd getDcalc ( StructureOptions, BasicInformation&,
  //      Eigen::MatrixXd& );
  Eigen::MatrixXd getSaupeTensor();
  Eigen::MatrixXd getSaupeEigenValues();
  Eigen::MatrixXd getSaupeEigenVectors();
  Eigen::MatrixXd getSaupeEulerAngles();

  /* Access structure nodes */
  Structure *getNext() const
  {
    return nextStruc;
  };
  Structure *getPrev() const
  {
    return prevStruc;
  };
  void setNext(Structure *next)
  {
    nextStruc = next;
  }
  void setPrev(Structure *prev)
  {
    prevStruc = prev;
  }

  /* Full initialization of opt Structure (including spherical harmonics and
   * atoms) and rings */
  Structure *
  addOptimizedStructure(Molecule &, BasicInformation &, Structure * = NULL);

  /* Access atoms of this Structure */
  unsigned int getNOA();
  void setTailAtom(Atom *t)
  {
    TailAtom = t;
  };
  Atom *setHeadAtom();
  void setHeadAtom(Atom *a)
  {
    HeadAtom = a;
  };
  Atom *getHeadAtom() const
  {
    return HeadAtom;
  };
  Atom *getTailAtom() const
  {
    return TailAtom;
  };
  Atom *getAtomByIndex(const unsigned int) const;

  /* Access spherical harmonics of this structure */
  SphericalHarmonicsMatrix *addHeadYmatrix();
  void setHeadYmatrix(SphericalHarmonicsMatrix *Y)
  {
    HeadYmatrix = Y;
  }
  SphericalHarmonicsMatrix *getHeadYmatrix() const
  {
    return HeadYmatrix;
  };
  SphericalHarmonicsMatrix *getTailYmatrix() const
  {
    return TailYmatrix;
  };
  bool getIsOptimized() const
  {
    return optimized;
  };
  Eigen::MatrixXd getDmaxMatrix(StructureOptions);

  /* Handling of inertia tensor pas */
  Eigen::Vector3d getCenterOfMass(StructureOptions);
  void Shift2CenterOfMass(StructureOptions);
  Eigen::MatrixXd Rotate2InertiaPAS(StructureOptions);
  Tensor *CalculateInertiaTensor(StructureOptions);
  Tensor *getInertiaTensor(StructureOptions);


  Eigen::MatrixXd getRDCDistances(StructureOptions);
  Eigen::MatrixXd getDistances(StructureOptions);
  unsigned int getNumberOfLongRangeRDCs();
  unsigned int getNOR();

  /* Converts coordinates to xyz structure */
  void retainCoordinates();
  void copyCoordinates(StructureOptions, Structure *, StructureOptions);
  void setUnfixed();

  // Monte Carlo Stuff
  int run_structure_MC(BasicInformation &, Flags &);
  int check_Stop(BasicInformation &, Flags &);
  MC_Output get_monte_carlo_output();
  MC_stop check_Q_diff(BasicInformation &, Flags &);
  MC_stop check_MC_A(BasicInformation &);
  MC_stop check_MC_p(BasicInformation &);

  int initializeRDCindex(Molecule &);
  void NBPP()
  {
    if (!optimized)
      ++NumberOfBonds;
  }
  unsigned int getNumberOfBonds() const
  {
    return NumberOfBonds;
  }
  void setNumberOfBonds(unsigned int b)
  {
    NumberOfBonds = b;
  }
  Bond **initializeListOfBonds();
  Bond **getListOfBonds() const
  {
    return ListOfBonds;
  }
  void count_distance_violations(StructureOptions);
  void
  check_inversion(StructureSimulator &, BasicInformation &, StructureOptions);
  void toxyz(StructureSimulator &,
             RedundantInternals &,
             BasicInformation &,
             Flags &,
             StructureOptions);
  void setPerformedRedundantSteps(const unsigned int steps)
  {
    performedRedundantSteps = steps;
  }
  unsigned int getPerformedRedundantSteps() const
  {
    return performedRedundantSteps;
  }
  void calculateVdW_Potential(StructureSimulator &, StructureOptions);
  void setRDC_rmsd(const double &rmsd)
  {
    rdc_rmsd = rmsd;
  }
  double getRDC_rmsd() const
  {
    return rdc_rmsd;
  }
  void Vector2Structure(StructureOptions);
  void rmsd2Structure(StructureOptions);
  double structure2rmsd(StructureOptions, unsigned int);
  Atom *getMaxStructure2rmsdAtom(unsigned int,
                                 unsigned int &,
                                 Eigen::VectorXd &,
                                 StructureOptions);
  Eigen::VectorXd Qfacs(Eigen::MatrixXd &,
                        Eigen::MatrixXd &,
                        StructureOptions,
                        Eigen::MatrixXd &,
                        Flags &);
  void save_q_factors(BasicInformation &,
                      Flags &,
                      StructureOptions options = StructureOptions::Optimized);
  void deleteBonds();

  void set_rdc_vector_sampling(double, double, StructureOptions);
  void get_rdc_vector_sampling(double &, double &, StructureOptions);
  double get_rdc_vector_sampling_r(StructureOptions);
  double get_rdc_vector_sampling_a(StructureOptions);

  double get_radius_of_Gyration(StructureOptions);

  double calculateAllAtomRMSDofStep();
  double calculateAllAtomRMSD2Input();
  Eigen::VectorXd coordinates2EigenVector(StructureOptions);
  static double calculateAllAtomRMSD(const Eigen::VectorXd,
                                     const Eigen::VectorXd);
};


/* Short function to initialize first spherical harmonics based on input */
int initializeSphericalHarmonics(Molecule &, Structure *);

/* LM function for determination of polar angles from spherical harmonics */
void Yopt(double *, double *, int, int, void *);
void jacYopt(double *, double *, int, int, void *);


#endif

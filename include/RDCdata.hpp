
#ifndef RDCDATA_HPP_
#define RDCDATA_HPP_

#include <eigen3/Eigen/Core>

#ifndef ATOM_HPP_
class Atom;
#endif

#ifndef RDCSET_HPP_
class RDCset;
#endif

#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
enum struct StructureOptions : unsigned int;
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

class RDCdata {
 private:
  unsigned int inputIndex; /* Input index based on first set in input file */
  Atom *rdcAtom1;          /* Link to respective atom */
  Atom *rdcAtom2;          /* Link to respective atom */
  double D;                /* rdc */
  double D_scaled;         /* Scaled rdc. Needed? */
  double D_norm;           /* Normalized rdc. Needed? */
  double deltaD;           /* Error of rdc */
  double weight;           /* hotFCHT weighting of rdc. */
  double effWeight;
  double kappa; /* Respectiv kappa (excluding r^-3 part!!!) */
  Eigen::MatrixXd weighting;
  RDCdata *prevData; /* Link to pervious rdc node */
  RDCdata *nextData; /* Link to next rdc node */
  RDCdata *prevSet;
  RDCdata *nextSet;
  RDCset *parent; /* RDC set which defines this rdc */
  bool defined;
  unsigned int range; /* Defines the amount of bonds between coupling spins */

 public:
  /* Standard constructor */
  RDCdata();

  /* Copy constructor to append rdc data in current set */
  RDCdata(RDCdata *);

  /* Standard destructor */
  ~RDCdata();

  /* Functions to handle nodes */
  RDCdata *appendData();
  RDCdata *getNext() const
  {
    return nextData;
  };
  RDCdata *getPrev() const
  {
    return prevData;
  };
  RDCdata *getNextSetData() const
  {
    return nextSet;
  };
  RDCdata *getPrevSetData() const
  {
    return prevSet;
  };
  void setNext(RDCdata *n)
  {
    nextData = n;
  };
  void setNextSetData(RDCdata *n)
  {
    nextSet    = n;
    n->prevSet = this;
  }
  RDCdata *copyRDC(RDCdata *);

  /* Functions to handle respective rdc sets */
  void setParent(RDCset *p)
  {
    parent = p;
  };
  RDCset *getParent() const
  {
    return parent;
  };

  /* Input index corresponding to first rdc set */
  void setIndex(int ii)
  {
    inputIndex = ii;
  }
  unsigned int getInputIndex() const
  {
    return inputIndex;
  }

  /* Functions to handle linking to rdc atoms */
  void setAtom1(Atom *atom1)
  {
    rdcAtom1 = atom1;
  }
  void setAtom2(Atom *atom2)
  {
    rdcAtom2 = atom2;
  }
  Atom *getAtom1() const
  {
    return rdcAtom1;
  };
  Atom *getAtom2() const
  {
    return rdcAtom2;
  };

  /* Functions to access rdc */
  void setD(double nD)
  {
    D = nD;
  }
  void setD_scaled(double nD)
  {
    D_scaled = nD;
  };
  void setUndefined()
  {
    D       = .0;
    defined = false;
  };
  bool isUndefined() const
  {
    return !defined;
  };
  void setD_norm(double nD)
  {
    D_norm = nD;
  }
  double getD() const
  {
    return D;
  };
  double getD_scaled() const
  {
    return D_scaled;
  };
  double getD_norm() const
  {
    return D_norm;
  };

  /* Functions to access error of rdc */
  void setDeltaD(double ndD)
  {
    deltaD = ndD;
  };
  double getDeltaD() const
  {
    return deltaD;
  };

  /* Functions to access weighting of rdc */
  void setWeight(double nW)
  {
    weight = nW;
  };
  double getWeight() const
  {
    return weight;
  };
  void determineEffectiveWeight(Flags &);
  double getEffectiveWeight() const
  {
    return effWeight;
  }
  Eigen::MatrixXd getVectorWeights();


  /* Functions to access coupling constants used for Dmax */
  void setKappa();
  double getKappa() const
  {
    return kappa;
  };
  double getDistance(StructureOptions) const;
  double getDmax(StructureOptions opt)
  {
    return (getKappa() / pow(getDistance(opt), 3));
  }


  void determineRange(Molecule &);
  unsigned int getRange() const
  {
    return range;
  }
  void setRange(unsigned int r)
  {
    range = r;
  }
};

void initializeRDCmatrix(Molecule &, BasicInformation &, Flags &);

int predictRDCs(Molecule &, BasicInformation &);

double Dcalc(Atom *, Atom *, const Eigen::MatrixXcd &);

double Dcalc(const Eigen::MatrixXcd &,
             const Eigen::MatrixXcd &,
             const double,
             const double);

Eigen::MatrixXd Dcalc(Eigen::MatrixXd &, Eigen::MatrixXd &, Eigen::MatrixXd &);

int linkCascade(RDCset *, RDCdata *, RDCset *, RDCdata *);

#endif

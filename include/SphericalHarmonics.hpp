
#ifndef SPHERICALHARMONICS_HPP_
#define SPHERICALHARMONICS_HPP_

#include <eigen3/Eigen/Core>

#ifndef ATOM_HPP_
class Atom;
#endif

#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef STRUCTURE_HPP_
class Structure;
#endif

#ifndef RDCDATA_HPP_
class RDCdata;
#endif

#ifndef SPHERICALHARMONICSMATRIX_HPP_
class SphericalHarmonicsMatrix;
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
enum struct StructureOptions : unsigned int;
#endif

class SphericalHarmonics {
 private:
  double theta;          /* Polar angle theta pre optimization*/
  double phi;            /* Polar angle phi pre optimization */
  double optTheta;       /* Polar angle theta after optimization */
  double optPhi;         /* Polar angle phi after optimization */
  double orderParameter; /* Order parameter S^2 */
  double S_axial;        /* Axial symmetrical order parameter */
  double eta;
  double aniso_theta;
  double *LM_info;
  SphericalHarmonics *prevEl; /* Next harmonics node */
  SphericalHarmonics *nextEl; /* Previous harmonics node */
  SphericalHarmonicsMatrix
      *parent; /* Harmonics matrix which defines this harmonic */
  Atom *atom1; /* rdc atom which defines this harmonic (see RDCdata.cpp) */
  Atom *atom2; /* rdc atom which defines this harmonic (see RDCdata.cpp) */
  unsigned inputIndex;      /* Indes of rdc which defines this harmonic (see
                               RDCdata.cpp) -> 1 based */
  Eigen::MatrixXd Brow;     /* Row of B (normalized cosine) matrix which defines
                               this harmonic */
  Eigen::MatrixXd optBrow;  /* Optimized version of Brow */
  Eigen::MatrixXcd optYrow; /* Vector containing all optimized harmonics */
  Eigen::MatrixXcd Yrow;    /* Vector containing all harmonics */
  RDCdata *rdc;             /* rdc which defines this harmonic */
  unsigned int range; /* Defines the amount of bonds between coupling spins */

 public:
  /* Standard constructor using polar angles as parameter */
  SphericalHarmonics(const double, const double);
  /* Copy constructor for initialization of optimized harmonics */
  SphericalHarmonics(SphericalHarmonics *,
                     SphericalHarmonicsMatrix *,
                     Structure *);

  /* Standard destructor */
  ~SphericalHarmonics();

  /* Functions to handle harmonics nodes */
  SphericalHarmonics *appendHarmonics(Molecule &);
  SphericalHarmonics *getNext() const
  {
    return nextEl;
  };
  SphericalHarmonics *getPrev() const
  {
    return prevEl;
  };
  SphericalHarmonics *initializeOptimizedHarmonics(SphericalHarmonicsMatrix *,
                                                   SphericalHarmonicsMatrix *);

  /* Handle (optimized) polar angles */
  double getTheta(StructureOptions) const;       // { return theta; };
  double getPhi(StructureOptions) const;         // { return phi; };
  void setTheta(const double, StructureOptions); // { theta = t; };
  void setPhi(const double, StructureOptions);   // { phi = p; };
  void setOptimizedAngles(double, double);
  void recalculateAngles(StructureOptions);

  void setS2rdc(double S2)
  {
    orderParameter = S2;
  }
  void setSaxial(double Sax)
  {
    S_axial = Sax;
  }
  void setetardc(double e)
  {
    eta = e;
  }
  void setAnisoTheta(double t)
  {
    aniso_theta = t;
  }
  double getS2rdc() const
  {
    return orderParameter;
  }
  double getSaxial() const
  {
    return S_axial;
  }
  double getetardc() const
  {
    return eta;
  }
  double getAnisoTheta() const
  {
    return aniso_theta;
  }

  /* Hande parent matrix */
  void setParent(SphericalHarmonicsMatrix *p)
  {
    parent = p;
  };

  /* Handle rdc atoms */
  void setAtoms(Atom *a1, Atom *a2)
  {
    atom1 = a1;
    atom2 = a2;
  };
  Atom *getAtom1() const
  {
    return atom1;
  };
  Atom *getAtom2() const
  {
    return atom2;
  };
  Atom *getPartner(Atom *a)
  {
    return (a == atom1 ? atom2 : atom1);
  }
  /* Handle indizes */
  void setInputIndex(const unsigned int i)
  {
    inputIndex = i;
  };
  unsigned int getInputIndex() const
  {
    return inputIndex;
  };

  /* Handle cosine/Y matrizes */
  Eigen::MatrixXd getBrow() const
  {
    return Brow;
  };

  Eigen::MatrixXd getOptimizedBrow() const
  {
    return Brow;
  };
  Eigen::MatrixXcd getYrow(StructureOptions); // const { return Yrow; };
  void initializeYrow();


  unsigned int getRange() const
  {
    return range;
  }
  void setRange(unsigned int r)
  {
    range = r;
  }

  RDCdata *getRDC() const
  {
    return rdc;
  }
  double getDmax(StructureOptions);
  double getD() const;
  double getDeltaD() const;
  void setRDC(RDCdata *r)
  {
    rdc = r;
  }
  double get_sigma_square();

  double *getLMinfo()
  {
    return LM_info;
  }
};
#endif

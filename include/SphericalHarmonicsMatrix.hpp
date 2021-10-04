
#ifndef SPHERICALHARMONICSMATRIX_HPP_
#define SPHERICALHARMONICSMATRIX_HPP_

#include <eigen3/Eigen/Core>

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

#ifndef SPHERICALHARMONICS_HPP_
class SphericalHarmonics;
#endif

#ifndef STRUCTURE_HPP_
class Structure;
#endif

#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
enum struct StructureOptions : unsigned int;
#endif

class SphericalHarmonicsMatrix {
 private:
  Eigen::MatrixXd SaupeTensor;
  Eigen::MatrixXd
      SaupeEulerAngles; /* Vector representation of positive euler angles */
  Eigen::MatrixXd SaupeEigenValues; /* Matrix representation of eigen values */
  Eigen::MatrixXd SaupeEigenVectors;
  Eigen::MatrixXd Bmatrix; /* Latest cosine matrix */
  Eigen::MatrixXcd
      Fmatrix; /* Latest F matrix (containing R, alpha, beta & gamma ) */
  Eigen::MatrixXcd
      Ymatrix; /* Latest matrix representation of all spherical harmonics */
  Eigen::MatrixXcd Yref; /* Refined structure model as Ymatrix */
  SphericalHarmonics
      *HeadEl; /* First spherical harmonic (via input index of rdcs) */
  SphericalHarmonics
      *TailEl;       /* Last spherical harmonic (via input index of rdcs) */
  Structure *parent; /* Structure of this harmonic matrix */
  SphericalHarmonicsMatrix *nextEl; /* Node for next harmonics matrix */
  SphericalHarmonicsMatrix *prevEl; /* Node for previous harmonics matrix */
  double Soverall;

 public:
  /* Standard constructor */
  SphericalHarmonicsMatrix();
  /* Copy constructor for optimized harmonics matrix */
  SphericalHarmonicsMatrix(SphericalHarmonicsMatrix *, Structure *);

  /* Standard destructor */
  ~SphericalHarmonicsMatrix();

  /* Handle harmonics  */
  SphericalHarmonics *addHeadElement();
  SphericalHarmonics *getHeadHarmonic() const
  {
    return HeadEl;
  };
  void setTailEl(SphericalHarmonics *t)
  {
    TailEl = t;
  };

  /* Handle structure defining this harmonics matrix */
  void setParent(Structure *p)
  {
    parent = p;
  };
  Structure *getParent() const
  {
    return parent;
  };

  /* Handle harmonics matrix nodes */
  void setNext(SphericalHarmonicsMatrix *n)
  {
    nextEl = n;
  };
  SphericalHarmonicsMatrix *getNext() const
  {
    return nextEl;
  };
  SphericalHarmonicsMatrix *getPrev() const
  {
    return prevEl;
  };


  /* Functions to handle structure/orientation describing matrizes */
  Eigen::MatrixXd determineBmatrix(Molecule &, Structure *, StructureOptions);
  Eigen::MatrixXd getBmatrix() const
  {
    return Bmatrix;
  };

  Eigen::MatrixXcd getYref() const
  {
    return Yref;
  };
  void setYref(Eigen::MatrixXcd Y)
  {
    Yref = Y;
  };
  void setBmatrix(Eigen::MatrixXd &B)
  {
    Bmatrix = B;
  };

  /* Eigen value decomposition of alignment matrix */
  void determineEuler(Molecule &,
                      Structure *,
                      BasicInformation &,
                      Flags &,
                      StructureOptions);
  Eigen::MatrixXd getSaupeTensor() const
  {
    return SaupeTensor;
  };
  Eigen::MatrixXd getSaupeEulerAngles() const
  {
    return SaupeEulerAngles;
  };
  Eigen::MatrixXd getSaupeEigenValues() const
  {
    return SaupeEigenValues;
  };
  Eigen::MatrixXd getSaupeEigenVectors() const
  {
    return SaupeEigenVectors;
  };

  /* Functions that handle matrix representation of F (orientational term) */
  Eigen::MatrixXcd determineFmatrix();
  Eigen::MatrixXcd getFmatrix() const
  {
    return Fmatrix;
  };

  void setSoverall(double S)
  {
    Soverall = S;
  }
  double getSoverall() const
  {
    return Soverall;
  }
};

#endif

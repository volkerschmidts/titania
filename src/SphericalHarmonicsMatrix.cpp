
#include <Atom.hpp>
#include <Declarations.hpp>
#include <LinAlg.hpp>
#include <Molecule.hpp>
#include <RDCdata.hpp>
#include <RDCset.hpp>
#include <SphericalHarmonics.hpp>
#include <SphericalHarmonicsMatrix.hpp>
#include <Structure.hpp>
#include <iostream>

SphericalHarmonicsMatrix::SphericalHarmonicsMatrix()
{
  SaupeEulerAngles = Eigen::MatrixXd::Zero(1, 1);
  SaupeEigenValues = Eigen::MatrixXd::Zero(1, 1);
  Bmatrix          = Eigen::MatrixXd::Zero(1, 1);
  Fmatrix          = Eigen::MatrixXd::Zero(1, 1);
  Ymatrix          = Eigen::MatrixXcd::Zero(1, 1);
  HeadEl = TailEl = NULL;
  nextEl = prevEl = NULL;
  parent          = NULL;
  Soverall        = .0;
}

SphericalHarmonicsMatrix::SphericalHarmonicsMatrix(
    SphericalHarmonicsMatrix *SHM,
    Structure *currStruc)
{
  SaupeEulerAngles = Eigen::MatrixXd::Zero(1, 1);
  SaupeEigenValues = Eigen::MatrixXd::Zero(1, 1);
  Bmatrix          = Eigen::MatrixXd::Zero(1, 1);
  Fmatrix          = Eigen::MatrixXd::Zero(1, 1);
  Ymatrix          = Eigen::MatrixXcd::Zero(1, 1);
  parent           = currStruc;
  HeadEl = new SphericalHarmonics(SHM->getHeadHarmonic(), this, parent);
  TailEl = HeadEl->initializeOptimizedHarmonics(SHM, this);
  nextEl = prevEl = NULL;
  Soverall        = .0;
}

SphericalHarmonicsMatrix::~SphericalHarmonicsMatrix()
{
  while (HeadEl->getNext())
  {
    HeadEl = HeadEl->getNext();
    delete HeadEl->getPrev();
  }

  delete HeadEl;
  TailEl = NULL;

  parent = NULL;
  prevEl = nextEl = NULL;
}

SphericalHarmonics *
SphericalHarmonicsMatrix::addHeadElement()
{
  RDCdata *CurrRDC = parent->getParent()->getHeadSet()->getHeadData();

  Eigen::Vector3d V =
      CurrRDC->getAtom1()->Coordinates2Eigen(StructureOptions::Initial) -
      CurrRDC->getAtom2()->Coordinates2Eigen(StructureOptions::Initial);
  V = V / V.norm();

  HeadEl =
      new SphericalHarmonics(acos(V(STANDARD_AXIS_Z_)),
                             atan2(V(STANDARD_AXIS_Y_), V(STANDARD_AXIS_X_)));
  HeadEl->setRDC(CurrRDC);
  HeadEl->setRange(CurrRDC->getRange());
  HeadEl->setAtoms(CurrRDC->getAtom1(), CurrRDC->getAtom2());
  HeadEl->getAtom1()->addBondharmonic(HeadEl);
  HeadEl->getAtom2()->addBondharmonic(HeadEl);
  HeadEl->setParent(this);
  HeadEl->setAtoms(CurrRDC->getAtom1(), CurrRDC->getAtom2());
  HeadEl->setInputIndex(1);
  CurrRDC = NULL;
  return HeadEl;
}

Eigen::MatrixXd
SphericalHarmonicsMatrix::determineBmatrix(Molecule &CurrMol,
                                           Structure *CurrStruc,
                                           StructureOptions opt)
{
  unsigned int NOR = CurrMol.getNOR();
  SphericalHarmonics *CurrY;
  unsigned int i;

  if (Bmatrix.rows() != NOR || Bmatrix.cols() != COSINE_ELEMENTS_)
    Bmatrix.resize(NOR, COSINE_ELEMENTS_);
  if (Ymatrix.rows() != NOR || Ymatrix.cols() != HARMONIC_ELEMENTS_)
    Ymatrix.resize(NOR, HARMONIC_ELEMENTS_);

  for (i = 0, CurrY = CurrStruc->getHeadYmatrix()->getHeadHarmonic(); CurrY;
       ++i, CurrY   = CurrY->getNext())
  {
    Bmatrix.row(i) =
        Sphere2B(CurrY->getTheta(opt), CurrY->getPhi(opt), .0, .0, .0);
    Ymatrix.row(i) = CurrY->getYrow(opt).row(0);
  }
  return Bmatrix;
}

void
SphericalHarmonicsMatrix::determineEuler(Molecule &CurrMol,
                                         Structure *CurrStruc,
                                         BasicInformation &baseInformation,
                                         Flags &flags,
                                         StructureOptions options)
{
  unsigned int NOR, NOS;
  NOR = CurrMol.getNOR();
  NOS = CurrMol.getNORsets();

  CurrStruc->updateRDCmatrix(baseInformation, flags, options);

  Eigen::MatrixXd D = CurrStruc->getRDCmatrix(rdcMatrixOptions::Scaled);
  SaupeEulerAngles  = Eigen::MatrixXd::Zero(NUM_EULER_ANGLES_, NOS);
  SaupeEigenValues  = Eigen::MatrixXd::Zero(NUM_EIGENVALUES_, NOS);
  SaupeEigenVectors = Eigen::MatrixXd::Zero(NUM_EIGENVECTORS_, NOS);
  SaupeTensor       = Eigen::MatrixXd::Zero(SAUPE_ELEMENTS_, NOS);

  Eigen::MatrixXd weights;
  Eigen::MatrixXd B = Bmatrix;

  if (flags.calculateFullMatrix)
    weights = Eigen::MatrixXd::Identity(NOR, NOR);
  else
    weights = CurrMol.getWmatrix();

  SaupeEigenSystems(B, weights, D, SaupeTensor, SaupeEigenValues,
                    SaupeEulerAngles, SaupeEigenVectors,
                    flags.calculateFullMatrix);
}

Eigen::MatrixXcd
SphericalHarmonicsMatrix::determineFmatrix()
{
  Fmatrix = S2F(SaupeEigenValues, SaupeEulerAngles);
  return Fmatrix;
}

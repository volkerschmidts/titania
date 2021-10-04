
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
#include <phobos.h>

SphericalHarmonics::SphericalHarmonics(const double theta = .0,
                                       const double phi   = .0)
{
  /* Standard values */
  prevEl = nextEl = NULL;
  parent          = NULL;
  atom1 = atom2 = NULL;
  inputIndex    = 0;
  rdc           = NULL;
  optPhi = optTheta = .0;
  orderParameter = S_axial = eta = aniso_theta = .0;
  LM_info = (double *) malloc(PHOBOS_INFOSIZE * sizeof(double));
  range   = 1;
  /* Parameters from arguments */
  this->theta = theta;
  this->phi   = phi;
  Brow        = Sphere2B(theta, phi, .0, .0, .0);
  Yrow        = Eigen::MatrixXcd::Zero(1, HARMONIC_ELEMENTS_);
  optYrow     = Eigen::MatrixXcd::Zero(1, HARMONIC_ELEMENTS_);
  this->initializeYrow();
}

SphericalHarmonics::SphericalHarmonics(SphericalHarmonics *SH,
                                       SphericalHarmonicsMatrix *p,
                                       Structure *CurrStruc)
{
  theta  = SH->optTheta;
  phi    = SH->optPhi;
  prevEl = nextEl = NULL;
  parent          = p;
  atom1           = CurrStruc->getAtomByIndex(SH->getAtom1()->getIndex());
  atom2           = CurrStruc->getAtomByIndex(SH->getAtom2()->getIndex());
  range = SH->range, rdc = SH->rdc;
  optPhi = optTheta = .0;
  orderParameter = S_axial = eta = aniso_theta = .0;
  LM_info = (double *) malloc(PHOBOS_INFOSIZE * sizeof(double));
  atom1->addBondharmonic(this);
  atom2->addBondharmonic(this);
  inputIndex = SH->inputIndex;
  Brow       = Sphere2B(theta, phi, .0, .0, .0);
  Yrow       = SH->getYrow(StructureOptions::Optimized);
  optYrow    = Eigen::MatrixXcd::Zero(1, HARMONIC_ELEMENTS_);
  this->initializeYrow();
}

SphericalHarmonics::~SphericalHarmonics()
{
  /* Standard values */
  rdc    = NULL;
  prevEl = nextEl = NULL;
  parent          = NULL;
  atom1 = atom2 = NULL;
  free(LM_info);
}

double
SphericalHarmonics::getTheta(enum StructureOptions options) const
{
  if (options & StructureOptions::Initial)
  {
    return theta;
  }
  else if (options & StructureOptions::Optimized)
  {
    return optTheta;
  }
  else
  {
    std::cerr << "ERROR:\tUnknown type of polar angle theta was "
                 "requested...\n\t\tReturning initial theta...\n";
    return theta;
  }
}

double
SphericalHarmonics::getPhi(enum StructureOptions options) const
{
  if (options & StructureOptions::Initial)
  {
    return phi;
  }
  else if (options & StructureOptions::Optimized)
  {
    return optPhi;
  }
  else
  {
    std::cerr << "ERROR:\tUnknown type of polar angle phi was "
                 "requested...\n\t\tReturning initial phi...\n";
    return phi;
  }
}

void
SphericalHarmonics::setTheta(const double t, enum StructureOptions options)
{
  if (options & StructureOptions::Initial)
  {
    theta = t;
  }
  else if (options & StructureOptions::Optimized)
  {
    optTheta = t;
  }
  else
  {
    std::cerr << "ERROR:\tUnknown type of polar angle theta was "
                 "requested...\n\t\tChanging optimized theta...\n";
    optTheta = t;
  }
}

void
SphericalHarmonics::setPhi(const double p, enum StructureOptions options)
{
  if (options & StructureOptions::Initial)
  {
    phi = p;
  }
  else if (options & StructureOptions::Optimized)
  {
    optPhi = p;
  }
  else
  {
    std::cerr << "ERROR:\tUnknown type of polar angle phi was "
                 "requested...\n\t\tChanging optimized phi...\n";
    optPhi = p;
  }
}

void
SphericalHarmonics::setOptimizedAngles(double optT, double optP)
{
  double iniT, iniP, invT, invP, optA, invA;

  iniT = theta;
  iniP = phi;
  invT = PI_ - optT;
  if (optP > .0)
    invP = optP - PI_;
  else
    invP = PI_ + optP;
  Eigen::Vector3d ini = Polar2Eigen(iniT, iniP);
  Eigen::Vector3d inv = Polar2Eigen(invT, invP);
  Eigen::Vector3d opt = Polar2Eigen(optT, optP);

  optA = ini.transpose() * opt;
  invA = ini.transpose() * inv;

  if (optA > .0)
  {
    optTheta = optT;
    optPhi   = optP;
    return;
  }
  else if (invA >= .0)
  {
    optTheta = invT;
    optPhi   = invP;
    return;
  }
  else
    std::cerr << "WARING:\tSomething went pretty wrong in with polar angles of "
              << atom1->getIdentifier() << " " << atom2->getIdentifier()
              << "...\n\t\tSetting the polar angles to initial state...\n";
  optTheta = theta;
  optPhi   = phi;
}

void
SphericalHarmonics::recalculateAngles(StructureOptions opt)
{
  Eigen::Vector3d V =
      atom1->Coordinates2Eigen(opt) - atom2->Coordinates2Eigen(opt);
  double tmp_T, tmp_P;

  Eigen2Polar(V, tmp_T, tmp_P);
  setTheta(tmp_T, opt);
  setPhi(tmp_P, opt);
}

SphericalHarmonics *
SphericalHarmonics::initializeOptimizedHarmonics(
    SphericalHarmonicsMatrix *preMatrix,
    SphericalHarmonicsMatrix *p)
{
  Structure *CurrStruc             = p->getParent();
  SphericalHarmonics *CurrHarmonic = preMatrix->getHeadHarmonic();

  SphericalHarmonics *newHarmonic;
  CurrHarmonic        = CurrHarmonic->getNext();
  newHarmonic         = new SphericalHarmonics(CurrHarmonic, p, p->getParent());
  this->nextEl        = newHarmonic;
  newHarmonic->prevEl = this;
  while (CurrHarmonic->getNext())
  {
    CurrHarmonic = CurrHarmonic->nextEl;
    newHarmonic->nextEl =
        new SphericalHarmonics(CurrHarmonic, p, p->getParent());
    newHarmonic->nextEl->prevEl = newHarmonic;
    newHarmonic                 = newHarmonic->nextEl;
  }

  CurrStruc    = NULL;
  CurrHarmonic = NULL;
  delete CurrStruc;
  delete CurrHarmonic;
  return newHarmonic;
}

SphericalHarmonics *
SphericalHarmonics::appendHarmonics(Molecule &CurrMol)
{
  RDCdata *CurrRDC = CurrMol.getHeadSet()->getHeadData();
  unsigned int i;
  for (i = 1; i <= inputIndex; ++i)
    CurrRDC = CurrRDC->getNext();

  Eigen::Vector3d V =
      CurrRDC->getAtom1()->Coordinates2Eigen(StructureOptions::Initial) -
      CurrRDC->getAtom2()->Coordinates2Eigen(StructureOptions::Initial);
  V = V / V.norm();

  SphericalHarmonics *tmp =
      new SphericalHarmonics(acos(V(2)), atan2(V(1), V(0)));
  tmp->rdc   = CurrRDC;
  tmp->atom1 = CurrRDC->getAtom1();
  tmp->atom2 = CurrRDC->getAtom2();
  tmp->getAtom1()->addBondharmonic(tmp);
  tmp->getAtom2()->addBondharmonic(tmp);
  tmp->prevEl     = this;
  tmp->parent     = parent;
  tmp->inputIndex = inputIndex + 1;
  tmp->range      = CurrRDC->getRange();
  nextEl          = tmp;
  parent->setTailEl(tmp);
  CurrRDC = NULL;

  return tmp;
}

void
SphericalHarmonics::initializeYrow()
{
  unsigned int i;
  for (i = 0; i < HARMONIC_ELEMENTS_; ++i)
    Yrow(0, i) = Ylm(2, i - 2, theta, phi);
}

double
SphericalHarmonics::getDmax(StructureOptions opt)
{
  return (rdc->getKappa() / pow((atom1->getDistance(atom2, opt) * 1e-10), 3.0));
}

double
SphericalHarmonics::getD() const
{
  return rdc->getD();
}
double
SphericalHarmonics::getDeltaD() const
{
  return rdc->getDeltaD();
}

double
SphericalHarmonics::get_sigma_square()
{
  double sigmasq = .0;
  RDCdata *RDC;
  for (RDC = rdc; RDC; RDC = RDC->getNextSetData())
  {
    sigmasq += (RDC->getDeltaD() * RDC->getDeltaD());
  }
  return sigmasq;
}

Eigen::MatrixXcd
SphericalHarmonics::getYrow(StructureOptions opt)
{
  for (int i = 0; i < HARMONIC_ELEMENTS_; ++i)
  {
    optYrow(0, i) = Ylm(2, i - 2, getTheta(opt), getPhi(opt));
  }
  return optYrow;
}

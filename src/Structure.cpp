
#include <Atom.hpp>
#include <Bond.hpp>
#include <Declarations.hpp>
#include <Eckart.hpp>
#include <Helper.hpp>
#include <LinAlg.hpp>
#include <Molecule.hpp>
#include <MonteCarloStatistic.hpp>
#include <Output.hpp>
#include <RDCdata.hpp>
#include <RDCset.hpp>
#include <RedundantInternals.hpp>
#include <SphericalHarmonics.hpp>
#include <SphericalHarmonicsMatrix.hpp>
#include <Structure.hpp>
#include <StructureSimulator.hpp>
#include <complex>
#include <iomanip>
#include <iostream>


Structure::Structure(Structure *s)
{
  label                   = "";
  inputIndex              = 0;
  performedRedundantSteps = 0;
  parent                  = s->parent;
  iniCosineMatrix         = Eigen::MatrixXd::Zero(1, 1);
  optCosineMatrix         = Eigen::MatrixXd::Zero(1, 1);
  prevStruc               = s;
  nextStruc               = NULL;
  TailAtom                = NULL;
  HeadAtom                = s->HeadAtom;
  optimized               = false;
  hasUndefinedRDCs        = false;
  radiusOfGyration        = .0;
  rdc_rmsd                = .0;
  initialInertiaTensor    = NULL;
  optimizedInertiaTensor  = NULL;
  NumberOfBonds           = 0;
  ListOfBonds             = NULL;
  HeadYmatrix             = NULL;
  TailYmatrix             = NULL;
  monte_carlo_bootstrap   = NULL;
}

Structure::Structure()
{
  label                   = "";
  inputIndex              = 0;
  performedRedundantSteps = 0;
  prevStruc               = NULL;
  nextStruc               = NULL;
  parent                  = NULL;
  HeadAtom                = NULL;
  TailAtom                = NULL;
  iniCosineMatrix         = Eigen::MatrixXd::Zero(1, 1);
  optCosineMatrix         = Eigen::MatrixXd::Zero(1, 1);
  optimized               = false;
  hasUndefinedRDCs        = false;
  radiusOfGyration        = .0;
  rdc_rmsd                = .0;
  initialInertiaTensor    = NULL;
  optimizedInertiaTensor  = NULL;
  NumberOfBonds           = 0;
  ListOfBonds             = NULL;
  HeadYmatrix             = NULL;
  TailYmatrix             = NULL;
  monte_carlo_bootstrap   = NULL;
}

Structure::~Structure()
{
  parent    = NULL;
  prevStruc = NULL;
  nextStruc = NULL;
  if (NumberOfBonds && ListOfBonds)
  {
    for (unsigned int b = 0; b < NumberOfBonds; ++b)
    {
      delete ListOfBonds[b];
    }
  }
  free(ListOfBonds);

  while (HeadAtom->getNext())
  {
    HeadAtom = HeadAtom->getNext();
    delete HeadAtom->getPrev();
  }
  TailAtom = NULL;
  delete HeadAtom;

  while (HeadYmatrix && HeadYmatrix->getNext())
  {
    HeadYmatrix = HeadYmatrix->getNext();
    delete HeadYmatrix->getPrev();
  }
  TailYmatrix = NULL;
  delete HeadYmatrix;

  delete initialInertiaTensor;
  delete optimizedInertiaTensor;
  delete monte_carlo_bootstrap;
}

Structure *
Structure::addOptimizedStructure(Molecule &CurrMol,
                                 BasicInformation &baseInformation,
                                 Structure *copy)
{
  static int optstep = 2;
  unsigned int b, bi;
  if (copy == NULL)
  {
    copy = this;
  }

  Structure *newS                   = new Structure();
  CurrMol.getTailStruc()->nextStruc = newS;
  CurrMol.raiseNumberOfStrucs();
  newS->prevStruc        = parent->getTailStruc();
  newS->optimized        = true;
  newS->hasUndefinedRDCs = copy->hasUndefinedRDCs;
  parent->setTailStruc(newS);
  newS->label = CurrMol.getHeadStruc()->getLabel() +
                baseInformation.structureLabelExtension +
                std::to_string(optstep++);
  newS->inputIndex               = copy->inputIndex + 1;
  newS->NumberOfBonds            = copy->NumberOfBonds;
  newS->rdcMatrix                = copy->rdcMatrix;
  newS->rdcMatrixScaled          = copy->rdcMatrixScaled;
  newS->rdcMatrixNorm            = copy->rdcMatrixNorm;
  newS->rdcVectorSampling.row(0) = copy->rdcVectorSampling.row(1);
  newS->HeadAtom                 = new Atom();
  newS->TailAtom = newS->HeadAtom->initializeOptimizedAtoms(copy, newS);

  Atom *CurrAtom = newS->HeadAtom;
  newS->initializeListOfBonds();
  newS->parent = copy->parent;
  while (CurrAtom)
  {
    for (b = 0; b < CurrAtom->getNumberOfBonds(); ++b)
    {
      for (bi = 0; bi < NumberOfBonds; ++bi)
      {
        if (CurrAtom->getBond(b) == newS->ListOfBonds[bi])
          break;
        else if (newS->ListOfBonds[bi] == NULL)
        {
          newS->ListOfBonds[bi] = CurrAtom->getBond(b);
          CurrAtom->getBond(b)->copyOrder();
          break;
        }
      }
    }
    CurrAtom = CurrAtom->getNext();
  }
  // For some reason the bond harmonics pointer in Atom does not stay
  // NULL for HeadAtom after initialization. By this the harmonics
  // cannot be linked correctly! To be sure all bond harmonics pointer
  // are cleared again.

  CurrAtom = newS->HeadAtom;
  while (CurrAtom)
  {
    CurrAtom->clearHarmonics();
    CurrAtom = CurrAtom->getNext();
  }

  newS->iniCosineMatrix = copy->optCosineMatrix;
  newS->optCosineMatrix = Eigen::MatrixXd::Zero(1, 1);
  newS->HeadYmatrix     = new SphericalHarmonicsMatrix(HeadYmatrix, newS);
  newS->TailYmatrix     = NULL;

  CurrAtom = NULL;
  delete CurrAtom;

  return newS;
}

unsigned int
Structure::getNOA()
{
  if (parent)
    return parent->getNOA();
  else
    return 0;
}

Atom *
Structure::setHeadAtom()
{
  Atom *tmp = new Atom();
  HeadAtom  = tmp;
  tmp->setParent(this);
  TailAtom = tmp;
  return tmp;
}

void
Structure::setRDCmatrix(enum rdcMatrixOptions rdcFlag)
{
  switch (static_cast<unsigned int>(rdcFlag))
  {
    case (1):
      rdcMatrix = parent->getRDCmatrix();
      break;
    default:
      std::cerr << "ERROR:\tTried to save a rdc matrix with wrong options to "
                   "Structure...\n";
      break;
  }
}

void
Structure::recalculateRDCs(StructureOptions options, Flags &flags)
{
  if (!hasUndefinedRDCs || !HeadYmatrix)
    return;
  RDCset *CurrSet = parent->getHeadSet();
  RDCdata *CurrData;

  unsigned int NOS         = parent->getNORsets();
  unsigned int NOR         = parent->getNOR();
  Eigen::MatrixXcd Yref    = HeadYmatrix->getYref();
  Eigen::MatrixXcd Fmatrix = HeadYmatrix->getFmatrix();
  Eigen::MatrixXd Saupe;
  Eigen::MatrixXd Evals = Eigen::MatrixXd::Zero(NUM_EIGENVALUES_, NOS);
  Eigen::MatrixXd Eang  = Eigen::MatrixXd::Zero(NUM_EULER_ANGLES_, NOS);
  Eigen::MatrixXd Evec  = Eigen::MatrixXd::Zero(NUM_EIGENVECTORS_, NOS);
  Eigen::MatrixXd back_calc;

  determineCosineMatrix(options);
  Eigen::MatrixXd C = getCosineMatrix(options);
  Eigen::MatrixXd w;
  if (flags.calculateFullMatrix)
    w = Eigen::MatrixXd::Identity(NOR, NOR);
  else
    w = parent->getWmatrix();

  Eigen::MatrixXd Dmax = getDmaxMatrix(options);
  SaupeEigenSystems(C, w, rdcMatrixScaled, Saupe, Evals, Eang, Evec,
                    flags.calculateFullMatrix);
  back_calc = Dcalc(Dmax, C, Saupe);

  *std::cin.tie() << "Information:\tApproximating RDCs for new structure:"
                  << std::endl;
  while (CurrSet)
  {
    CurrData = CurrSet->getHeadData();
    while (CurrData)
    {
      if (CurrData->isUndefined())
      {
        std::string rdc_ident = (CurrData->getAtom1()->getIdentifier() + "-" +
                                 CurrData->getAtom2()->getIdentifier());

        *std::cin.tie() << std::string(15, ' ') << std::setw(10) << std::left
                        << rdc_ident;
        *std::cin.tie() << "[" << std::setw(15) << CurrSet->getLabel() << "]: ";
        rdcMatrix(CurrData->getInputIndex(), CurrSet->getIndex() - 1) =
            back_calc(CurrData->getInputIndex(), CurrSet->getIndex() - 1);
        *std::cin.tie() << std::right << std::fixed << std::setw(10)
                        << std::setprecision(4)
                        << rdcMatrix(CurrData->getInputIndex(),
                                     CurrSet->getIndex() - 1)
                        << std::endl;
      }
      CurrData = CurrData->getNext();
    }
    CurrSet = CurrSet->getNext();
  }
}

void
Structure::updateRDCmatrix(BasicInformation &baseInformation,
                           Flags &flags,
                           StructureOptions options)
{
  const unsigned int NOR = rdcMatrix.rows(), NOS = rdcMatrix.cols();
  unsigned int r, s;
  if (NOR <= 1 || NOS <= 1)
  {
    std::cerr << "ERROR:\tFunction Structure::updateRDCmatrix was called "
                 "without a properly defined RDC matrix...\n"
              << "\tKilling request until RDC matrix was set...\n";
    return;
  }
  double scalingFac;

  RDCdata *CurrRDC = parent->getHeadSet()->getHeadData();
  if (flags.recalculateRDCs)
    this->recalculateRDCs(options, flags);
  rdcMatrixScaled = rdcMatrixNorm = rdcMatrix;

  for (r = 0; r < NOR; ++r)
  {
    scalingFac = pow(CurrRDC->getDistance(options), 3.0) / CurrRDC->getKappa();
    for (s = 0; s < NOS; ++s)
    {
      rdcMatrixScaled(r, s) = rdcMatrix(r, s) * scalingFac;
      rdcMatrixNorm(r, s) = rdcMatrixScaled(r, s) * baseInformation.normFactRDC;
    }
    CurrRDC = CurrRDC->getNext();
  }
  CurrRDC = NULL;
}

Eigen::MatrixXd
Structure::getRDCmatrix(rdcMatrixOptions rdcFlag)
{
  switch (rdcFlag)
  {
    case (rdcMatrixOptions::Unscaled):
      return rdcMatrix; /* Unscaled */
    case (rdcMatrixOptions::Scaled):
      return rdcMatrixScaled; /* Scaled   */
    case (rdcMatrixOptions::Norm):
      return rdcMatrixNorm; /* Norm     */
    default:
      std::cerr << "ERROR:\tTried to load a rdc matrix with wrong options from "
                   "Structure...\n";
      return Eigen::MatrixXd::Zero(1, 1);
  }
}

Eigen::MatrixXd
Structure::getSaupeTensor()
{
  return HeadYmatrix->getSaupeTensor();
}
Eigen::MatrixXd
Structure::getSaupeEigenValues()
{
  return HeadYmatrix->getSaupeEigenValues();
}
Eigen::MatrixXd
Structure::getSaupeEigenVectors()
{
  return HeadYmatrix->getSaupeEigenVectors();
}
Eigen::MatrixXd
Structure::getSaupeEulerAngles()
{
  return HeadYmatrix->getSaupeEulerAngles();
}

SphericalHarmonicsMatrix *
Structure::addHeadYmatrix()
{
  HeadYmatrix = new SphericalHarmonicsMatrix();
  HeadYmatrix->setParent(this);
  TailYmatrix = HeadYmatrix;
  return HeadYmatrix;
}

Eigen::MatrixXd
Structure::getDmaxMatrix(StructureOptions opt)
{
  Eigen::MatrixXd Dmax =
      Eigen::MatrixXd::Zero(parent->getNOR(), parent->getNOR());
  SphericalHarmonics *Y;
  unsigned int i;
  for (i = 0, Y = HeadYmatrix->getHeadHarmonic(); Y; Y = Y->getNext(), ++i)
  {
    Dmax(i, i) =
        Y->getRDC()->getKappa() / pow((Y->getRDC()->getDistance(opt)), 3.0);
  }
  return Dmax;
}

Atom *
Structure::getAtomByIndex(const unsigned int index) const
{
  if (index == 1)
    return HeadAtom;
  Atom *tmp = HeadAtom;
  while (tmp)
  {
    if (tmp->getIndex() != index)
    {
      tmp = tmp->getNext();
      continue;
    }
    else
      break;
  }
  return tmp;
}

Eigen::Vector3d
Structure::getCenterOfMass(enum StructureOptions options)
{
  centerOfMass   = Eigen::Vector3d::Zero(3);
  Atom *CurrAtom = HeadAtom;
  while (CurrAtom)
  {
    centerOfMass +=
        (CurrAtom->getMass() * CurrAtom->Coordinates2Eigen(options));
    CurrAtom = CurrAtom->getNext();
  }
  centerOfMass /= parent->getMolecularMass();
  CurrAtom = NULL;
  return centerOfMass;
}

void
Structure::Shift2CenterOfMass(enum StructureOptions options)
{
  centerOfMass = -1.0 * getCenterOfMass(options);
  parent->setInputTransformation(centerOfMass);
  Atom *CurrAtom = HeadAtom;
  while (CurrAtom)
  {
    CurrAtom->shiftCoordinates(centerOfMass, options);
    CurrAtom = CurrAtom->getNext();
  }
  centerOfMass = getCenterOfMass(options);
  CurrAtom     = NULL;
}

Eigen::MatrixXd
Structure::Rotate2InertiaPAS(enum StructureOptions options)
{
  if (!((options & StructureOptions::Initial) ||
        (options & StructureOptions::Optimized)))
  {
    std::cerr << "ERROR:\nUnknown type of coordinates requested in "
                 "Rotate2InertiaPAS...\n\t\tUsing initial coordinates...\n";
    options = StructureOptions::Initial;
  }

  Eigen::Matrix3d Rotation =
      CalculateInertiaTensor(options)->getEvecs().transpose();
  parent->setInputTransformation(Rotation);
  parent->setInputInertiaTensor(initialInertiaTensor->getTensor());
  parent->setInputInertiaTensorEigVal(initialInertiaTensor->getEvals());
  parent->setInputInertiaTensorEigVec(initialInertiaTensor->getEvecs());
  Atom *CurrAtom;
  for (CurrAtom = HeadAtom; CurrAtom; CurrAtom = CurrAtom->getNext())
  {
    CurrAtom->rotateCoordinates(Rotation, options);
  }
  CalculateInertiaTensor(options);
  return Rotation;
}

Tensor *
Structure::CalculateInertiaTensor(enum StructureOptions options)
{
  Shift2CenterOfMass(options);
  Atom *CurrAtom;
  Coordinates *C;
  double m;
  Eigen::Matrix3d ITensor = Eigen::Matrix3d::Zero(3, 3);

  if (!((options & StructureOptions::Initial) ||
        (options & StructureOptions::Optimized)))
  {
    std::cerr
        << "ERROR:\tUnknown type of coordinates requested for "
           "CalculateInertiaTensor...\n\t\tUsing initial coordinates...\n";
    options = StructureOptions::Initial;
  }

  for (CurrAtom = HeadAtom; CurrAtom; CurrAtom = CurrAtom->getNext())
  {
    C = CurrAtom->getCoordinates(options);
    m = CurrAtom->getMass();
    ITensor(STANDARD_AXIS_X_, STANDARD_AXIS_X_) +=
        (m * (pow(C->y, 2.0) + pow(C->z, 2.0)));
    ITensor(STANDARD_AXIS_Y_, STANDARD_AXIS_Y_) +=
        (m * (pow(C->x, 2.0) + pow(C->z, 2.0)));
    ITensor(STANDARD_AXIS_Z_, STANDARD_AXIS_Z_) +=
        (m * (pow(C->x, 2.0) + pow(C->y, 2.0)));

    ITensor(STANDARD_AXIS_X_, STANDARD_AXIS_Y_) -= (m * C->x * C->y);
    ITensor(STANDARD_AXIS_X_, STANDARD_AXIS_Z_) -= (m * C->x * C->z);
    ITensor(STANDARD_AXIS_Y_, STANDARD_AXIS_Z_) -= (m * C->y * C->z);
  }
  ITensor(STANDARD_AXIS_Y_, STANDARD_AXIS_X_) =
      ITensor(STANDARD_AXIS_X_, STANDARD_AXIS_Y_);
  ITensor(STANDARD_AXIS_Z_, STANDARD_AXIS_X_) =
      ITensor(STANDARD_AXIS_X_, STANDARD_AXIS_Z_);
  ITensor(STANDARD_AXIS_Z_, STANDARD_AXIS_Y_) =
      ITensor(STANDARD_AXIS_Y_, STANDARD_AXIS_Z_);

  C        = NULL;
  CurrAtom = NULL;

  if (options & StructureOptions::Initial)
  {
    if (!initialInertiaTensor)
      initialInertiaTensor = new Tensor(ITensor);
    else
      initialInertiaTensor->setTensor(ITensor);
    return initialInertiaTensor;
  }
  else // if ( options & StructureOptions::Optimized )
  {
    if (!optimizedInertiaTensor)
    {
      optimizedInertiaTensor = new Tensor(ITensor);
    }
    else
      optimizedInertiaTensor->setTensor(ITensor);
    return optimizedInertiaTensor;
  }
}

Tensor *
Structure::getInertiaTensor(StructureOptions options)
{
  if (options & StructureOptions::Initial)
  {
    return initialInertiaTensor;
  }
  else
  {
    return optimizedInertiaTensor;
  }
}

Eigen::MatrixXd
Structure::getRDCDistances(StructureOptions options)
{
  unsigned int NOR    = parent->getNOR();
  Eigen::MatrixXd dis = Eigen::MatrixXd::Zero(NOR, NOR);

  SphericalHarmonics *Y;

  unsigned int i;

  for (i = 0, Y = HeadYmatrix->getHeadHarmonic(); Y; Y = Y->getNext(), ++i)
  {
    dis(i, i) = Y->getRDC()->getDistance(options);
  }
  return dis;
}

Eigen::MatrixXd
Structure::getDistances(StructureOptions options)
{
  unsigned int NOA    = parent->getNOA();
  Eigen::MatrixXd dis = Eigen::MatrixXd::Zero(NOA, NOA);

  unsigned int a, b;
  Atom *A, *B;

  for (a = 0, A = HeadAtom; A && a < NOA; A = A->getNext(), ++a)
  {
    for (b = (a + 1), B = A->getNext(); B && b < NOA; B = B->getNext(), ++b)
    {
      dis(a, b) = A->getDistance(B, options);
    }
  }
  return dis;
}

unsigned int
Structure::getNumberOfLongRangeRDCs()
{
  if (parent)
    return parent->getNumberOfLongRangeRDCs();
  else
    return 0;
}

unsigned int
Structure::getNOR()
{
  if (parent)
    return parent->getNOR();
  else
    return 0;
}

void
Structure::retainCoordinates()
{
  Atom *CurrAtom;
  for (CurrAtom = HeadAtom; CurrAtom; CurrAtom = CurrAtom->getNext())
  {
    CurrAtom->retainCoordinates();
  }
}

void
Structure::copyCoordinates(StructureOptions CurrOpt,
                           Structure *RefStruc,
                           StructureOptions RefOpt)
{
  Atom *CurrAtom, *RefAtom;
  for (CurrAtom = HeadAtom, RefAtom = RefStruc->HeadAtom; CurrAtom && RefAtom;
       CurrAtom = CurrAtom->getNext(), RefAtom = RefAtom->getNext())
  {
    CurrAtom->setCoordinates(RefAtom->Coordinates2Eigen(RefOpt), CurrOpt);
  }
}

void
Structure::setUnfixed()
{
  Atom *CurrAtom;
  for (CurrAtom = HeadAtom; CurrAtom; CurrAtom = CurrAtom->getNext())
  {
    CurrAtom->setUnfixed();
  }
}

int
Structure::run_structure_MC(BasicInformation &baseInformation, Flags &flags)
{
  int steps = 0, i;
  MC_Run_Information tmp;
  tmp.console_output = false;
  tmp.convergence    = 0.05;
  tmp.max_steps      = 5000;
  if (baseInformation.phobos_workspace_Sphericals[0] == NULL)
  {
    for (i = 0; i < baseInformation.numOfThreadsUsed; ++i)
    {
#pragma omp task firstprivate(i) shared(baseInformation)
      baseInformation.phobos_workspace_Sphericals[i] = (double *) malloc(
          baseInformation.phobos_worksize_Sphericals * sizeof(double));
    }
  }
  if (prevStruc)
  {
#pragma omp task
    get_radius_of_Gyration(StructureOptions::Optimized);
#pragma omp task
    prevStruc->get_radius_of_Gyration(StructureOptions::Optimized);
#pragma omp taskwait
    if ((radiusOfGyration / prevStruc->radiusOfGyration) > 1e3)
    {
      flags.monteCarloBootstrapping = false;
      flags.monteCarloOutput        = false;
      flags.plotMonteCarlo          = false;
      baseInformation.errorKey +=
          ("ERROR:\tRadius of gyration expanded from " +
           std::to_string(prevStruc->radiusOfGyration) + " Ang^2 to " +
           std::to_string(radiusOfGyration) + " Ang^2.\n");
      baseInformation.state = ERROR_STRUCTURE_EXPLOSION;
      return 0;
    }
  }

  monte_carlo_bootstrap = new MonteCarloStatistic();
  steps = monte_carlo_bootstrap->startMonteCarlo(*parent, this, baseInformation,
                                                 flags, tmp);
  monte_carlo_output = monte_carlo_bootstrap->getStatistics();
#pragma omp task shared(monte_carlo_bootstrap)
  monte_carlo_bootstrap->cleanup();
  return steps;
}

int
Structure::check_Stop(BasicInformation &baseInformation, Flags &flags)
{
  // Check Q-factor difference.
  if (!prevStruc || !flags.monteCarloBootstrapping)
    return 0;

  baseInformation.MC_stop_reason = baseInformation.MC_stop_reason |
                                   this->check_Q_diff(baseInformation, flags);
  baseInformation.MC_stop_reason =
      baseInformation.MC_stop_reason | this->check_MC_p(baseInformation);
  baseInformation.MC_stop_reason =
      baseInformation.MC_stop_reason | this->check_MC_A(baseInformation);

  if (baseInformation.MC_stop_reason == MC_stop::undefined)
    return 0;
  else if (flags.converged)
    return 0;
  else
  {
    flags.converged                         = true;
    baseInformation.structureLabelExtension = "_overOpt_";
    baseInformation.over_titania_iterations =
        baseInformation.numOfOptSteps + baseInformation.overOptimization;
    return 1;
  }
}

MC_Output
Structure::get_monte_carlo_output()
{
  return monte_carlo_output;
}

MC_stop
Structure::check_Q_diff(BasicInformation &baseInformation, Flags &flags)
{
  if (!prevStruc)
    return MC_stop::undefined;
  Eigen::MatrixXd rdcCalc;
  Eigen::MatrixXd rdcMatrix = parent->getRDCmatrix();
  Eigen::MatrixXd w;
  double dNOS = ((double) parent->getNORsets());
  if (flags.calculateFullMatrix)
    w = Eigen::MatrixXd::Identity(parent->getNOR(), parent->getNOR());
  else
    w = parent->getWmatrix();
  Eigen::VectorXd QFac =
      this->Qfacs(rdcMatrix, rdcCalc, StructureOptions::Optimized, w, flags);
  Eigen::VectorXd preFac = prevStruc->Qfacs(
      rdcMatrix, rdcCalc, StructureOptions::Optimized, w, flags);
  Eigen::VectorXd Q_diff                         = QFac - preFac;
  baseInformation.stop_crit.Q_factor_convergence = Q_diff.norm() / sqrt(dNOS);
  if (baseInformation.stop_crit.Q_factor_convergence <
      baseInformation.limits.Q_factor_convergence)
    return MC_stop::QFac;
  else
    return MC_stop::undefined;
}

MC_stop
Structure::check_MC_A(BasicInformation &baseInformation)
{
  MC_stop out = MC_stop::undefined;
  if (!prevStruc)
    return out;
  Eigen::MatrixXd A_mean_diff =
      this->monte_carlo_output.Aligns - prevStruc->monte_carlo_output.Aligns;
  Eigen::MatrixXd A_sigm_diff = this->monte_carlo_output.Aligns_sigm -
                                prevStruc->monte_carlo_output.Aligns_sigm;

  const unsigned int NOS = parent->getNORsets();
  double dNOS            = ((double) NOS);

  baseInformation.stop_crit.alignment_mean_convergence =
      A_mean_diff.norm() / sqrt(5.0 * dNOS);
  baseInformation.stop_crit.alignment_sigm_convergence =
      A_sigm_diff.norm() / sqrt(5.0 * dNOS);

  if (baseInformation.stop_crit.alignment_mean_convergence <
      baseInformation.limits.alignment_mean_convergence)
    out = out | MC_stop::AMean; // percentage mean change
  if (baseInformation.stop_crit.alignment_sigm_convergence <
      baseInformation.limits.alignment_sigm_convergence)
    out = out | MC_stop::ASigm; // percentage mean change

  return out;
}

MC_stop
Structure::check_MC_p(BasicInformation &baseInformation)
{
  MC_stop out = MC_stop::undefined;
  if (!prevStruc)
    return out;
  Eigen::MatrixXd p_mean_diff = this->monte_carlo_output.p_trig_mean -
                                prevStruc->monte_carlo_output.p_trig_mean;
  Eigen::MatrixXd p_sigm_diff =
      this->monte_carlo_output.p_sigm - prevStruc->monte_carlo_output.p_sigm;
  Eigen::MatrixXd R_mean_diff =
      this->monte_carlo_output.R_mean - prevStruc->monte_carlo_output.R_mean;

  const unsigned int NOR = parent->getNOR();
  double dNOR            = ((double) NOR);


  baseInformation.stop_crit.sphericals_mean_convergence =
      p_mean_diff.norm() / sqrt(4.0 * dNOR);
  baseInformation.stop_crit.sphericals_sigm_convergence =
      p_sigm_diff.norm() / sqrt(2.0 * dNOR);
  baseInformation.stop_crit.sphericals_spread_convergence =
      R_mean_diff.norm() / sqrt(dNOR);

  if (baseInformation.stop_crit.sphericals_mean_convergence <
      baseInformation.limits.sphericals_mean_convergence)
    out = out | MC_stop::pMean; // percentage mean change
  if (baseInformation.stop_crit.sphericals_sigm_convergence <
      baseInformation.limits.sphericals_sigm_convergence)
    out = out | MC_stop::pSigm; // percentage mean change
  if (baseInformation.stop_crit.sphericals_spread_convergence <
      baseInformation.limits.sphericals_spread_convergence)
    out = out | MC_stop::RMean;

  return out;
}

int
Structure::initializeRDCindex(Molecule &CurrMol)
{
  const unsigned int NOA = CurrMol.getNOA();
  unsigned int i, j, bestVal, testVal, currVal, currI;

  Atom *CurrAtom, *TestAtom, *BestAtom;

  Atom **BestList = (Atom **) malloc(NOA * sizeof(Atom *));
  for (i = 0; i < NOA; ++i)
    BestList[i] = NULL;

  CurrAtom = HeadAtom;
  bestVal = testVal = 0;
  currI             = 0;
  while (CurrAtom)
  {
    CurrAtom->setUnfixed();
    if (CurrAtom->getZ() == 1)
    {
      CurrAtom = CurrAtom->getNext();
      continue;
    }
    currVal = CurrAtom->getNumberOfHarmonics();
    if (currVal > bestVal)
    {
      bestVal     = currVal;
      BestList[0] = CurrAtom;
    }
    else if (currVal && currVal == bestVal)
    {
      currVal = 0;
      testVal = 0;
      for (i = 0; i < CurrAtom->getNumberOfBonds(); ++i)
      {
        currVal += CurrAtom->getBond(i)
                       ->getBondpartner(CurrAtom)
                       ->getNumberOfHarmonics();
      }
      for (i = 0; i < BestList[0]->getNumberOfBonds(); ++i)
      {
        testVal += BestList[0]
                       ->getBond(i)
                       ->getBondpartner(BestList[0])
                       ->getNumberOfHarmonics();
      }
      if (currVal > testVal)
      {
        bestVal     = CurrAtom->getNumberOfHarmonics();
        BestList[0] = CurrAtom;
      }
    }
    CurrAtom = CurrAtom->getNext();
  }

  currI = 0;

  BestAtom = CurrAtom = BestList[0];
  CurrAtom->setFixed();
  bestVal = 0;

  for (currI = 1; currI < NOA; ++currI)
  {
    bestVal = testVal = 0;
    BestAtom          = NULL;

    for (j = 0; j < currI; ++j)
    {
      CurrAtom = BestList[j];
      for (i = 0; i < CurrAtom->getNumberOfBonds(); ++i)
      {
        TestAtom = CurrAtom->getBond(i)->getBondpartner(CurrAtom);
        if (TestAtom->getZ() == 1)
          continue; // First prefere Atoms that build the backbone
        if (TestAtom->isFixed())
          continue;
        currVal = TestAtom->getNumberOfHarmonics();
        if (currVal > bestVal)
        {
          BestAtom = TestAtom;
          bestVal  = currVal;
        }
      }
    }
    if (BestAtom == NULL)
    {
      for (j = 0; j < currI; ++j)
      {
        CurrAtom = BestList[j];
        for (i = 0; i < CurrAtom->getNumberOfBonds(); ++i)
        {
          TestAtom = CurrAtom->getBond(i)->getBondpartner(CurrAtom);
          if (TestAtom->isFixed())
            continue;
          currVal = TestAtom->getNumberOfHarmonics();
          if (currVal >= bestVal)
          {
            BestAtom = TestAtom;
            bestVal  = currVal;
          }
        }
      }
    }
    BestList[currI] = BestAtom;
    BestAtom->setFixed();
  }

  for (i = 0; i < NOA; ++i)
  {
    BestList[i]->setUnfixed();
    BestList[i]->setrdcIndex(i);
    BestList[i] = NULL;
    delete BestList[i];
  }

  CurrAtom = TestAtom = BestAtom = NULL;
  delete CurrAtom;
  delete TestAtom;
  delete BestAtom;
  delete BestList;
  return 0;
}

void
Structure::deleteBonds()
{
  Atom *CurrAtom, *BondAtom;
  Bond *CurrBond;
  unsigned int i, j;
  for (CurrAtom = HeadAtom; CurrAtom;
       CurrAtom = CurrAtom->getNext()) // Go for all atoms
  {
    for (i = 0; i < CurrAtom->getNumberOfBonds(); ++i) // Cycle all bonds
    {
      CurrBond = CurrAtom->getBond(i);
      if (CurrBond == NULL)
        continue;                                    // Check if it still exists
      BondAtom = CurrBond->getBondpartner(CurrAtom); // Go for the bonded atom
      for (j = 0; j < BondAtom->getNumberOfBonds();
           ++j) // Cycle all bond of the partner
      {
        if (CurrBond == BondAtom->getBond(j))
        {
          BondAtom->setBond(NULL, j); // Since it still has to exist
          break;                      // at this point, unlink the bond
        }                             // and stop searching
        else
          continue;
      }
      delete CurrBond;            // Delete the bond
      CurrAtom->setBond(NULL, i); // and unlink it at the current atom
    }
  }
  CurrAtom = NULL;
  BondAtom = NULL;
  CurrBond = NULL;
}

Bond **
Structure::initializeListOfBonds()
{
  ListOfBonds = (Bond **) malloc(NumberOfBonds * sizeof(Bond *));
  unsigned int b;
  for (b = 0; b < NumberOfBonds; ++b)
    ListOfBonds[b] = NULL;
  return ListOfBonds;
}

// void count_distance_violations ( StructureOptions );

void
Structure::check_inversion(StructureSimulator &MainSimulation,
                           BasicInformation &baseInformation,
                           StructureOptions options)
{
  /*
   *  To early use of the inversion leads to
   *  inconsistencies in the algorithm.
   */
  if (baseInformation.numOfOptSteps < 10)
    return;
  Atom *CurrAtom, *MaxAtom;
  double CurrPot, MaxPot;

  int min_violations = (baseInformation.lowerInversionAfter <= 0 ?
                            5 :
                            (baseInformation.numOfOptSteps <
                                     baseInformation.lowerInversionAfter ?
                                 5 :
                                 10));
  unsigned int inversionCenter, min_index = getNOA();
  CurrPot = MaxPot = .0;
  for (CurrAtom = HeadAtom, MaxAtom = NULL; CurrAtom;
       CurrAtom = CurrAtom->getNext())
  {
    if (CurrAtom->get_number_of_distance_violations() >= min_violations)
    {
      CurrPot = ((double) CurrAtom->get_number_of_distance_violations() *
                 CurrAtom->calculateVdW_Potential(MainSimulation, options));
      if (CurrAtom->is_invertable())
        for (unsigned int i = 0; i < CurrAtom->getNumberOfBonds(); ++i)
          CurrPot += (CurrAtom->getBondpartner(i)->getA() == 1 ?
                          ((double) CurrAtom->getBondpartner(i)
                               ->get_number_of_distance_violations() *
                           CurrAtom->getBondpartner(i)->calculateVdW_Potential(
                               MainSimulation, options)) :
                          .0);
      if (CurrPot > MaxPot)
      {
        MaxAtom = CurrAtom;
        MaxPot  = CurrPot;
      }
    }
  }
  if (MaxAtom && MaxAtom->getA() == 1 &&
      MaxAtom->getBondpartner(0)->is_invertable())
  {
    MaxAtom = MaxAtom->getBondpartner(0);
    for (unsigned int i = 0; i < MaxAtom->getNumberOfBonds(); ++i)
    {
      if (MaxAtom->getBondpartner(i)->getZ() == 1)
        MaxAtom->getBondpartner(i)->reset_distance_violations();
    }
  }
  for (unsigned int i = inversionCenter = 0;
       MaxAtom && i < MaxAtom->getNumberOfBonds(); ++i)
  {
    if (min_index < MaxAtom->getBondpartner(i)->getrdcIndex())
    {
      min_index       = MaxAtom->getBondpartner(i)->getrdcIndex();
      inversionCenter = i;
    }
  }
  if (MaxAtom)
  {
    Eigen::Vector3d saved_coord = MaxAtom->Coordinates2Eigen(options);
    MaxAtom->setCoordinates((2.0 * MaxAtom->getBondpartner(inversionCenter)
                                       ->Coordinates2Eigen(options) -
                             MaxAtom->Coordinates2Eigen(options)),
                            options);
    CurrPot = ((double) MaxAtom->get_number_of_distance_violations() *
               MaxAtom->calculateVdW_Potential(MainSimulation, options));
    if (MaxAtom->is_invertable())
      for (unsigned int i = 0; i < MaxAtom->getNumberOfBonds(); ++i)
        CurrPot += (MaxAtom->getBondpartner(i)->getA() == 1 ?
                        ((double) MaxAtom->getBondpartner(i)
                             ->get_number_of_distance_violations() *
                         MaxAtom->getBondpartner(i)->getVdW_Potential()) :
                        .0);
    if (CurrPot > MaxPot)
    {
      MaxAtom->setCoordinates(saved_coord, options);
      MaxAtom->reset_distance_violations();
    }
    else
      for (CurrAtom = HeadAtom; CurrAtom; CurrAtom = CurrAtom->getNext())
        CurrAtom->reset_distance_violations();
  }
}

void
Structure::toxyz(StructureSimulator &MainSimulation,
                 RedundantInternals &MainRedundants,
                 BasicInformation &baseInformation,
                 Flags &flags,
                 StructureOptions options)
{
  addTimeSlot(baseInformation, baseInformation.SCRMtime);

  if (prevStruc)
    MainSimulation.newStructure(this, options);
  addTimeSlot(baseInformation, baseInformation.Systemtime);

  calculateVdW_Potential(MainSimulation, StructureOptions::Initial);
  addTimeSlot(baseInformation, baseInformation.MMFF94time);

  if (options == StructureOptions::Optimized)
    retainCoordinates();

  // If redundant internal coordinates took a very bad step
  // turn them off for this iteration step
  if (baseInformation.state == ERROR_INVALID_REDUNDANTS_STEP)
    flags.useRedundants = false;

  /*
   * Run the full "non-redundant" algorithm if you are still in an
   * requesting this behaviour or redundant internal coordinates
   * are turned off (on this run or this iteration[see above]).
   */
  if (baseInformation.numOfOptSteps <= baseInformation.useRedundantsOnlyAfter ||
      !flags.useRedundants)
  {
    Vector2Structure(options);
    rmsd2Structure(options);
  }
  addTimeSlot(baseInformation, baseInformation.structureTime);

  if (flags.useRedundants)
  {
    // Check the number of distance violations and invert
    // the most violoating atom(-group).
    if ((baseInformation.redundants_distance_optimization == FULL_DISTANCES_) ||
        (baseInformation.redundants_distance_optimization ==
         INVERSION_DISTANCES_))
      check_inversion(MainSimulation, baseInformation, options);

    // Use the redundant internal coordinates to generate
    // new coordinates.
    MainRedundants.newStructure(this, &MainSimulation,
                                StructureOptions::Optimized);
    //      RedundantInternals Red( &MainSimulation, baseInformation, flags );
    baseInformation.state = MainRedundants.S2x(baseInformation, flags);

    /*
     * If redundant internal coordinates took a very bad step
     * turn them off for this iteration step (see above) and
     * restart the structure generation.
     */
    if (baseInformation.state == ERROR_INVALID_REDUNDANTS_STEP)
    {
      addTimeSlot(baseInformation, baseInformation.structureTime);
      toxyz(MainSimulation, MainRedundants, baseInformation, flags, options);
    }

    // Do the eckart tansformation
    if (!flags.skipEckart)
    {
      Rotate_2_Eckart_Frame(parent->getHeadStruc(), this, baseInformation,
                            options);
      if (baseInformation.state == ERROR_ECKART_TRANSFORMATION)
      {
        flags.skipEckart      = true;
        baseInformation.state = GOOD_STATE;
      }
      else
        outputXYZ(HeadAtom, baseInformation, StructureOptions::Optimized,
                  "redundants-in-Eckart");
    }

    /*
     * Recalculate all polar angles since they might
     * have changed due to the  Eckart transformation
     * what will lead to undefined behaviour in
     * the following optimization steps.
     */
    for (SphericalHarmonics *CurrHarmonic = HeadYmatrix->getHeadHarmonic();
         CurrHarmonic; CurrHarmonic       = CurrHarmonic->getNext())
    {
#pragma omp task firstprivate(CurrHarmonic)
      CurrHarmonic->recalculateAngles(StructureOptions::Optimized);
    }
  }
  else
  {
    MainSimulation.startSimplexMinimizer(
        (MainSimulation.getNumberToOpt() * 1000), true);
    addTimeSlot(baseInformation, baseInformation.MMFF94time);

    for (Atom *CurrAtom = HeadAtom; CurrAtom; CurrAtom = CurrAtom->getNext())
    {
      if (CurrAtom->getNumberOfHarmonics())
      {
#pragma omp task firstprivate(CurrAtom) shared(MainSimulation)
        MainSimulation.AtomEnergy(CurrAtom, options);
      }
    }
#pragma omp taskwait
    if (baseInformation.numOfOptSteps > 1 &&
        !maxEnergyAtom(HeadAtom, MainSimulation))
    {
      MainSimulation.startSimplexMinimizer(
          (MainSimulation.getNumberToOpt() * 1000), true);
      addTimeSlot(baseInformation, baseInformation.MMFF94time);
    }
#pragma omp task shared(HeadAtom, baseInformation, options)
    outputXYZ(HeadAtom, baseInformation, options, "MMFF94-opt");

    /*
     * If redundant internal coordinates took a very bad step
     * and were turned off for this iteration step (see above),
     * they will be turned on again after using the "non-redundant"
     * algorithm on this step.
     */
    if (baseInformation.state == ERROR_INVALID_REDUNDANTS_STEP)
      flags.useRedundants = true;
  }
#pragma omp taskwait

  addTimeSlot(baseInformation, baseInformation.structureTime);
  TODO(Check if this van der Waals potential has to be
           calculated.Especially using initial coordinates.)
  /*   calculateVdW_Potential ( MainSimulation, StructureOptions::Initial );
     addTimeSlot( baseInformation, baseInformation.MMFF94time );*/
}

void
Structure::calculateVdW_Potential(StructureSimulator &MainSim,
                                  StructureOptions opt)
{
  for (Atom *CurrAtom = HeadAtom; CurrAtom; CurrAtom = CurrAtom->getNext())
  {
    CurrAtom->calculateVdW_Potential(MainSim, opt);
  }
}

void
Structure::Vector2Structure(StructureOptions options)
{
  unsigned int a, b;
  double rnext, aopt;

  Atom *CurrAtom, *BondAtom;
  Bond *CurrBond;
  Eigen::MatrixXd central, preC, currC, bo;
  Eigen::Vector3d vecPre, vecBond;

  setUnfixed();
  const unsigned int NOA = getNOA();

  CurrAtom = HeadAtom;
  atomByrdcIndex(&CurrAtom, 0);

  CurrAtom->setFixed();
  CurrAtom = HeadAtom;

  for (a = 0; a < NOA; ++a)
  {
    atomByrdcIndex(&CurrAtom, a);
    if (!CurrAtom->isFixed())
      continue;
    central = CurrAtom->Coordinates2Eigen(options);
    for (b = 0; b < CurrAtom->getNumberOfBonds(); ++b)
    {
      CurrBond = CurrAtom->getBond(b);
      BondAtom = CurrBond->getBondpartner(CurrAtom);
      if (BondAtom->isFixed())
        continue;
      else if (!CurrBond->getHarmonic())
      {
        BondAtom->setFixed();
        continue;
      }
      preC    = BondAtom->Coordinates2Eigen(StructureOptions::Initial);
      vecBond = Polar2Eigen(CurrBond->getHarmonic()->getTheta(options),
                            CurrBond->getHarmonic()->getPhi(options));

      vecPre = preC - central;
      vecPre = vecPre / vecPre.norm();

      rnext = CurrBond->getLength();
      aopt  = vecBond.transpose() * vecPre;

      if (aopt < .0)
      {
        vecBond *= -1.0;
      }

      vecBond = central + rnext * vecBond;
      BondAtom->setCoordinates(vecBond, options);

      BondAtom->setFixed();
    }
  }
  CurrAtom = NULL;
  BondAtom = NULL;
  CurrBond = NULL;
}

void
Structure::rmsd2Structure(StructureOptions options)
{
  unsigned int a, b, preA, runs;
  const unsigned int NOA = getNOA();
  double r, rmsd1, rmsd2, mean, sigma, rezNORA, maxrmsd, maxscat;

  rezNORA = (1.0 / ((double) parent->getNORA()));

  Eigen::VectorXd rmsd;

  Atom *CurrAtom = HeadAtom;
  Atom *partner;
  SphericalHarmonics *CurrHarmonic;
  Eigen::Vector3d CurrCoord, PartnerCoord, HarmonicsVec;

  CurrAtom = getMaxStructure2rmsdAtom(NOA, a, rmsd, options);
  preA     = a;

  mean = rmsd.mean();
  for (sigma = 0, a = 0; a < NOA; ++a)
  {
    if (rmsd(a) > .0) // Just monitor those who contribute to rdc vectors
    {
      sigma += (pow((rmsd(a) - mean), 2.0));
    }
  }
  sigma = sqrt(sigma * rezNORA);

  maxscat = mean - 0.5 * sigma;
  runs    = 0;
  while (true)
  {
    rmsd.maxCoeff(&a);
    atomByrdcIndex(&CurrAtom, a);

    for (b = 0; b < CurrAtom->getNumberOfHarmonics(); ++b)
    {
      if (CurrAtom->getNumberOfHarmonics() == 1)
        break;
      CurrHarmonic = CurrAtom->getHarmonic(b);
      if (CurrHarmonic->getRange() == 1)
      {
        if (CurrAtom->getHarmonic(b)->getAtom1() == CurrAtom)
          partner = CurrAtom->getHarmonic(b)->getAtom2();
        else
          partner = CurrAtom->getHarmonic(b)->getAtom1();

        PartnerCoord = partner->Coordinates2Eigen(options);
        HarmonicsVec = Polar2Eigen(CurrHarmonic->getTheta(options),
                                   CurrHarmonic->getPhi(options));
        r            = CurrAtom->getBond(partner)->getLength();

        CurrCoord = PartnerCoord + r * HarmonicsVec;
        CurrAtom->setCoordinates(CurrCoord, options);
        rmsd1 = structure2rmsd(options, a);

        CurrCoord = PartnerCoord - r * HarmonicsVec;
        CurrAtom->setCoordinates(CurrCoord, options);
        rmsd2 = structure2rmsd(options, a);


        if (rmsd1 < rmsd2)
        {
          CurrCoord = PartnerCoord + r * HarmonicsVec;
          CurrAtom->setCoordinates(CurrCoord, options);
        }
        break;
      }
    }
    rmsd(preA) = .0;
    maxrmsd    = rmsd.maxCoeff(&a);
    if (preA == a)
      break;
    preA = a;
    if (maxrmsd < maxscat)
      break;
    if (runs++ > parent->getNORA())
      break;
  }

  CurrAtom     = NULL;
  partner      = NULL;
  CurrHarmonic = NULL;
}

Eigen::VectorXd
Structure::Qfacs(Eigen::MatrixXd &rdcs,
                 Eigen::MatrixXd &rdcCalc,
                 StructureOptions opt,
                 Eigen::MatrixXd &w,
                 Flags &flags)
{
  Eigen::MatrixXd C    = HeadYmatrix->determineBmatrix(*parent, this, opt);
  Eigen::MatrixXd Dmax = getDmaxMatrix(opt);

  return (getQfacs(C, rdcs, rdcCalc, Dmax, w, flags));
}

Atom *
Structure::getMaxStructure2rmsdAtom(unsigned int NOA,
                                    unsigned int &a,
                                    Eigen::VectorXd &rmsd,
                                    StructureOptions options)
{
  Atom *CurrAtom = HeadAtom;
  rmsd           = Eigen::VectorXd::Zero(NOA);
  for (a = 0; a < NOA; ++a)
    rmsd(a) = structure2rmsd(options, a);
  rmsd.maxCoeff(&a);
  atomByrdcIndex(&CurrAtom, a);
  return CurrAtom;
}

double
Structure::structure2rmsd(enum StructureOptions options, unsigned int a)
{
  double rmsd = .0;
  unsigned int s;
  Atom *A1;
  Eigen::Vector3d C1, R, dH;

  SphericalHarmonics *CurrHarmonic;

  A1 = HeadAtom;
  atomByrdcIndex(&A1, a);

  C1 = A1->Coordinates2Eigen(options);

  for (s = 0; s < A1->getNumberOfHarmonics(); ++s)
  {
    CurrHarmonic = A1->getHarmonic(s);
    R  = C1 - CurrHarmonic->getPartner(A1)->Coordinates2Eigen(options);
    R  = R / R.norm();
    dH = Polar2Eigen(CurrHarmonic, options);

    if (R.transpose() * dH < .0)
      dH *= -1.0;
    R -= dH;
    rmsd += R.norm();
  }

  if (rmsd)
    rmsd = rmsd / ((double) A1->getNumberOfHarmonics());

  A1           = NULL;
  CurrHarmonic = NULL;
  if (rmsd)
    return rmsd;
  else
    return -1.0;
}

void
Structure::determineCosineMatrix(StructureOptions options)
{
  const double sr3 = sqrt(3.0);

  Atom *A1 = HeadAtom, *A2 = HeadAtom;
  Eigen::Vector3d V;
  RDCdata *CurrRDC;
  Eigen::MatrixXd cosineMatrix = Eigen::MatrixXd::Zero(parent->getNOR(), 5);

  unsigned int i;
  for (i = 0, CurrRDC = parent->getHeadSet()->getHeadData(); CurrRDC;
       ++i, CurrRDC   = CurrRDC->getNext())
  {
    atomByrdcIndex(&A1, CurrRDC->getAtom1()->getrdcIndex());
    atomByrdcIndex(&A2, CurrRDC->getAtom2()->getrdcIndex());

    V = A1->Coordinates2Eigen(options) - A2->Coordinates2Eigen(options);
    V = V / V.norm();

    cosineMatrix(i, SAUPE_ZZ_) =
        (3.0 * pow(V(STANDARD_AXIS_Z_), 2.0) - 1.0) / 2.0;
    cosineMatrix(i, SAUPE_XX_YY_) =
        sr3 * (pow(V(STANDARD_AXIS_X_), 2.0) - pow(V(STANDARD_AXIS_Y_), 2.0)) /
        2.0;
    cosineMatrix(i, SAUPE_XY_) =
        sr3 * V(STANDARD_AXIS_X_) * V(STANDARD_AXIS_Y_);
    cosineMatrix(i, SAUPE_XZ_) =
        sr3 * V(STANDARD_AXIS_X_) * V(STANDARD_AXIS_Z_);
    cosineMatrix(i, SAUPE_YZ_) =
        sr3 * V(STANDARD_AXIS_Y_) * V(STANDARD_AXIS_Z_);
  }
  if (options & StructureOptions::Initial)
  {
    iniCosineMatrix = cosineMatrix;
  }
  else if (options & StructureOptions::Optimized)
  {
    optCosineMatrix = cosineMatrix;
  }
  else
  {
    std::cerr << "ERROR:\tUnkonwn type of StructureOptions requested in "
                 "determineCosineMatrix of structure "
              << label
              << "...\n"
                 "\t\tThe request will be skipped...\n";
  }
  A1      = NULL;
  A2      = NULL;
  CurrRDC = NULL;
}

Eigen::MatrixXd
Structure::getCosineMatrix(StructureOptions options)
{
  unsigned int NOR = parent->getNOR();

  if (options & StructureOptions::Initial)
  {
    if ((iniCosineMatrix.rows() != NOR) ||
        (iniCosineMatrix == Eigen::MatrixXd::Zero(NOR, 5)))
    {
      determineCosineMatrix(options);
    }
    return iniCosineMatrix;
  }
  else if (options & StructureOptions::Optimized)
  {
    if ((optCosineMatrix.rows() != NOR) ||
        (optCosineMatrix == Eigen::MatrixXd::Zero(NOR, 5)))
    {
      determineCosineMatrix(options);
    }
    return optCosineMatrix;
  }
  // This statement can just be reached if StructureOption is unknown
  std::cerr << "ERROR:\tUnknown type of StructureOptions was requested for the "
               "cosine matrix of "
            << label << "...\n"
            << "\t\tReturning null matrix...\n";
  return Eigen::MatrixXd::Zero(NOR, COSINE_ELEMENTS_);
}

void
Structure::set_rdc_vector_sampling(double r, double a, StructureOptions opt)
{
  int index                   = 0;
  index                       = (opt == StructureOptions::Initial ? 0 : 1);
  rdcVectorSampling(index, 0) = r;
  rdcVectorSampling(index, 1) = a;
}

void
Structure::get_rdc_vector_sampling(double &r, double &a, StructureOptions opt)
{
  int index = 0;
  index     = (opt == StructureOptions::Initial ? 0 : 1);
  r         = rdcVectorSampling(index, 0);
  a         = rdcVectorSampling(index, 1);
}

double
Structure::get_rdc_vector_sampling_r(StructureOptions opt)
{
  return (opt == StructureOptions::Initial ? rdcVectorSampling(0, 0) :
                                             rdcVectorSampling(1, 0));
}

double
Structure::get_rdc_vector_sampling_a(StructureOptions opt)
{
  return (opt == StructureOptions::Initial ? rdcVectorSampling(0, 1) :
                                             rdcVectorSampling(1, 1));
}

double
Structure::calculateAllAtomRMSDofStep()
{
  return calculateAllAtomRMSD(
      coordinates2EigenVector(StructureOptions::Initial),
      coordinates2EigenVector(StructureOptions::Optimized));
}
double
Structure::calculateAllAtomRMSD2Input()
{
  return calculateAllAtomRMSD(
      parent->getHeadStruc()->coordinates2EigenVector(
          StructureOptions::Initial),
      coordinates2EigenVector(StructureOptions::Optimized));
}

Eigen::VectorXd
Structure::coordinates2EigenVector(StructureOptions options)
{
  unsigned int numberOfCoordinates = NUMBER_OF_AXIS_ * parent->getNOA();
  Eigen::VectorXd coordinates      = Eigen::VectorXd::Zero(numberOfCoordinates);
  unsigned int coordinatesIndex, atomIndex;
  atomIndex = 0;
  for (Atom *CurrAtom = HeadAtom; CurrAtom; CurrAtom = CurrAtom->getNext())
  {
    Eigen::Vector3d current = CurrAtom->Coordinates2Eigen(options);
    for (coordinatesIndex = 0; coordinatesIndex < NUMBER_OF_AXIS_;
         ++coordinatesIndex)
    {
      coordinates[atomIndex * NUMBER_OF_AXIS_ + coordinatesIndex] =
          current(coordinatesIndex);
    }
    ++atomIndex;
  }
  return coordinates;
}

double
Structure::calculateAllAtomRMSD(const Eigen::VectorXd referenceCoordinates,
                                const Eigen::VectorXd currentCoordinates)
{
  double rmsd = .0;
  for (int index = 0; index < referenceCoordinates.size(); ++index)
  {
    rmsd += pow((referenceCoordinates[index] - currentCoordinates[index]), 2.0);
  }
  rmsd = sqrt(rmsd / (referenceCoordinates.size() / 3.0));
  return rmsd;
}


int
initializeSphericalHarmonics(Molecule &CurrMol, Structure *CurrStruc)
{
  SphericalHarmonics *CurrY;
  SphericalHarmonicsMatrix *CurrYmat;
  CurrYmat       = CurrStruc->addHeadYmatrix();
  Atom *CurrAtom = CurrStruc->getHeadAtom();
  Bond *CurrBond;
  unsigned int i, j;
  while (CurrAtom)
  {
    CurrAtom->initializeHarmonics();
    CurrAtom = CurrAtom->getNext();
  }

  for (i = 0; i < CurrMol.getNOR(); ++i)
  {
    if (i == 0)
    {
      CurrY = CurrYmat->addHeadElement();
      CurrY->setRange(CurrY->getRDC()->getRange());
      continue;
    }
    CurrY = CurrY->appendHarmonics(CurrMol);
    CurrY->setRange(CurrY->getRDC()->getRange());
  }
  CurrY = CurrYmat->getHeadHarmonic();
  while (CurrY)
  {
    CurrY = CurrY->getNext();
  }
  CurrAtom = CurrMol.getHeadStruc()->getHeadAtom();
  while (CurrAtom)
  {
    for (i = 0; i < CurrAtom->getNumberOfBonds(); ++i)
    {
      CurrBond = CurrAtom->getBond(i);
      for (j = 0; j < CurrAtom->getNumberOfHarmonics(); ++j)
      {
        CurrY = CurrAtom->getBondharmonic(j);
        if (CurrY && (CurrY->getAtom1() == CurrBond->getBondpartner(CurrAtom) ||
                      CurrY->getAtom2() == CurrBond->getBondpartner(CurrAtom)))
        {
          CurrBond->setHarmonic(CurrY);
        }
      }
    }
    CurrAtom = CurrAtom->getNext();
  }
  CurrY    = NULL;
  CurrYmat = NULL;
  CurrAtom = NULL;
  CurrBond = NULL;
  delete CurrAtom;
  delete CurrYmat;
  delete CurrY;
  delete CurrBond;
  return 0;
}


double
Structure::get_radius_of_Gyration(StructureOptions opt)
{
  radiusOfGyration    = .0;
  Eigen::Vector3d COM = getCenterOfMass(opt);
  Eigen::Vector3d r;
  const double NOA = (double) parent->getNOA();
  for (Atom *CurrAtom = HeadAtom; CurrAtom; CurrAtom = CurrAtom->getNext())
  {
    r = CurrAtom->Coordinates2Eigen(opt) - COM;
    radiusOfGyration += (pow(r.norm(), 2.0));
  }
  radiusOfGyration /= NOA;
  return radiusOfGyration;
}

void
Yopt(double *p, double *x, int m, int n, void *data)
{
  if (n != m)
    std::cerr << "ERROR:\tStrange parameters in function Yopt "
                 "(double*,double*,int,int,void*)...\n\t\tContact your trusted "
                 "nerd...";
  struct Y_parameters *dptr;
  int j;
  dptr = (struct Y_parameters *) data;
  std::complex<double> tmp;
  tmp = 0.0;
  for (j = 0; j < HARMONIC_ELEMENTS_; ++j)
  {
    tmp += (D2Mm(j - 2, 0, p[0], p[1], 0.0) * dptr->Yref(dptr->rdc, j));
    //      tmp += ( sqrt(4.0*PI_/5.0)*std::conj(Ylm ( 2, j-2, p[1], p[0] )) *
    //      dptr->Yref(dptr->rdc,j) );
  }
  x[0] = 1.0 - tmp.real();
  x[1] = tmp.imag();

  dptr = NULL;
  delete dptr;
}

void
jacYopt(double *p, double *jac, int m, int n, void *data)
{
  if (n != m)
    std::cerr << "ERROR:\tStrange parameters in function Yopt "
                 "(double*,double*,int,int,void*)...\n\t\tContact your trusted "
                 "nerd...";
  struct Y_parameters *dptr;
  dptr = (struct Y_parameters *) data;
  std::complex<double> dd_dt, dd_dp;
  static const std::complex<double> im(.0, 1.0);
  double M;
  int j;
  dd_dt = dd_dp = 0.0;
  for (j = 0, M = -2.0; j < HARMONIC_ELEMENTS_; ++j, M += 1.0)
  {
    dd_dp -= (im * M * D2Mm(j - 2, 0, p[PHOBOS_PHI], p[PHOBOS_THETA], 0.0) *
              dptr->Yref(dptr->rdc, j));
    dd_dt += (std::exp(-M * im * p[PHOBOS_PHI]) *
              dd2Mmdb(j - 2, 0, p[PHOBOS_THETA]) * dptr->Yref(dptr->rdc, j));
  }
  jac[3] = jac[1] = -dd_dt.real();
  jac[0] = jac[2] = -dd_dp.real();

  dptr = NULL;
  delete dptr;
}

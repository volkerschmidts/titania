
#include <Atom.hpp>
#include <Declarations.hpp>
#include <Helper.hpp>
#include <LinAlg.hpp>
#include <Molecule.hpp>
#include <MonteCarloStatistic.hpp>
#include <RDCdata.hpp>
#include <RDCset.hpp>
#include <SphericalHarmonics.hpp>
#include <SphericalHarmonicsMatrix.hpp>
#include <Structure.hpp>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <phobos.h>


MonteCarloStep::MonteCarloStep(Molecule &CurrMol,
                               Structure *CurrStruc,
                               Flags &flags)
{
  /*
   * Constructor for the first initialized
   * Monte-Carlo bootstrap step.
   *
   * The bootstrap is heavy on computation
   * time or memory allocation. Since time
   * is the limiting factor we allocate
   * all matrizes needed, save the Monte-
   * Carlo values until we use the matrizes
   * and free the space after all is done
   * while TITANIA keeps on running the non
   * Monte-Carlo stuff.
   *
   * 100 TITANIA Iterations (performing 103
   * thousand MC-Steps) allocated 4 GB of
   * memory in total. With the cleanup
   * call its ~1 GB of memory for the same
   * calculation.
   */
  step = 0;
  NOS  = CurrMol.getNORsets();
  NOR  = CurrMol.getNOR();

  numericalDeviation = flags.numericalGradients;
  useMatrixFormalism = flags.calculateFullMatrix;

  RDCmatrix       = CurrMol.getRDCmatrix();
  RDCmatrixScaled = RDCmatrix;
  Dmax            = CurrStruc->getDmaxMatrix(StructureOptions::Optimized);
  DeltaRDC        = Eigen::MatrixXd::Zero(NOR, NOS);
  DeltaRDCScaled  = Eigen::MatrixXd::Zero(NOR, NOS);
  RDCVar          = Eigen::MatrixXd::Zero(NOR, NOS);
  RDCVarScaled    = Eigen::MatrixXd::Zero(NOR, NOS);

  A           = Eigen::MatrixXd::Zero(SAUPE_ELEMENTS_, NOS);
  S           = Eigen::MatrixXd::Zero(FULL_SAUPE_ELEMENTS_, NOS);
  EulerAngles = Eigen::MatrixXd::Zero(NUM_EULER_ANGLES_, NOS);
  EigenVals   = Eigen::MatrixXd::Zero(NUM_EIGENVALUES_, NOS);
  EigenVecs   = Eigen::MatrixXd::Zero(NUM_EIGENVECTORS_, NOS);

  p         = Eigen::MatrixXd::Zero(NUMBER_OF_SPHERICAL_COORDINATES, NOR);
  direction = Eigen::MatrixXd::Zero(NUMBER_OF_AXIS_, NOR);
  C         = CurrStruc->getCosineMatrix(StructureOptions::Optimized);
  if (useMatrixFormalism)
    weights = Eigen::MatrixXd::Identity(NOR, NOR);
  else
    weights = CurrMol.getWmatrix();

  next = NULL;
  prev = NULL;

  /*
   * Get a normal distributed matrix with
   * <x> = 0 and s(x) = 1
   *
   * Multiplying each element with the
   * corresponding delta(RDC) results in
   * a matrix with <x> = 0 and s(x) = d(RDC)
   *
   * Summation of each element with the
   * corresponding RDC results in a matrix
   * with <x> = RDC and s(x) = d(RDC)
   */

  randMat = MonteCarloStep::MC_randMatrix(NOR, NOS);
  SetupMatrizes(CurrStruc);
  VaryMatrix();
}

MonteCarloStep::MonteCarloStep(MonteCarloStep *previousStep)
{
  prev       = previousStep;
  prev->next = this;

  numericalDeviation = previousStep->numericalDeviation;
  useMatrixFormalism = previousStep->useMatrixFormalism;

  /*
   * MonteCarloStatistics initializes all steps
   * based on the HeadStep (preavious) using one
   * thread. This constructor is used to do so.
   *
   * Matrizes which are the same (RDCs, delta RDCs,
   * Dmax, etc.) are just copied. Everything else
   * is done like in standard constructor.
   */
  step = (prev->step + 1);
  NOS  = prev->NOS;
  NOR  = prev->NOR;

  RDCmatrix       = prev->RDCmatrix;
  RDCmatrixScaled = prev->RDCVarScaled;
  DeltaRDC        = prev->DeltaRDC;
  DeltaRDCScaled  = prev->DeltaRDCScaled;
  Dmax            = prev->Dmax;
  RDCVar          = Eigen::MatrixXd::Zero(NOR, NOS);
  RDCVarScaled    = Eigen::MatrixXd::Zero(NOR, NOS);

  A           = Eigen::MatrixXd::Zero(5, NOS);
  S           = Eigen::MatrixXd::Zero(6, NOS);
  EulerAngles = Eigen::MatrixXd::Zero(3, NOS);
  EigenVals   = Eigen::MatrixXd::Zero(3, NOS);
  EigenVecs   = Eigen::MatrixXd::Zero(9, NOS);

  p         = Eigen::MatrixXd::Zero(2, NOR);
  direction = Eigen::MatrixXd::Zero(NUMBER_OF_AXIS_, NOR);
  C         = prev->C;
  weights   = prev->weights;
  randMat   = MonteCarloStep::MC_randMatrix(NOR, NOS);

  next = NULL;
  VaryMatrix();
}

MonteCarloStep::~MonteCarloStep()
{
  prev = next = NULL;
}

void
MonteCarloStep::SetupMatrizes(Structure *CurrStruc)
{
  unsigned int i, j;
  RDCset *CurrSet;
  RDCdata *RDC;

  for (i = 0; i < NOR; ++i)
  {
    RDCmatrixScaled.row(i) = (RDCmatrix.row(i) / Dmax(i, i));
  }

  i = j = 0;
  for (CurrSet = CurrStruc->getParent()->getHeadSet(); (CurrSet && j < NOS);
       CurrSet = CurrSet->getNext(), ++j)
  {
    for (RDC = CurrSet->getHeadData(), i = 0; (RDC && i < NOR);
         RDC = RDC->getNext(), ++i)
    {
      DeltaRDC(i, j)       = RDC->getDeltaD();
      DeltaRDCScaled(i, j) = (DeltaRDC(i, j) / Dmax(i, i));
    }
  }
}

void
MonteCarloStep::VaryMatrix()
{
  unsigned int r, c;
  double one_over_Dmax = .0;
  for (r = 0; r < NOR; ++r)
  {
    one_over_Dmax = 1.0 / Dmax(r, r);
    for (c = 0; c < NOS; ++c)
    {
      RDCVar(r, c)       = RDCmatrix(r, c) + randMat(r, c) * DeltaRDC(r, c);
      RDCVarScaled(r, c) = RDCVar(r, c) * one_over_Dmax;
    }
  }
}

Eigen::MatrixXd
MonteCarloStep::MC_randMatrix(unsigned int &r, unsigned int &c)
{
  Eigen::MatrixXd randMat(r, c);
  unsigned int i, j;
  for (i = 0; i < r; ++i)
  {
    for (j = 0; j < c; ++j)
      randMat(i, j) = normRand(.0, 1.0);
  }
  return randMat;
}

MonteCarloStep *
MonteCarloStep::addNextStep()
{
  return new MonteCarloStep(this);
}

void
MonteCarloStep::run(BasicInformation &baseInformation,
                    Eigen::MatrixXd &InitialP)
{
  unsigned int i;
  /*
   * Solve the Saupe Eigen System which results
   * in the Saupe vectors of all rdc sets.
   */
  SaupeEigenSystems(C, weights, RDCVarScaled, A, EigenVals, EulerAngles,
                    EigenVecs, useMatrixFormalism);

  Eigen::MatrixXd S_tmp;

  /*
   * Transform the Saupe vectors in the
   * respective Saupe tensors.
   */

  for (i = 0; i < NOS; ++i)
  {
    S_tmp   = Sv2St(A.col(i));
    S(0, i) = S_tmp(STANDARD_AXIS_X_, STANDARD_AXIS_X_);
    S(1, i) = S_tmp(STANDARD_AXIS_Y_, STANDARD_AXIS_Y_);
    S(2, i) = S_tmp(STANDARD_AXIS_Z_, STANDARD_AXIS_Z_);
    S(3, i) = S_tmp(STANDARD_AXIS_X_, STANDARD_AXIS_Y_);
    S(4, i) = S_tmp(STANDARD_AXIS_X_, STANDARD_AXIS_Z_);
    S(5, i) = S_tmp(STANDARD_AXIS_Y_, STANDARD_AXIS_Z_);
  }

  RDC_Calc = Dcalc(Dmax, C, A);

  /*
   * Run a non object based TITANIA iteration.
   */
  F = S2F(EigenVals, EulerAngles);

  for (i = 0; i < NOS; ++i)
    RDCVarScaled.col(i) /= EigenVals(2, i);

  Y = back_calc_Y_ref(weights, RDCVarScaled, F, useMatrixFormalism);

  /*
   * getMonteCarloAngles performes the respective
   * Levenberg Marquardt optimizations of polar
   * angles.
   *
   * This function runs serial! Running this
   * in parallel would need a tracking of
   * the optimized and unoptimized angles!
   * Additionally the workspace would have to
   * be spread even more.
   */
  getMonteCarloAngles(InitialP, p, Y, NOR, baseInformation, numericalDeviation);

  /*
   * Encrypt transforms the optimized angles
   * in sin(theta), sin(phi), cos(theta) and
   * cos(phi). By this simple statistic can
   * be performed. For more information check
   * K.V. Mardia, Journal of the Royal
   * Statistics Society. Series B, Vol. 37,
   * No. 3 (1975), pp. 349-393:
   * Statistics of directional Data.
   */
  encrypt();
}

Eigen::MatrixXd
MonteCarloStep::getStatP()
{
  unsigned int r;
  for (r = 0; r < NOR; ++r)
  {
    while (p(1, r) < .0)
      p(1, r) += PI_TIMES_TWO_;
  }
  return p;
}

/*
 * All matrizes which are "trivial" and
 * just used to boost performance are
 * cleaned right after the Monte-Carlo
 * data were collected.
 */

void
MonteCarloStep::cleanup()
{
  MonteCarloStep *CurrS;
  for (CurrS = this; CurrS; CurrS = CurrS->getNext())
  {
    CurrS->RDCmatrix.resize(0, 0);
    CurrS->RDCmatrixScaled.resize(0, 0);
    CurrS->DeltaRDC.resize(0, 0);
    CurrS->DeltaRDCScaled.resize(0, 0);
    CurrS->RDCVar.resize(0, 0);
    CurrS->RDCVarScaled.resize(0, 0);
    CurrS->Dmax.resize(0, 0);

    CurrS->F_i.resize(0, 0);

    CurrS->C.resize(0, 0);
    CurrS->weights.resize(0, 0);
  }
}

MonteCarloStatistic::MonteCarloStatistic()
{
  mc_steps   = 0;
  stan_steps = MC_MINSTEPS_TITANIA;
  NOS        = 0;
  NOR        = 0;
  NOR2       = 0;
  MCSd       = 1.0;
  rmsd       = .0;
  HeadStep   = NULL;
}

MonteCarloStatistic::~MonteCarloStatistic()
{
  while (HeadStep && HeadStep->getNext())
  {
    HeadStep = HeadStep->getNext();
    delete HeadStep->getPrev();
  }
  delete HeadStep;
}

/***********************************************
 *                                             *
 *            How to run prallel               *
 *                                             *
 * Generate HeadStep                           *
 *   __________|                               *
 *   |                                         *
 *   V                                         *
 * pragma omp single                          *
 * ++|++++++++++++++++++++++++++++++++++++++++ *
 * + |_________________                      + *
 * + |                 |                     + *
 * + V                 |                     + *
 * + Add new Step -------                    + *
 * +                  |||___________________ + *
 * + TASK: Distribute ||| Data to Threads  | + *
 * + ###############+ ||| +##############  | + *
 * + #                 |                #  | + *
 * + #                 |                #  | + *
 * + #                 V                #  | + *
 * + #             run Calc             #  | + *
 * + #                 |___________________| + *
 * + #                 V                #  | + *
 * + #            save Values           #  | + *
 * + #                 |                #  | + *
 * + ################+ | +###############  | + *
 * +                   |____free Thread____| + *
 * +                                         + *
 * +++++++++++++++++++++++++++++++++++++++++++ *
 *                                             *
 ***********************************************/

int
MonteCarloStatistic::startMonteCarlo(Molecule &CurrMol,
                                     Structure *CurrStruc,
                                     BasicInformation &baseInformation,
                                     Flags &flags,
                                     MC_Run_Information &information)
{
  omp_init_lock(&A_lock);
  omp_init_lock(&A_sigm_lock);
  omp_init_lock(&S_lock);
  omp_init_lock(&S_sigm_lock);
  omp_init_lock(&D_lock);
  omp_init_lock(&D_sigm_lock);
  omp_init_lock(&D_calc_lock);
  omp_init_lock(&D_calc_sigm_lock);
  omp_init_lock(&p_lock);
  omp_init_lock(&p_sigm_lock);
  unsigned int i;

  if (information.console_output)
  {
#pragma omp critical
    *std::cin.tie()
        << "\nInformation:\tStarting Monte-Carlo bootstrapping...\n";
  }
  NOS  = CurrMol.getNORsets();
  NOR  = CurrMol.getNOR();
  NOR2 = NOR * 2;
  MCPi = Eigen::MatrixXd::Zero(NUMBER_OF_SPHERICAL_COORDINATES, NOR);

  SphericalHarmonics
      *CurrY; // = CurrStruc->getHeadYmatrix()->getHeadHarmonic();

  /*
   * Save the initial values of the optimized
   * polar coordinates. Use the standard
   * conventions ( theta=[0,pi], phi=[-pi,pi]).
   */
  for (i = 0, CurrY = CurrStruc->getHeadYmatrix()->getHeadHarmonic(); CurrY;
       CurrY = CurrY->getNext(), ++i)
  {
    MCPi(STANDARD_THETA_, i) = CurrY->getTheta(StructureOptions::Optimized);
    MCPi(STANDARD_PHI_, i)   = CurrY->getPhi(StructureOptions::Optimized);
    A2goodA(MCPi(STANDARD_THETA_, i), MCPi(STANDARD_PHI_, i));
  }

  A              = Eigen::MatrixXd::Zero(SAUPE_ELEMENTS_, NOS);
  S              = Eigen::MatrixXd::Zero(FULL_SAUPE_ELEMENTS_, NOS);
  D_calc         = Eigen::MatrixXd::Zero(NOR, NOS);
  D              = Eigen::MatrixXd::Zero(NOR, NOS);
  p              = Eigen::MatrixXd::Zero(NUMBER_OF_SPHERICAL_COORDINATES, NOR);
  direction_mean = Eigen::MatrixXd::Zero(NUMBER_OF_AXIS_, NOR);
  p_trig   = Eigen::MatrixXd::Zero(2 * NUMBER_OF_SPHERICAL_COORDINATES, NOR);
  p_sigm   = Eigen::MatrixXd::Zero(NUMBER_OF_SPHERICAL_COORDINATES, NOR);
  R_mean   = Eigen::VectorXd::Zero(NOR);
  R_2_mean = Eigen::VectorXd::Zero(NOR);
  circular_sigma = Eigen::VectorXd::Zero(NOR);

  /*
   * Add the first step. This is in omp single
   * section (no change to the thread usage
   * of standard TITANIA call).
   */
  HeadStep              = new MonteCarloStep(CurrMol, CurrStruc, flags);
  MonteCarloStep *CurrS = HeadStep;

  /*
   * Save the delta(RDC) values from HeadStep. This
   * values are needed to check the convergence
   * of the bootstrap.
   */

  dD = CurrS->getDeltaD();

  /*
   * Run the first Step and fetch the values.
   */
  CurrS->run(baseInformation, MCPi);
  fetchValues(CurrS);
  ++mc_steps;
  do
  {
    while (mc_steps < stan_steps)
    {
      /*
       * Master thread initializes all new
       * Monte-Carlo Steps and sends them as
       * task to the slaves.
       */
      CurrS = CurrS->addNextStep();
      ++mc_steps;
#pragma omp task firstprivate(CurrS, mc_steps) shared(MCPi, baseInformation)
      {
        /*
         * Slave steps run the Monte-Carlo
         * step and sends the result to a
         * sub slave who saves them.
         */
        CurrS->run(baseInformation, MCPi);
        if (baseInformation.numOfThreadsUsed > 1)
        {
#pragma omp task firstprivate(CurrS)
          fetchValues(CurrS);
        }
        else
          fetchValues(CurrS);
      }
    }
    /*
     * Master thread waits for all (sub)
     * slaves to be finished and run
     * checkConsistency.
     */
  } while (checkConsistency(information));
  /*(STANDARD_AXIS_X_, i) = p_trig(0,i)*p_trig(2,i); // sin(theta)*cos(phi)
   * Since checkConsistency does not run
   * in parallel we do not have to wait
   * for anything and reduceStatistics
   * can do its stuff.
   *
   * This function itself does run in parallel.
   * Since we have to destroy the locks in
   * the end we have to wait for this function
   * do go on. This means that reduceStatistics
   * has implicit barriers.
   */
  reduceStatistic();

  if (information.console_output)
    *std::cin.tie() << "\t\tMonte-Carlo converged in " << mc_steps
                    << " steps...\n";

  omp_destroy_lock(&A_lock);
  omp_destroy_lock(&A_sigm_lock);
  omp_destroy_lock(&S_lock);
  omp_destroy_lock(&S_sigm_lock);
  omp_destroy_lock(&D_lock);
  omp_destroy_lock(&D_sigm_lock);
  omp_destroy_lock(&D_calc_lock);
  omp_destroy_lock(&D_calc_sigm_lock);
  omp_destroy_lock(&p_lock);
  omp_destroy_lock(&p_sigm_lock);

  information.steps = mc_steps;
  return mc_steps;
}

void
MonteCarloStatistic::fetchValues(MonteCarloStep *CurrStep)
{
#pragma omp task
  {
    omp_set_lock(&A_lock);
    A = A + CurrStep->getA();
    omp_unset_lock(&A_lock);
  }
#pragma omp task
  {
    omp_set_lock(&S_lock);
    S = S + CurrStep->getS();
    omp_unset_lock(&S_lock);
  }
#pragma omp task
  {
    omp_set_lock(&D_calc_lock);
    D_calc = D_calc + CurrStep->getD_calc();
    omp_unset_lock(&D_calc_lock);
  }
#pragma omp task
  {
    omp_set_lock(&D_lock);
    D = D + CurrStep->getD();
    omp_unset_lock(&D_lock);
  }
#pragma omp task
  {
    omp_set_lock(&p_lock);
    direction_mean = direction_mean + CurrStep->getStat_direction();
    omp_unset_lock(&p_lock);
  }
}

/*************************************************************************************
 *                                                                                   *
 *                                  PARALLIZATION SCHEME *
 *                                                                                   *
 *  ++++++++++++++++ *
 *  +              + *
 *  +  HEADTHREAD  + *
 *  +      |       + *
 *  +++++# | #++++++ * | * STEP  | <-------- * |         | *
 *         |_________| *
 *        ||| *
 *        ||| *
 *        |||=============================================================== *
 *        |||                  |||                  |||                  ||| *
 *  ++++# ||| #+++++     ++++# ||| #+++++     ++++# ||| #+++++     ++++# |||
 *#+++++  *
 *  +     vvv      +     +     VVV      +     +     VVV      +     +     VVV + *
 *  +  HEADTHREAD  +     +  HEADTHREAD  +     +  HEADTHREAD  +     +  HEADTHREAD
 *+  *
 *  +              +     +              +     +              +     + +  *
 *  ++++++++++++++++     ++++++++++++++++     ++++++++++++++++ ++++++++++++++++
 **
 *                                                                                   *
 *************************************************************************************/

void
MonteCarloStatistic::reduceStatistic()
{
  A              = A * MCSd;
  S              = S * MCSd;
  D_calc         = D_calc * MCSd;
  direction_mean = direction_mean * MCSd;

  A_sigm      = Eigen::MatrixXd::Zero(SAUPE_ELEMENTS_, NOS);
  S_sigm      = Eigen::MatrixXd::Zero(FULL_SAUPE_ELEMENTS_, NOS);
  D_calc_sigm = Eigen::MatrixXd::Zero(NOR, NOS);
  direction_cov =
      Eigen::MatrixXd::Zero((NUMBER_OF_AXIS_ * NOR), (NUMBER_OF_AXIS_ * NOR));
  p_sigm = Eigen::MatrixXd::Zero(NUMBER_OF_SPHERICAL_COORDINATES, NOR);

  MonteCarloStep *CurrS = HeadStep;

/*
 * Master thread just starts with initializing
 * all matrizes and variables. Afterwards
 * the Masterthread distributes the HeadStep
 * to the slaves which run a full analysis
 * of the respective matrix they received.
 *
 * This scheme leeds to the use of fiewer
 * locks since only one thread uses one
 * matrix in one task!
 */
#pragma omp task firstprivate(CurrS)
  {
    while (CurrS)
    {
      addD_calc_sigm(CurrS);
      CurrS = CurrS->getNext();
    }
  }
#pragma omp task firstprivate(CurrS)
  {
    while (CurrS)
    {
      addA_sigm(CurrS);
      CurrS = CurrS->getNext();
    }
  }
#pragma omp task firstprivate(CurrS)
  {
    while (CurrS)
    {
      addp_sigm(CurrS);
      CurrS = CurrS->getNext();
    }
  }
#pragma omp task
  {
    for (unsigned int i = 0; i < NOR; ++i)
    {
      R_2_mean(i) =
          pow(direction_mean(STANDARD_AXIS_X_, i),
              2.0) + // sin(theta)*cos(phi) = x
          pow(direction_mean(STANDARD_AXIS_Y_, i),
              2.0) + // sin(theta)*sin(phi) = y
          pow(direction_mean(STANDARD_AXIS_Z_, i), 2.0); // cos(theta) = z
      R_mean(i)         = sqrt(R_2_mean(i));
      circular_sigma(i) = sqrt(-log2(R_2_mean(i)));
    }
  }
#pragma omp taskwait
  reduceSigmas();
  decrypt();
#pragma omp taskwait
}

void
MonteCarloStatistic::addA_sigm(MonteCarloStep *CurrS)
{
  Eigen::MatrixXd deltaA = A - CurrS->getA();
  Eigen::MatrixXd deltaS = S - CurrS->getS();
  unsigned int r, c;
  for (r = 0; r < 5; ++r)
  {
    for (c = 0; c < NOS; ++c)
    {
      deltaA(r, c) *= deltaA(r, c);
      deltaS(r, c) *= deltaS(r, c);
    }
  }
  for (c = 0; c < NOS; ++c)
    deltaS(5, c) *= deltaS(5, c);

  /*
   * The locks are not needed anymore since
   * only one thread is calculating the
   * sum of the squarded deviations.
   */
  A_sigm = A_sigm + deltaA;
  S_sigm = S_sigm + deltaS;
}

void
MonteCarloStatistic::addD_sigm(MonteCarloStep *CurrS)
{
  Eigen::MatrixXd delta = D - CurrS->getD();
  unsigned int r, c;
  for (r = 0; r < NOR; ++r)
  {
    for (c = 0; c < NOS; ++c)
      delta(r, c) *= delta(r, c);
  }
  D_sigm = D_sigm + delta;
}

void
MonteCarloStatistic::addD_calc_sigm(MonteCarloStep *CurrS)
{
  Eigen::MatrixXd delta = D_calc - CurrS->getD_calc();
  unsigned int r, c;
  for (r = 0; r < NOR; ++r)
  {
    for (c = 0; c < NOS; ++c)
      delta(r, c) *= delta(r, c);
  }
  D_calc_sigm = D_calc_sigm + delta;
}

void
MonteCarloStatistic::addp_sigm(MonteCarloStep *CurrS)
{
  // Diviation for current data point
  Eigen::MatrixXd delta = CurrS->getStat_direction() - direction_mean;

  // Covariance elements for current data point
  Eigen::MatrixXd cov_direction =
      Eigen::MatrixXd::Zero((NUMBER_OF_AXIS_ * NOR), (NUMBER_OF_AXIS_ * NOR));

  unsigned int BI; // Block matrix base indezes.
  unsigned int r, fx, fy;

  // Go for all rdc vectors (we only want the covariance within
  // the respective vectors). There should be no usefull covariance
  // between different vectors.

  for (r = 0; r < NOR; ++r)
  {
    BI = NUMBER_OF_AXIS_ * r;
    for (fx = 0; fx < NUMBER_OF_AXIS_; ++fx) // select axis 1
    {
      for (fy = 0; fy < NUMBER_OF_AXIS_; ++fy) // select axis 2
      {
        if (fy < fx)
        {
          cov_direction(BI + fx, BI + fy) = cov_direction(BI + fy, BI + fx);
          continue;
        }
        cov_direction(BI + fx, BI + fy) = delta(fx, r) * delta(fy, r);
      }
    }
  }
  direction_cov = direction_cov + cov_direction;
}

void
MonteCarloStatistic::reduceSigmas()
{
  A_sigm *= MCSd;
  S_sigm *= MCSd;
  D_calc_sigm *= MCSd;
  direction_cov *= MCSd;

  Eigen::MatrixXd weights = HeadStep->getWeights();
  unsigned int r, c;

  for (r = 0; r < NOR; ++r)
  {
    for (c = 0; c < NOS; ++c)
    {
      D_calc_sigm(r, c) = sqrt(D_calc_sigm(r, c));
      if (!HeadStep->useMatrizes() && weights(r, c) == .0)
        D_calc(r, c) = D_calc_sigm(r, c) = .0;
    }
  }

  for (r = 0; r < SAUPE_ELEMENTS_; ++r)
  {
    for (c = 0; c < NOS; ++c)
    {
      A_sigm(r, c) = sqrt(A_sigm(r, c));
      S_sigm(r, c) = sqrt(S_sigm(r, c));
    }
  }
  for (c = 0; c < NOS; ++c)
    S_sigm(5, c) = sqrt(S_sigm(5, c));

  // We need the covariance of sin(angle)/cos(angle) and not the standard
  // deviation like analogue The reason for this is, that coupling terms in
  // gauss propagation of errors need the covariance. Diagonal terms use the
  // standard deviation and square it -> rebuild the variance.
}

Eigen::VectorXd
d_angle_d_direction(unsigned int angle, double x, double y, double z)
{
  switch (angle)
  {
    case STANDARD_THETA_:
      return dtheta_dxi(x, y, z);
    case STANDARD_PHI_:
      return dphi_dxi(x, y, z);
  }
  return Eigen::VectorXd::Zero(NUMBER_OF_AXIS_);
}

void
MonteCarloStatistic::decrypt()
{
  Eigen::VectorXd gradatan2_ji; // Deviations of atan2. grad(0) = d(atan2)/d(x),
                                // grad(1) = d(atan2)/d(y)
  Eigen::VectorXd grad_da_dd;
  unsigned int BI; // base index of the block matrix containing the covariance
                   // of vector r in p_trig_cov
  unsigned int i, j, g1, g2;
  double cov, grad1, grad2; // current gradients and (co)variance
  p = Eigen::MatrixXd::Zero(NUMBER_OF_SPHERICAL_COORDINATES, NOR);
  for (i = 0; i < NOR; ++i)
  {
    BI = NUMBER_OF_AXIS_ * i;

    p(STANDARD_THETA_, i) =
        acos(direction_mean(STANDARD_AXIS_Z_, i) / R_mean(i));
    p(STANDARD_PHI_, i) =
        atan2(direction_mean(STANDARD_AXIS_Y_, i) / R_mean(i),
              direction_mean(STANDARD_AXIS_X_, i) / R_mean(i));
    for (int a = STANDARD_THETA_; a < NUMBER_OF_SPHERICAL_COORDINATES; ++a)
    {
      p_trig(a * 2, i)     = sin(p(a, i));
      p_trig(a * 2 + 1, i) = cos(p(a, i));
    }
    for (j = 0; j < NUMBER_OF_SPHERICAL_COORDINATES;
         ++j) // j -> current angle ( 0 = theta, 1 = phi )
    {
      grad_da_dd = d_angle_d_direction(
          j, direction_mean(STANDARD_AXIS_X_, i) / R_mean(i),
          direction_mean(STANDARD_AXIS_Y_, i) / R_mean(i),
          direction_mean(STANDARD_AXIS_Z_, i) / R_mean(i));

      for (g1 = 0; g1 < NUMBER_OF_AXIS_; ++g1)
      {
        grad1 = grad_da_dd(g1);
        for (g2 = 0; g2 < NUMBER_OF_AXIS_; ++g2)
        {
          grad2 = grad_da_dd(g2);
          cov   = direction_cov(BI + g1, BI + g2);
          p_sigm(j, i) += (cov * grad1 * grad2);
        }
      }
      p_sigm(j, i) = sqrt(p_sigm(j, i));
      continue;
    }
  }
}

int
MonteCarloStatistic::checkConsistency(MC_Run_Information &information)
{
  MCSd = (1.0 / ((double) mc_steps));
#pragma omp taskwait
  unsigned int r, c;
  Eigen::MatrixXd D_check = D * MCSd;
  D_sigm                  = Eigen::MatrixXd::Zero(NOR, NOS);
  rmsd                    = .0;

  MonteCarloStep *CurrS = HeadStep;

  while (CurrS)
  {
    {
      Eigen::MatrixXd delta = D_check - CurrS->getD();
      for (r = 0; r < NOR; ++r)
      {
        for (c = 0; c < NOS; ++c)
          delta(r, c) *= delta(r, c);
      }
      D_sigm = D_sigm + delta;
    }
    CurrS = CurrS->getNext();
  }
  D_sigm *= MCSd;

  for (c = 0; c < NOS; ++c)
  {
    for (r = 0; r < NOR; ++r)
    {
      D_sigm(r, c) = sqrt(D_sigm(r, c));
      rmsd += (pow((dD(r, c) - D_sigm(r, c)), 2.0));
    }
  }
  rmsd = rmsd / ((double) (NOR * NOS));
  rmsd = sqrt(rmsd);

  if (rmsd > information.convergence && stan_steps < information.max_steps)
  {
    stan_steps += 500;
    return 1;
  }
  else
  {
    D = D_check;
    return 0;
  }
}

MC_Output
MonteCarloStatistic::getStatistics()
{
  MC_Output o;
  o.p_mean      = p;
  o.p_trig_mean = p_trig;
  o.p_sigm      = p_sigm;

  o.R_2_mean       = R_2_mean;
  o.R_mean         = R_mean;
  o.circular_sigma = circular_sigma;

  o.D_mean = D;
  o.D_sigm = D_sigm;

  o.D_calc_mean = D_calc;
  o.D_calc_sigm = D_calc_sigm;

  o.Aligns      = A;
  o.Aligns_sigm = A_sigm;

  o.Saupe_tensor      = S;
  o.Saupe_tensor_sigm = S_sigm;

  o.steps = mc_steps;
  return o;
}

void
getMonteCarloAngles(Eigen::MatrixXd &Pi,
                    Eigen::MatrixXd &Po,
                    Eigen::MatrixXcd &MCY,
                    unsigned int NOR,
                    BasicInformation &baseInformation,
                    bool numericalGrad)
{
  double *p_Y, *x_Y, *work, *covar;

  Y_parameters Ydata;
  Ydata.Yref = MCY;

  int id = omp_get_thread_num();
  double info[10];
  unsigned int i;

  p_Y   = baseInformation.phobos_workspace_Sphericals[id];
  x_Y   = p_Y + NUMBER_OF_SPHERICAL_COORDINATES;
  work  = x_Y + NUMBER_OF_SPHERICAL_FUNCTIONS;
  covar = p_Y + baseInformation.phobos_worksize_Sphericals -
          PHOBOS_SPHERICAL_COVAR_SIZE;

  for (i = 0; i < NOR; ++i)
  {
    Ydata.rdc         = i;
    p_Y[PHOBOS_THETA] = Pi(STANDARD_THETA_, i);
    p_Y[PHOBOS_PHI]   = Pi(STANDARD_PHI_, i);
    x_Y[0] = x_Y[1] = .0;

    if (numericalGrad)
      phobos(Yopt, NULL, p_Y, x_Y, NUMBER_OF_SPHERICAL_COORDINATES,
             NUMBER_OF_SPHERICAL_FUNCTIONS,
             baseInformation.limits.max_lm_iterations, NULL, info, work, covar,
             (void *) &Ydata);
    else
      phobos(Yopt, jacYopt, p_Y, x_Y, NUMBER_OF_SPHERICAL_COORDINATES,
             NUMBER_OF_SPHERICAL_FUNCTIONS,
             baseInformation.limits.max_lm_iterations, NULL, info, work, covar,
             (void *) &Ydata);

    if (fabs(getPolar2Beta(p_Y[PHOBOS_THETA], p_Y[PHOBOS_PHI],
                           Pi(STANDARD_THETA_, i), Pi(STANDARD_PHI_, i))) >
        PI_HALF_)
    {
      p_Y[PHOBOS_THETA] += PI_;
    }
    A2goodA(p_Y[PHOBOS_THETA], p_Y[PHOBOS_PHI]);
    Po(STANDARD_THETA_, i) = p_Y[PHOBOS_THETA];
    Po(STANDARD_PHI_, i)   = p_Y[PHOBOS_PHI];
  }
  p_Y = x_Y = work = covar = NULL;
}

double
phiOutput(double p)
{
  if (p < .0)
  {
    return (360.0 + p);
  }
  else
    return p;
}

void
MonteCarloStep::encrypt()
{
  p_trig = Eigen::MatrixXd::Zero(NUMBER_OF_TRIG_COORDINATES_, NOR);
  unsigned int i, j;
  for (j = 0; j < NUMBER_OF_SPHERICAL_COORDINATES; ++j)
  {
    for (i = 0; i < NOR; ++i)
      p_trig((2 * j), i) = sin(p(j, i));
    for (i = 0; i < NOR; ++i)
      p_trig((2 * j + 1), i) = cos(p(j, i));
  }
  for (i = 0; i < NOR; ++i)
  {
    direction(STANDARD_AXIS_X_, i) =
        p_trig(0, i) * p_trig(3, i); // sin(theta)*cos(phi)
    direction(STANDARD_AXIS_Y_, i) =
        p_trig(0, i) * p_trig(2, i);               // sin(theta)*sin(phi)
    direction(STANDARD_AXIS_Z_, i) = p_trig(1, i); // cos(theta)
  }
}

Eigen::MatrixXd
estimateSaupeErrors(Eigen::MatrixXd A_Vec,
                    Eigen::MatrixXd &Saupe,
                    Eigen::MatrixXd &Evec,
                    Eigen::MatrixXd &Eval,
                    Eigen::Vector3d &Euler,
                    Eigen::MatrixXd A_Vec_sigm,
                    // Eigen::MatrixXd &Saupe_sigm, // This value is computed
                    // (not estimated) in MC part
                    Eigen::MatrixXd &Evec_sigm,
                    Eigen::MatrixXd &Eval_sigm,
                    Eigen::Vector3d &Euler_sigm)
{
  double fac = 1e-4;
  double h, g, sigm_Si;

  int i, m, n;

  Eigen::MatrixXd Evec_tmp   = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd S_MC_Gauss = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd Saupe_tmp  = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd Eval_tmp   = Eigen::MatrixXd::Zero(3, 1);
  Eigen::Vector3d Euler_tmp  = Eigen::Vector3d::Zero(3);

  Evec = Evec_sigm = Evec_tmp;
  Eval = Eval_sigm = Eval_tmp;
  Euler = Euler_sigm = Euler_tmp;

  Euler = Tensor2Euler(Saupe, Eval, Evec, 0, true, true);

  for (i = 0; i < 5; ++i)
  {
    h = A_Vec(i, 0) * fac;
    A_Vec(i, 0) += h;
    Saupe_tmp = Sv2St(A_Vec);
    A_Vec(i, 0) -= h;
    Euler_tmp = Tensor2Euler(Saupe_tmp, Eval_tmp, Evec_tmp, 0, true, true);
    Eval_tmp  = Eval_tmp - Eval;
    Evec_tmp  = Evec_tmp - Evec;
    Euler_tmp = Euler_tmp - Euler;
    Saupe_tmp = Saupe_tmp - Saupe;
    sigm_Si   = A_Vec_sigm(i, 0);
    for (n = 0; n < 3; ++n)
    {
      g = Eval_tmp(n) / h;
      Eval_sigm(n) += (pow((g * sigm_Si), 2.0));
      g = Euler_tmp(n) / h;
      Euler_sigm(n) += (pow((g * sigm_Si), 2.0));
      for (m = 0; m < 3; ++m)
      {
        g = Evec_tmp(n, m) / h;
        Evec_sigm(n, m) += (pow((g * sigm_Si), 2.0));
        g = Saupe_tmp(n, m) / h;
        S_MC_Gauss(n, m) += (pow((g * sigm_Si), 2.0));
      }
    }
  }
  for (m = 0; m < 3; ++m)
  {
    Eval_sigm(m)  = sqrt(Eval_sigm(m));
    Euler_sigm(m) = sqrt(Euler_sigm(m));
    for (n = 0; n < 3; ++n)
    {
      Evec_sigm(m, n)  = sqrt(Evec_sigm(m, n));
      S_MC_Gauss(m, n) = sqrt(S_MC_Gauss(m, n));
    }
  }
  return S_MC_Gauss;
}

void
MonteCarloStatistic::cleanup()
{
  HeadStep->cleanup();
}

void
MonteCarloStatistic::plot_mc_angles(Molecule &CurrMol,
                                    Structure *CurrStruc,
                                    BasicInformation &baseInformation,
                                    Flags &flags)
{
  unsigned int NOS = CurrMol.getNORsets();
  unsigned int i, j, set;
  Eigen::MatrixXd eigenVecs;
  RDCset *CurrSet, *HeadSet;
  HeadSet = CurrMol.getHeadSet();

  set = 0;
  i   = 0;

  // Maybe one could save the file-pointer somewhere or
  // even keep them open after ali-files were generated.
  // But fastest implementation is to reopen them with the
  // "a" option.
  FILE **outputs = (FILE **) malloc(NOS * sizeof(FILE *));
  char *buf      = (char *) malloc(STANDARD_BUFFER_SIZE * sizeof(char));

  for (CurrSet = HeadSet; CurrSet; CurrSet = CurrSet->getNext())
  {
    snprintf(buf, STANDARD_BUFFER_SIZE * sizeof(char), "%s.%s.ali",
             baseInformation.outputFileName.c_str(),
             CurrSet->getLabel().c_str());
    outputs[set++] = fopen(buf, "a");
  }

  for (MonteCarloStep *currStep = HeadStep; currStep;
       currStep                 = currStep->getNext())
  {
    eigenVecs = currStep->getEigenVecs();
    for (set = 0; set < NOS; ++set)
    {
      for (j = 0; j < NUM_EIGENVECTORS_; ++j)
      {
        fprintf(outputs[set], " %7.4f", eigenVecs(j, set));
      }
      fprintf(outputs[set], "\n");
    }
  }

  for (i = 0; i < NOS; ++i)
  {
    fclose(outputs[i]);
  }

  free(outputs);
  free(buf);
  HeadSet = NULL;
  CurrSet = NULL;

  flags.plotMonteCarlo = false;
  return;
// TODO update plot_mc_angles
#ifdef DEACTIVATE
  if (baseInformation.MC_stop_reason == MC_stop::MaxIter)
  {
    flags.plotMonteCarlo = false;
    fprintf(baseInformation.output, "Since the stop reason for TITANIA was max "
                                    "iterations Monte Carlo is skiped! \n");
    fprintf(
        baseInformation.output,
        "Rerun TITANIA with more iterations allowed or check your data! \n");
    return;
  }
  Eigen::MatrixXd out;
  baseInformation.output << "Monte Carlo output:\n";
  unsigned int i;
  double theta, phi;
  for (SphericalHarmonics *CurrY =
           CurrStruc->getHeadYmatrix()->getHeadHarmonic();
       CurrY; CurrY = CurrY->getNext())
  {
    baseInformation.output << std::left << std::setw(9)
                           << (CurrY->getAtom1()->getIdentifier() +
                               CurrY->getAtom2()->getIdentifier())
                           << " ";
    baseInformation.output << std::left << std::setw(9)
                           << (CurrY->getAtom1()->getIdentifier() +
                               CurrY->getAtom2()->getIdentifier())
                           << " ";
  }
  baseInformation.output << std::endl;
  for (MonteCarloStep *CurrS = HeadStep; CurrS; CurrS = CurrS->getNext())
  {
    out = CurrS->getStatP();
    for (i = 0; i < NOR; ++i)
    {
      theta = out(0, i);
      phi   = out(1, i);
      A2goodA(theta, phi);
      baseInformation.output << std::left << std::setw(10)
                             << std::setprecision(4) << theta << std::left
                             << std::setw(10) << std::setprecision(4) << phi;
    }
    baseInformation.output << std::endl;
  }
#endif
}

void
MonteCarloStatistic::print_mc_polar_angles(Structure *CurrStruc,
                                           BasicInformation &baseInformation)
{
  Eigen::MatrixXd out;
  unsigned int i;
  double theta, phi;

  fprintf(baseInformation.output, "Monte Carlo output:\n");

  for (SphericalHarmonics *CurrY =
           CurrStruc->getHeadYmatrix()->getHeadHarmonic();
       CurrY; CurrY = CurrY->getNext())
  {
    fprintf(baseInformation.output, "t(%9s) ",
            (CurrY->getAtom1()->getIdentifier() +
             CurrY->getAtom2()->getIdentifier())
                .c_str());
    fprintf(baseInformation.output, "p(%9s) ",
            (CurrY->getAtom1()->getIdentifier() +
             CurrY->getAtom2()->getIdentifier())
                .c_str());
  }
  fprintf(baseInformation.output, "\n");
  for (MonteCarloStep *CurrS = HeadStep; CurrS; CurrS = CurrS->getNext())
  {
    out = CurrS->getStatP();
    for (i = 0; i < NOR; ++i)
    {
      theta = out(0, i);
      phi   = out(1, i);
      A2goodA(theta, phi);
      fprintf(baseInformation.output, "%9.4f ", theta);
      fprintf(baseInformation.output, "%9.4f ", phi);
    }
    fprintf(baseInformation.output, "\n");
  }
}

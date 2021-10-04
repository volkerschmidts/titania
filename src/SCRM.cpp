
#include <Atom.hpp>
#include <Declarations.hpp>
#include <Helper.hpp>
#include <LinAlg.hpp>
#include <Molecule.hpp>
#include <RDCdata.hpp>
#include <RDCset.hpp>
#include <SCRM.hpp>
#include <SphericalHarmonics.hpp>
#include <SphericalHarmonicsMatrix.hpp>
#include <Structure.hpp>
#include <omp.h>
#include <phobos.h>

int
computeAngles(Molecule &CurrMol,
              Structure *CurrStruc,
              BasicInformation &baseInformation,
              Flags &flags,
              Y_parameters &Ydata)
{
  baseInformation.state =
      computeYref(CurrMol, CurrStruc, baseInformation, flags, Ydata,
                  StructureOptions::Initial, false);
  if (baseInformation.state)
    return baseInformation.state;

  baseInformation.state = SCRM(CurrStruc, baseInformation, flags, Ydata, false);
  if (baseInformation.state || !flags.scaleWithSoverall)
    return baseInformation.state;

  baseInformation.state =
      computeYref(CurrMol, CurrStruc, baseInformation, flags, Ydata,
                  StructureOptions::Initial, true);
  if (baseInformation.state)
    return baseInformation.state;

  baseInformation.state = SCRM(CurrStruc, baseInformation, flags, Ydata, true);
  return baseInformation.state;
}

/*
 * Performs a MFA step according to N.-A. Lakomek, K. F. A. Walter, C. Farès, O.
 * F. Lange, B. L. de Groot, H. Grubmüller, R. Brüschweiler, A. Munk, S. Becker,
 * J. Meiler, C. Griesinger, Journal of Biomolecular NMR 2008, 41, 139,
 * DOI 10.1007/s10858-008-9244-4. In contrast to the actual SCRM implementation
 * just one MFA step is performed.
 */

int
SCRM(Structure *CurrStruc,
     BasicInformation &baseInformation,
     Flags &flags,
     Y_parameters &Ydata,
     bool Sscaling)
{
  SphericalHarmonics *CurrY;
  StructureOptions opt;
  if (Sscaling)
    opt = StructureOptions::Optimized;
  else
    opt = StructureOptions::Initial;

  double *par_Y, *mea_Y, *work, *covar;
  int rdccounter = 0, i, id;

  /*
   * Check if phobos workspace was allready allocated.
   * This should be done as late as possible and in a parallel
   * section, since this will lead to a large displacement
   * of memory locations! If the memory locacations are
   * connected one might run into false sharing problems.
   */
  if (baseInformation.phobos_workspace_Sphericals[0] == NULL)
  {
    for (i = 0; i < baseInformation.numOfThreadsUsed; ++i)
    {
#pragma omp task firstprivate(i) shared(baseInformation)
      baseInformation.phobos_workspace_Sphericals[i] = (double *) malloc(
          baseInformation.phobos_worksize_Sphericals * sizeof(double));
    }
/*
 * Make sure that every task has allocated the respective
 * memory before starting the optimization.
 */
#pragma omp taskwait
  }

  for (CurrY = CurrStruc->getHeadYmatrix()->getHeadHarmonic();
       CurrY && !Sscaling; CurrY = CurrY->getNext())
  {
    Ydata.rdc = rdccounter++;
/*
 * Every task needs a copy of the variables. Additionally
 * a copy of the Ydata and current spherical harmonics
 * are needed. BasicInformation and Flags can be shared.
 */
#pragma omp task private(par_Y, mea_Y, work, covar, id, i) \
    firstprivate(Ydata, CurrY) shared(baseInformation, flags)
    {
      // Check out the thread id
      id = omp_get_thread_num();

      // Clean the memory of for the thread
      for (i = 0; i < baseInformation.phobos_worksize_Sphericals; ++i)
        baseInformation.phobos_workspace_Sphericals[id][i] = .0;

      // Devide the memory in the small pieces.
      par_Y = baseInformation.phobos_workspace_Sphericals[id];
      mea_Y = par_Y + NUMBER_OF_SPHERICAL_COORDINATES;
      work  = mea_Y + NUMBER_OF_SPHERICAL_FUNCTIONS;
      covar = par_Y + baseInformation.phobos_worksize_Sphericals -
              PHOBOS_SPHERICAL_COVAR_SIZE;

      // Define the starting point.
      par_Y[PHOBOS_THETA] = CurrY->getTheta(opt);
      par_Y[PHOBOS_PHI]   = CurrY->getPhi(opt);

      // Run the phobos optimization on numerical or analytical gradients
      if (flags.numericalGradients)
      {
        phobos(Yopt, NULL, par_Y, mea_Y, NUMBER_OF_SPHERICAL_COORDINATES,
               NUMBER_OF_SPHERICAL_FUNCTIONS,
               baseInformation.limits.max_lm_iterations,
               baseInformation.limits.phobos_opts_sphericals,
               CurrY->getLMinfo(), work, covar, (void *) &Ydata);
      }
      else
      {
        phobos(Yopt, jacYopt, par_Y, mea_Y, NUMBER_OF_SPHERICAL_COORDINATES,
               NUMBER_OF_SPHERICAL_FUNCTIONS,
               baseInformation.limits.max_lm_iterations,
               baseInformation.limits.phobos_opts_sphericals,
               CurrY->getLMinfo(), work, covar, (void *) &Ydata);
      }

      /*
       *  Save the optimized angles and analyze the
       *  motional information on the RDCs based
       *  on the refined spherical harmonics and
       *  optimized spherical angles.
       */
      CurrY->setOptimizedAngles(par_Y[PHOBOS_THETA], par_Y[PHOBOS_PHI]);

      par_Y = mea_Y = work = covar = NULL;
    }
  }
  CurrY = NULL;

// Make sure that all motional information were collected.
#pragma omp taskwait

  for (CurrY = CurrStruc->getHeadYmatrix()->getHeadHarmonic(), rdccounter = 0;
       CurrY; CurrY = CurrY->getNext(), ++rdccounter)
    analyzeYtensor(CurrY, Ydata.Yref.row(rdccounter));
  if (Sscaling)
    return 0;

  // Determine S(overall) for the structure.
  SrdcBySoverall(CurrStruc->getHeadYmatrix());
  return 0;
}

int
computeYref(Molecule &CurrMol,
            Structure *CurrStruc,
            BasicInformation &baseInformation,
            Flags &flags,
            Y_parameters &Ydata,
            StructureOptions opt,
            bool Sscaling)
{
  unsigned int NOS = CurrMol.getNORsets();
  unsigned int i;
  Eigen::MatrixXcd Fmatrix, Fmatrix_i;
  Eigen::MatrixXd SaupeEigenValues, rdcMatrix, w;

  SphericalHarmonicsMatrix *Ymatrix = CurrStruc->getHeadYmatrix();
  Eigen::MatrixXd weights           = CurrMol.getWmatrix();

  /*
   * Make sure the latest (normalized) cosine matrix
   * is used for the determination of the orientation.
   */

  Eigen::MatrixXd C = Ymatrix->determineBmatrix(CurrMol, CurrStruc, opt);

  /*
   * Collect the orientational information on the
   * current structure and decompose the resulting
   * order tensors via eigenvalue decomposition.
   * Build the F-matrix (Wigner rotations) based
   * on the result and build its pseudo inverse
   * (= Moore-Penrose inverse).
   */
  Ymatrix->determineEuler(CurrMol, CurrStruc, baseInformation, flags, opt);

  Fmatrix = Ymatrix->determineFmatrix();

  SaupeEigenValues = Ymatrix->getSaupeEigenValues();

  /*
   * Build the latest scaled RDC matrix. Rows
   * (1 spin pair in all media) are scaled with
   * Dmax and columns (all spin pairs of 1 medium)
   *  are scaled with the S(zz) eigenvalue
   */
  rdcMatrix = CurrStruc->getRDCmatrix(rdcMatrixOptions::Scaled);

  for (i = 0; i < NOS; ++i)
    rdcMatrix.col(i) /= (SaupeEigenValues(2, i));
  if (Sscaling)
    rdcMatrix *= Ymatrix->getSoverall();

  /*
   * Calculate the refined spherical harmonics
   * <Y(ref)>.
   * IMPORTANT: CurrMol.getWmatrix() returns
   * the unit matrix. In the beginning of the
   * programming an weighting was used which
   * is now turned of.
   */
  Ydata.Yref =
      back_calc_Y_ref(weights, rdcMatrix, Fmatrix, flags.calculateFullMatrix);
  Ymatrix->setYref(Ydata.Yref);

  Ymatrix = NULL;
  return 0;
}

void
analyzeYtensor(SphericalHarmonics *CurrY, Eigen::MatrixXcd Y)
{
  double Srdc, Sax, eta, phi_aniso;

  /*
   * Calculate the axial component S^2(rdc),
   * rhombicity eta and orientation of the
   * anisotropic motion phi_aniso. Addionally
   * approximated value for S^2 is calculated.
   * This value uses the assumption that no
   * anisotropic motion is pressent and thus
   * the cosine tensors are symmetrical.
   * All values are calculated in the RDC
   * vector frame!
   */
  motionalAnalysisOfSphericalHarmonics(
      Y, CurrY->getTheta(StructureOptions::Optimized),
      CurrY->getPhi(StructureOptions::Optimized), Srdc, Sax, eta, phi_aniso);

  CurrY->setS2rdc(Srdc);
  CurrY->setSaxial(Sax);
  CurrY->setetardc(eta);
  CurrY->setAnisoTheta(phi_aniso);
}

void
SrdcBySoverall(SphericalHarmonicsMatrix *YM)
{
  /*
   * The S^2 parameter is in average equal to one!
   * Neglecting the rare case where S^2_i = 1 for
   * all i we thus have some S^2_i > 1.
   *
   * To counter act this problem all S^2(rdc) are
   * scaled to the largest value. For more information
   * check Peti et. al., JACS 2002.
   * DOI:10.1021/ja011883c
   *
   * CAUTION: S(overall) = sqrt( 1 / S^2(rdc,max) )
   *
   * S^2(rdc,scaled) = S(overall)^2 * S^2
   *
   * For reasons of conventions S(overall) is not
   * squared!
   */

  SphericalHarmonics *CurrY;
  double S = .0, maxSrdc = 1.0;

  for (CurrY = YM->getHeadHarmonic(); CurrY; CurrY = CurrY->getNext())
  {
    if (CurrY->getS2rdc() > maxSrdc)
      maxSrdc = CurrY->getS2rdc();
  }

  maxSrdc = 1.0 / maxSrdc;
  YM->setSoverall(sqrt(maxSrdc));

  for (CurrY = YM->getHeadHarmonic(); CurrY; CurrY = CurrY->getNext())
  {
    S = CurrY->getS2rdc() * maxSrdc;
    CurrY->setS2rdc(S);
    S = CurrY->getSaxial() * maxSrdc;
    CurrY->setSaxial(S);
  }
}

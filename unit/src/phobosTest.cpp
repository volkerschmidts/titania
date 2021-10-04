
#include <main.h>
#include <math.h>
#include <phobos.h>
#include <phobosTest.h>
// there is still some weird include ordering dependence hidden in
// Declarations.hpp
// clang-format off
#include <Declarations.hpp>
// clang-format on

/******************************************
 *                                        *
 * Define some variables used for testing *
 *        phobos on some functions.       *
 *                                        *
 ******************************************/

#ifndef PHOBOS_ROSENBROCK_PARAMETERS
#define PHOBOS_ROSENBROCK_PARAMETERS 2
#define PHOBOS_ROSENBROCK_A          1
#define PHOBOS_ROSENBROCK_B          283.51
#define PHOBOS_ROSENBROCK_NPAR       2
#define PHOBOS_ROSENBROCK_NMEA       2
#define PHOBOS_ROSENBROCK_COVAR      4
#define PHOBOS_ROSENBROCK_X_START    -1.2
#define PHOBOS_ROSENBROCK_Y_START    1.0
#define PHOBOS_ROSENBROCK_X_END      1.0
#define PHOBOS_ROSENBROCK_Y_END      1.0
#endif

#ifndef PHOBOS_RASTRIGIN_PARAMETERS
#define PHOBOS_RASTRIGIN_PARAMETERS 1
#define PHOBOS_RASTRIGIN_A          10
#define PHOBOS_RASTRIGIN_NPAR       2
#define PHOBOS_RASTRIGIN_NMEA       2
#define PHOBOS_RASTRIGIN_COVAR      4
#define PHOBOS_RASTRIGIN_X_START    0.4
#define PHOBOS_RASTRIGIN_Y_START    -0.4
#define PHOBOS_RASTRIGIN_X_END      0.0
#define PHOBOS_RASTRIGIN_Y_END      0.0
#endif

#ifndef PHOBOS_BEALE_PARAMETERS
#define PHOBOS_BEALE_PARAMETERS 3
#define PHOBOS_BEALE_A          1.5
#define PHOBOS_BEALE_B          2.25
#define PHOBOS_BEALE_C          2.625
#define PHOBOS_BEALE_NPAR       2
#define PHOBOS_BEALE_NMEA       2
#define PHOBOS_BEALE_COVAR      4
#define PHOBOS_BEALE_X_START    10.0
#define PHOBOS_BEALE_Y_START    10.0
#define PHOBOS_BEALE_X_END      3.0
#define PHOBOS_BEALE_Y_END      0.5
#endif

void
Rosenbrock_4_phobos(double *p, double *x, int m, int n, void *data)
{
  int i;
  double *Ros_Par = (double *) data;
  for (i = 0; i < n; ++i)
  {
    x[i] = (Ros_Par[0] - p[0]) * (Ros_Par[0] - p[0]) +
           Ros_Par[1] * (p[1] - p[0] * p[0]) * (p[1] - p[0] * p[0]);
  }
  //   std::cout << "f(" << p[0] << "," << p[1] << ") = " << x[0] << std::endl;
  Ros_Par = NULL;
}

void
jac_Rosenbrock_4_phobos(double *p, double *jac, int m, int n, void *data)
{
  int i;
  double *Ros_Par = (double *) data;
  for (i = 0; i < n * m;)
  {
    jac[i++] = (-2.0 * (Ros_Par[0] - p[0]) -
                4.0 * p[0] * Ros_Par[1] * (p[1] - p[0] * p[0]));
    jac[i++] = (2.0 * Ros_Par[1] * (p[1] - p[0] * p[0]));
  }

  Ros_Par = NULL;
}

void
Rastrigin_4_phobos(double *p, double *x, int m, int n, void *data)
{
  int j, i;
  double A = *((double *) data);
  for (j = 0; j < n; ++j)
  {
    x[j] = A * 2.0;
    for (i = 0; i < m; ++i)
    {
      x[j] += (p[i] * p[i] - A * cos(2.0 * PI_ * p[i]));
    }
  }
  //   std::cout << "f(" << p[0] << "," << p[1] << ") = " << x[0] << std::endl;
}

void
Beale_4_phobos(double *p, double *x, int m, int n, void *data)
{
  int i;
  double *Bea_Par = (double *) data;
  for (i = 0; i < n; ++i)
  {
    x[i] =
        ((Bea_Par[0] - p[0] + p[0] * p[1]) * (Bea_Par[0] - p[0] + p[0] * p[1]) +
         (Bea_Par[1] - p[0] + p[0] * p[1] * p[1]) *
             (Bea_Par[1] - p[0] + p[0] * p[1] * p[1]) +
         (Bea_Par[2] - p[0] + p[0] * p[1] * p[1] * p[1]) *
             (Bea_Par[2] - p[0] + p[0] * p[1] * p[1] * p[1]));
  }
  //   std::cout << "f(" << p[0] << "," << p[1] << ") = " << x[0] << std::endl;
  Bea_Par = NULL;
}

void
jac_Beale_4_phobos(double *p, double *jac, int m, int n, void *data)
{
  double *Bea_Par = (double *) data;
  int i;
  for (i = 0; i < m * n;)
  {
    jac[i++] = ( // df/dx
        2.0 * (p[1] - 1.0) * (Bea_Par[0] - p[0] + p[0] * p[1]) +
        2.0 * (p[1] * p[1] - 1.0) * (Bea_Par[1] - p[0] + p[0] * p[1] * p[1]) +
        2.0 * (p[1] * p[1] * p[1] - 1.0) *
            (Bea_Par[2] - p[0] + p[0] * p[1] * p[1] * p[1]));
    jac[i++] = ( // df/dy
        2.0 * (p[0]) * (Bea_Par[0] - p[0] + p[0] * p[1]) +
        2.0 * (2.0 * p[1] * p[0]) * (Bea_Par[1] - p[0] + p[0] * p[1] * p[1]) +
        2.0 * (3.0 * p[1] * p[1] * p[0]) *
            (Bea_Par[2] - p[0] + p[0] * p[1] * p[1] * p[1]));
  }
  Bea_Par = NULL;
}

CPPUNIT_TEST_SUITE_REGISTRATION(phobosTest);

phobosTest::phobosTest()
{
  ;
}

void
phobosTest::testPhobosRosenbrock()
{
  double *par, *mea, *Ros_Par, *work, *covar;
  int npar, nmea, worksize;

  // Calculate worksize
  npar     = PHOBOS_ROSENBROCK_NPAR;
  nmea     = PHOBOS_ROSENBROCK_NMEA;
  worksize = PHOBOS_WORKSIZE(npar, nmea);

  // Allocate and clean the workspace
  double *phobos_workspace = (double *) malloc(worksize * sizeof(double));
  for (int i = 0; i < worksize; ++i)
    phobos_workspace[i] = .0;
  double *opts = (double *) malloc(4 * sizeof(double));
  double *info = (double *) malloc(PHOBOS_INFOSIZE * sizeof(double));
  Ros_Par = (double *) malloc(PHOBOS_ROSENBROCK_PARAMETERS * sizeof(double));

  // Choose values to force LM to the optimum.
  opts[0] = 1e-10;
  opts[1] = opts[2] = opts[3] = 7e-27;

  // Distribute the workspace
  par   = phobos_workspace;
  mea   = par + npar;
  work  = mea + nmea;
  covar = par + worksize - PHOBOS_ROSENBROCK_COVAR;

  // Fill the needed function values
  Ros_Par[0] = PHOBOS_ROSENBROCK_A;
  Ros_Par[1] = PHOBOS_ROSENBROCK_B;
  par[0]     = PHOBOS_ROSENBROCK_X_START;
  par[1]     = PHOBOS_ROSENBROCK_Y_START;

  // Run and evaluate LM
  phobos(Rosenbrock_4_phobos, jac_Rosenbrock_4_phobos, par, mea, npar, nmea,
         10000, opts, info, work, covar, (void *) Ros_Par);
  CPPUNIT_ASSERT(fabs(par[0] - PHOBOS_ROSENBROCK_X_END) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(par[1] - PHOBOS_ROSENBROCK_Y_END) < ALLOWED_ERROR);

  // Clean the workspace
  par = mea = work = covar = NULL;
  free(Ros_Par);
  free(phobos_workspace);
  free(opts);
  free(info);
}

void
phobosTest::testPhobosRastriginNumerical()
{
  double *par, *mea, *Ras_Par, *work, *covar;
  int npar, nmea, worksize;

  // Calculate worksize
  npar     = PHOBOS_RASTRIGIN_NPAR;
  nmea     = PHOBOS_RASTRIGIN_NMEA;
  worksize = PHOBOS_WORKSIZE(npar, nmea);

  // Allocate and clean the workspace
  double *phobos_workspace = (double *) malloc(worksize * sizeof(double));
  for (int i = 0; i < worksize; ++i)
    phobos_workspace[i] = .0;
  double *opts = (double *) malloc(4 * sizeof(double));
  double *info = (double *) malloc(PHOBOS_INFOSIZE * sizeof(double));
  Ras_Par = (double *) malloc(PHOBOS_RASTRIGIN_PARAMETERS * sizeof(double));

  // Choose values to force LM to the optimum.
  opts[0] = 1e-10;
  opts[1] = opts[2] = opts[3] = 7e-27;

  // Distribute the workspace
  par   = phobos_workspace;
  mea   = par + npar;
  work  = mea + nmea;
  covar = par + worksize - PHOBOS_RASTRIGIN_COVAR;

  // Fill the needed function values
  Ras_Par[0] = PHOBOS_RASTRIGIN_A;

  par[0] = PHOBOS_RASTRIGIN_X_START;
  par[1] = PHOBOS_RASTRIGIN_Y_START;

  // Run and evaluate LM
  phobos(Rastrigin_4_phobos, NULL, par, mea, npar, nmea, 10000, opts, info,
         work, covar, (void *) Ras_Par);
  CPPUNIT_ASSERT(fabs(par[0] - PHOBOS_RASTRIGIN_X_END) <
                 ALLOWED_NUMERICAL_ERROR);
  CPPUNIT_ASSERT(fabs(par[1] - PHOBOS_RASTRIGIN_Y_END) <
                 ALLOWED_NUMERICAL_ERROR);

  // Clean the workspace
  par = mea = work = covar = NULL;
  free(Ras_Par);
  free(phobos_workspace);
  free(opts);
  free(info);
}

void
phobosTest::testPhobosBeale()
{
  double *par, *mea, *Bea_Par, *work, *covar;
  int npar, nmea, worksize;

  // Calculate worksize
  npar     = PHOBOS_BEALE_NPAR;
  nmea     = PHOBOS_BEALE_NMEA;
  worksize = PHOBOS_WORKSIZE(npar, nmea);

  // Allocate and clean the workspace
  double *phobos_workspace = (double *) malloc(worksize * sizeof(double));
  for (int i = 0; i < worksize; ++i)
    phobos_workspace[i] = .0;
  double *opts = (double *) malloc(4 * sizeof(double));
  double *info = (double *) malloc(PHOBOS_INFOSIZE * sizeof(double));
  Bea_Par      = (double *) malloc(PHOBOS_BEALE_PARAMETERS * sizeof(double));

  // Choose values to force LM to the optimum.
  opts[0] = 1e-10;
  opts[1] = opts[2] = opts[3] = 7e-27;

  // Distribute the workspace
  par   = phobos_workspace;
  mea   = par + npar;
  work  = mea + nmea;
  covar = par + worksize - PHOBOS_BEALE_COVAR;

  // Fill the needed function values
  Bea_Par[0] = PHOBOS_BEALE_A;
  Bea_Par[1] = PHOBOS_BEALE_B;
  Bea_Par[2] = PHOBOS_BEALE_C;

  par[0] = PHOBOS_BEALE_X_START;
  par[1] = PHOBOS_BEALE_Y_START;

  // Run and evaluate LM
  phobos(Beale_4_phobos, jac_Beale_4_phobos, par, mea, npar, nmea, 10000, opts,
         info, work, covar, (void *) Bea_Par);
  CPPUNIT_ASSERT(fabs(par[0] - PHOBOS_BEALE_X_END) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(par[1] - PHOBOS_BEALE_Y_END) < ALLOWED_ERROR);

  // Clean the workspace
  par = mea = work = covar = NULL;
  free(Bea_Par);
  free(phobos_workspace);
  free(opts);
  free(info);
}

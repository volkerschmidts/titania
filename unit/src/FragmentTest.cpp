
#include <Fragment.hpp>
#include <FragmentTest.h>
// there is still some weird include ordering dependence hidden in
// Declarations.hpp
// clang-format off
#include <Declarations.hpp>
// clang-format on

#ifndef ALLOWED_ERROR
#define ALLOWED_ERROR 1e-3
#endif

#ifndef FRAGMENT_CALCULATE_VALUES // calculate_x and calculate_x_derivative
                                  // values
#define FRAGMENT_CALCULATE_VALUES
// Vales are the respective coordinates of Fusicoccin_A
double C1[3] = {5.66321, 3.26927, 4.61031};
double C2[3] = {4.66908, 5.57021, 5.09261};
double C3[3] = {5.00112, 4.52926, 4.04117};
double C4[3] = {3.56977, 4.19106, 3.55960};
double C5[3] = {2.70972, 4.31542, 4.80645};
double C6[3] = {3.40229, 5.42375, 5.53272};
double C7[3] = {2.67138, 6.18540, 6.62672};

double r_C1C3 = 1.5329;
double r_C2C3 = 1.5164;
double r_C3C4 = 1.5476;
double r_C4C5 = 1.5198;
double r_C5C6 = 1.4952;
double r_C6C7 = 1.5203;

double a_C1C3C2 = 1.9838;
double a_C2C3C4 = 1.7348;
double a_C3C4C5 = 1.8237;
double a_C4C5C6 = 1.7691;
double a_C5C6C7 = 2.0924;

// Arbitrary angles in radians
// This angles are used as polar coordinates for rdcs
double deg_n20 = -0.3491;
double deg_0   = 0.0;
double deg_5   = 0.0873;
double deg_20  = 0.3491;
double deg_90  = 1.5708;
double deg_160 = 2.7925;
double deg_175 = 3.0543;

double rdc_C1C3_0_0     = 1.1904;
double rdc_C1C3_90_90   = 0.6060;
double rdc_C5C6_160_n20 = 1.1655;
#endif

#ifndef FRAGMENT_CALCULATE_DERIVATIVES
#define FRAGMENT_CALCULATE_DERIVATIVES
// Diviations obtained by Matlab (numerical deviations)
double drC1C3dx3 = -0.4319;
double drC1C3dz1 = 0.3713;

double drC2C3dx2 = -0.2190;
double drC2C3dy3 = -0.6865;

double drC3C4dz3 = 0.3112;
double drC3C4dy4 = -0.2185;

double drC4C5dx4 = 0.5659;
double drC4C5dy4 = -0.0818;

double dC1C3C2dr[9] = {0.0325, 0.2152,  -0.2477, -0.2540, -0.1395,
                       0.3934, -0.6000, 1.0677,  -0.4677};

double dC5C6C7dr[9] = {0.5489, -1.0820, 0.5331, -0.1016, -0.2714,
                       0.3730, -0.3684, 0.2719, 0.0965};


double dC5C6_0_0dr[6] = {0.1127, -0.1127, -0.2144, 0.2144, -0.6057, 0.6057};

double dC5C6_90_90dr[6] = {0.4066, -0.4066, 0.3716, -0.3715, 0.3495, -0.3495};

double dC5C6_160_n20dr[6] = {0.3668, -0.3668, 0.1276, -0.1276, -0.5445, 0.5445};

double dC5C6_20_160dr[6] = {-0.3668, 0.3668, -0.1276, 0.1276, 0.5445, -0.5445};

#endif

CPPUNIT_TEST_SUITE_REGISTRATION(FragmentTest);

FragmentTest::FragmentTest()
{
  ;
}

void
FragmentTest::testcalculate_distance()
{
  CPPUNIT_ASSERT(fabs(Fragment::calculate_distance(&(C1[0]), &(C3[0])) -
                      r_C1C3) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(Fragment::calculate_distance(&(C2[0]), &(C3[0])) -
                      r_C2C3) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(Fragment::calculate_distance(&(C3[0]), &(C4[0])) -
                      r_C3C4) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(Fragment::calculate_distance(&(C4[0]), &(C5[0])) -
                      r_C4C5) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(Fragment::calculate_distance(&(C5[0]), &(C6[0])) -
                      r_C5C6) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(Fragment::calculate_distance(&(C6[0]), &(C7[0])) -
                      r_C6C7) < ALLOWED_ERROR);
}

void
FragmentTest::testcalculate_distance_derivative()
{
  std::cout << "   - Subruns: Running 8 deviations on 4 vectors...\n";

  CPPUNIT_ASSERT(fabs(Fragment::calculate_distance_derivative(
                          &(C1[0]), &(C3[0]), &(C3[0]), STANDARD_AXIS_X_) -
                      drC1C3dx3) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(Fragment::calculate_distance_derivative(
                          &(C1[0]), &(C3[0]), &(C1[0]), STANDARD_AXIS_Z_) -
                      drC1C3dz1) < ALLOWED_ERROR);

  CPPUNIT_ASSERT(fabs(Fragment::calculate_distance_derivative(
                          &(C2[0]), &(C3[0]), &(C2[0]), STANDARD_AXIS_X_) -
                      drC2C3dx2) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(Fragment::calculate_distance_derivative(
                          &(C2[0]), &(C3[0]), &(C3[0]), STANDARD_AXIS_Y_) -
                      drC2C3dy3) < ALLOWED_ERROR);

  CPPUNIT_ASSERT(fabs(Fragment::calculate_distance_derivative(
                          &(C3[0]), &(C4[0]), &(C3[0]), STANDARD_AXIS_Z_) -
                      drC3C4dz3) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(Fragment::calculate_distance_derivative(
                          &(C3[0]), &(C4[0]), &(C4[0]), STANDARD_AXIS_Y_) -
                      drC3C4dy4) < ALLOWED_ERROR);

  CPPUNIT_ASSERT(fabs(Fragment::calculate_distance_derivative(
                          &(C4[0]), &(C5[0]), &(C4[0]), STANDARD_AXIS_X_) -
                      drC4C5dx4) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(Fragment::calculate_distance_derivative(
                          &(C4[0]), &(C5[0]), &(C4[0]), STANDARD_AXIS_Y_) -
                      drC4C5dy4) < ALLOWED_ERROR);
}

void
FragmentTest::testcalculate_angle()
{
  CPPUNIT_ASSERT(fabs(Fragment::calculate_angle(&(C1[0]), &(C3[0]), &(C2[0])) -
                      a_C1C3C2) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(Fragment::calculate_angle(&(C2[0]), &(C3[0]), &(C4[0])) -
                      a_C2C3C4) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(Fragment::calculate_angle(&(C3[0]), &(C4[0]), &(C5[0])) -
                      a_C3C4C5) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(Fragment::calculate_angle(&(C4[0]), &(C5[0]), &(C6[0])) -
                      a_C4C5C6) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(Fragment::calculate_angle(&(C5[0]), &(C6[0]), &(C7[0])) -
                      a_C5C6C7) < ALLOWED_ERROR);
}

void
FragmentTest::testcalculate_angle_derivative()
{
  std::cout
      << "   - Subruns: Running 18 deviations on 4 vectors in 2 angles...\n";
  unsigned int i;
  for (i = 0; i < 3; ++i)
  {
    CPPUNIT_ASSERT(fabs(Fragment::calculate_angle_derivative(
                            &(C1[0]), &(C3[0]), &(C2[0]), &(C1[0]), i) -
                        dC1C3C2dr[3 * i]));
    CPPUNIT_ASSERT(fabs(Fragment::calculate_angle_derivative(
                            &(C1[0]), &(C3[0]), &(C2[0]), &(C3[0]), i) -
                        dC1C3C2dr[3 * i + 1]));
    CPPUNIT_ASSERT(fabs(Fragment::calculate_angle_derivative(
                            &(C1[0]), &(C3[0]), &(C2[0]), &(C2[0]), i) -
                        dC1C3C2dr[3 * i + 2]));
  }
  for (i = 0; i < 3; ++i)
  {
    CPPUNIT_ASSERT(fabs(Fragment::calculate_angle_derivative(
                            &(C5[0]), &(C6[0]), &(C7[0]), &(C5[0]), i) -
                        dC5C6C7dr[3 * i]));
    CPPUNIT_ASSERT(fabs(Fragment::calculate_angle_derivative(
                            &(C5[0]), &(C6[0]), &(C7[0]), &(C6[0]), i) -
                        dC5C6C7dr[3 * i + 1]));
    CPPUNIT_ASSERT(fabs(Fragment::calculate_angle_derivative(
                            &(C5[0]), &(C6[0]), &(C7[0]), &(C7[0]), i) -
                        dC5C6C7dr[3 * i + 2]));
  }
}

void
FragmentTest::testcalculate_rdc_angle()
{
  double theta, phi;
  theta = deg_0;
  phi   = deg_0;
  CPPUNIT_ASSERT(
      fabs(Fragment::calculate_rdc_angle(&(C1[0]), &(C3[0]), theta, phi) -
           rdc_C1C3_0_0) < ALLOWED_ERROR);
  theta = deg_90;
  phi   = deg_90;
  CPPUNIT_ASSERT(
      fabs(Fragment::calculate_rdc_angle(&(C1[0]), &(C3[0]), theta, phi) -
           rdc_C1C3_90_90) < ALLOWED_ERROR);
  theta = deg_160;
  phi   = deg_n20;
  CPPUNIT_ASSERT(
      fabs(Fragment::calculate_rdc_angle(&(C5[0]), &(C6[0]), theta, phi) -
           rdc_C5C6_160_n20) < ALLOWED_ERROR);
  theta = deg_20;
  phi   = deg_160;
  CPPUNIT_ASSERT(
      fabs(Fragment::calculate_rdc_angle(&(C5[0]), &(C6[0]), theta, phi) -
           rdc_C5C6_160_n20) < ALLOWED_ERROR);
}

void
FragmentTest::testcalculate_rdc_angle_derivative()
{
  double theta, phi;
  unsigned int i;
  theta = deg_0;
  phi   = deg_0;
  for (i = 0; i < 3; ++i)
  {
    CPPUNIT_ASSERT(fabs(Fragment::calculate_rdc_angle_derivative(
                            &(C1[0]), &(C3[0]), &(C1[0]), theta, phi, i) -
                        dC5C6_0_0dr[2 * i]) < ALLOWED_ERROR);
    CPPUNIT_ASSERT(fabs(Fragment::calculate_rdc_angle_derivative(
                            &(C1[0]), &(C3[0]), &(C3[0]), theta, phi, i) -
                        dC5C6_0_0dr[2 * i + 1]) < ALLOWED_ERROR);
  }

  theta = deg_90;
  phi   = deg_90;
  for (i = 0; i < 3; ++i)
  {
    CPPUNIT_ASSERT(fabs(Fragment::calculate_rdc_angle_derivative(
                            &(C1[0]), &(C3[0]), &(C1[0]), theta, phi, i) -
                        dC5C6_90_90dr[2 * i]) < ALLOWED_ERROR);
    CPPUNIT_ASSERT(fabs(Fragment::calculate_rdc_angle_derivative(
                            &(C1[0]), &(C3[0]), &(C3[0]), theta, phi, i) -
                        dC5C6_90_90dr[2 * i + 1]) < ALLOWED_ERROR);
  }

  theta = deg_160;
  phi   = deg_n20;
  for (i = 0; i < 3; ++i)
  {
    CPPUNIT_ASSERT(fabs(Fragment::calculate_rdc_angle_derivative(
                            &(C5[0]), &(C6[0]), &(C5[0]), theta, phi, i) -
                        dC5C6_20_160dr[2 * i]) < ALLOWED_ERROR);
    CPPUNIT_ASSERT(fabs(Fragment::calculate_rdc_angle_derivative(
                            &(C5[0]), &(C6[0]), &(C6[0]), theta, phi, i) -
                        dC5C6_20_160dr[2 * i + 1]) < ALLOWED_ERROR);
  }

  theta = deg_20;
  phi   = deg_160;
  for (i = 0; i < NUMBER_OF_AXIS_; ++i)
  {
    CPPUNIT_ASSERT(fabs(Fragment::calculate_rdc_angle_derivative(
                            &(C5[0]), &(C6[0]), &(C5[0]), theta, phi, i) -
                        dC5C6_20_160dr[2 * i]) < ALLOWED_ERROR);
    CPPUNIT_ASSERT(fabs(Fragment::calculate_rdc_angle_derivative(
                            &(C5[0]), &(C6[0]), &(C6[0]), theta, phi, i) -
                        dC5C6_20_160dr[2 * i + 1]) < ALLOWED_ERROR);
  }
}


double A[15] = {-4.828, 9.752, -0.182, 7.918,  -1.931, 1.293, -5.216, 1.743,
                4.341,  0.136, -3.549, -1.287, 0.123,  0.123, 0.123};

double dMPI_5x3[15] = {0.0157,  0.1084, -0.0426, -0.0367, 0.0036,
                       0.0974,  0.0444, -0.0240, -0.0499, 0.0029,
                       -0.0345, 0.0994, 0.1815,  -0.0595, 0.0064};

void
FragmentTest::mpi_ff_test()
{
  double ainv[15];
  svdMoorePenroseInverse(A, 5, 3, ainv);
  for (int j = 0; j < 15; ++j)
  {
    //      std::cout << ainv[j*5+i] << "   " << dMPI_5x3[j*5+i] << std::endl;
    CPPUNIT_ASSERT(fabs(ainv[j] - dMPI_5x3[j]) < ALLOWED_ERROR);
  }
}

void
FragmentTest::testperform_optimization_step()
{
  ;
}

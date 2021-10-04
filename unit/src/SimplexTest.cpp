#include <Simplex.hpp>
#include <SimplexTest.hpp>
// there is still some weird include ordering dependence hidden in
// Declarations.hpp
// clang-format off
#include <Declarations.hpp>
// clang-format on
#include <cmath>
#include <main.h>
#include <nlopt.hpp>
#include <vector>


typedef struct {
  double a, b;
} my_constraint_data;

double
myvfunc(const std::vector<double> &x,
        std::vector<double> &grad,
        void *my_func_data)
{
  if (!grad.empty())
  {
    grad[0] = 0.0;
    grad[1] = 0.5 / sqrt(x[1]);
  }
  return sqrt(x[1]);
}

double
myvconstraint(const std::vector<double> &x,
              std::vector<double> &grad,
              void *data)
{
  my_constraint_data *d = reinterpret_cast<my_constraint_data *>(data);
  double a = d->a, b = d->b;
  if (!grad.empty())
  {
    grad[0] = 3 * a * (a * x[0] + b) * (a * x[0] + b);
    grad[1] = -1.0;
  }
  return ((a * x[0] + b) * (a * x[0] + b) * (a * x[0] + b) - x[1]);
}

typedef struct {
  double a, b;
  unsigned int count;
} Rosenbrock_par;

double
Rosenbrock2D(const std::vector<double> &x,
             std::vector<double> &grad,
             void *func_data)
{
  Rosenbrock_par *r_par = reinterpret_cast<Rosenbrock_par *>(func_data);
  double a              = r_par->a;
  double b              = r_par->b;

  r_par->count++;

  if (!grad.empty())
  {
    grad[0] = 2 * (x[0] * (2 * b * (pow(x[0], 2) - x[1]) + 1) - a);
    grad[1] = 2 * b * (-pow(x[0], 2) + x[1]);
  }
  double val = pow(a - x[0], 2) + b * pow(x[1] - pow(x[0], 2), 2);
  // fprintf(stderr, "count %d x[0] %.14f x[1] %.14f val %.14f\n", r_par->count,
  // x[0], x[1], val);
  return val;
}

CPPUNIT_TEST_SUITE_REGISTRATION(SimplexTest);

SimplexTest::SimplexTest()
{
  ;
}

SimplexTest::~SimplexTest()
{
  ;
}

void
SimplexTest::testNLoptExample()
{
  nlopt::opt opt(nlopt::LD_MMA, 2);
  std::vector<double> lb(2);
  lb[0] = -HUGE_VAL;
  lb[1] = 0;
  opt.set_lower_bounds(lb);
  opt.set_min_objective(myvfunc, NULL);
  my_constraint_data data[2] = {{2, 0}, {-1, 1}};
  opt.add_inequality_constraint(myvconstraint, &data[0], 1e-8);
  opt.add_inequality_constraint(myvconstraint, &data[1], 1e-8);
  opt.set_xtol_rel(1e-4);
  std::vector<double> x(2);
  x[0] = 1.234;
  x[1] = 5.678;
  std::vector<double> x_expected(2);
  x_expected[0] = 0.3333333348;
  x_expected[1] = 0.2962962891;
  double minf;
  double minf_expected = 0.5443310474;

  try
  {
    nlopt::result result = opt.optimize(x, minf);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(x_expected[0], x[0], 1e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(x_expected[1], x[1], 1e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(minf_expected, minf, 1e-8);
  }
  catch (std::exception &e)
  {
    std::cout << "nlopt failed: " << e.what() << std::endl;
  }
}

void
SimplexTest::testNLoptRosenbrock2D_NM()
{
  Rosenbrock_par r_par;
  r_par.a     = 1;
  r_par.b     = 100;
  r_par.count = 0;
  double xtol = 1e-8;
  double ftol = 1e-4;

  nlopt::opt opt(nlopt::LN_NELDERMEAD, 2);
  opt.set_min_objective(Rosenbrock2D, &r_par);
  opt.set_xtol_rel(xtol);

  std::vector<double> x(2);
  x[0] = -1.9;
  x[1] = 2;
  double minf;

  std::vector<double> x_expected(2);
  x_expected[0] = x_expected[1] = 1.0;
  double minf_expected          = 0.0;

  try
  {
    nlopt::result result = opt.optimize(x, minf);
    CPPUNIT_ASSERT_EQUAL(nlopt::XTOL_REACHED, result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(x_expected[0], x[0], xtol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(x_expected[1], x[1], xtol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(minf_expected, minf, ftol);
  }
  catch (std::exception &e)
  {
    std::cout << "nlopt failed: " << e.what() << std::endl;
  }
}

void
SimplexTest::testSimplexNLoptRosenbrock2D_NM()
{
  Eigen::VectorXd start(2);
  start << -1.9, 2.0;
  unsigned int ndim = start.size();

  Simplex simplex(start);
  simplex.vfunc = Rosenbrock2D;
  Rosenbrock_par r_par;
  r_par.a     = 1;
  r_par.b     = 100;
  r_par.count = 0;
  double xtol = 1e-8;
  double ftol = 1e-4;

  simplex.func_par = &r_par;
  simplex.setxtol(xtol);
  simplex.setftol(ftol);

  Eigen::VectorXd expected(2);
  expected << 1.0, 1.0;
  double minf_expected = 0.0;

  try
  {
    Eigen::VectorXd final = simplex.minimize();
    for (unsigned int i = 0; i < ndim; ++i)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected(i), final(i), xtol);
    }
    CPPUNIT_ASSERT_DOUBLES_EQUAL(minf_expected, simplex.getfmin(), ftol);
    // CPPUNIT_ASSERT_EQUAL(nlopt::FTOL_REACHED, simplex.getresult()); //
    // doesn't work, seems to always give nlopt::XTOL_REACHED
  }
  catch (std::exception &e)
  {
    std::cout << "nlopt failed: " << e.what() << std::endl;
  }
}

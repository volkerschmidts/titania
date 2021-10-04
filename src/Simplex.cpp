#include <Declarations.hpp>
#include <Simplex.hpp>
#include <StructureSimulator.hpp>
#include <iostream>
#include <nlopt.hpp>
#include <vector>

Simplex::Simplex(Eigen::VectorXd &point)
{
  ndim         = point.size();
  ftol         = 5e-5;
  fmin         = HUGE_VAL;
  xtol         = 1e-10;
  NMAX         = 1000;
  nfunc        = 0;
  initial_step = 0.01;
  x            = std::vector<double>(point.data(), point.data() + point.size());
}

Eigen::VectorXd
Simplex::minimize()
{
  nlopt::opt opt(nlopt::LN_NELDERMEAD, ndim);
  opt.set_min_objective(vfunc, func_par);
  opt.set_xtol_rel(xtol);
  opt.set_ftol_rel(ftol);
  opt.set_maxeval(NMAX);
  opt.set_initial_step(initial_step);

  try
  {
    result = opt.optimize(x, fmin);
  }
  catch (std::exception &e)
  {
    std::cout << "nlopt failed: " << e.what() << std::endl;
  }
  nfunc = opt.get_numevals();

  return Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(x.data(), x.size());
}

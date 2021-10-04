#ifndef SIMPLEX_HPP_
#define SIMPLEX_HPP_

#include <eigen3/Eigen/Core>
#include <nlopt.hpp>
#include <vector>

#ifndef STRUCTURE_SIMULATOR_HPP_
class StructureSimulator;
#endif

#ifndef MOLECULE_HPP_
class Molecule;
#endif

class Simplex {
 public:
  Simplex(Eigen::VectorXd &);
  Eigen::VectorXd minimize();

  void setNMAX(unsigned int nmx)
  {
    NMAX = nmx;
  }
  unsigned int getNMAX() const
  {
    return NMAX;
  }

  unsigned int getnfunc() const
  {
    return nfunc;
  }

  void setftol(double tol)
  {
    ftol = tol;
  }
  double getftol() const
  {
    return ftol;
  }

  void setxtol(double tol)
  {
    xtol = tol;
  }
  double getxtol() const
  {
    return xtol;
  }

  void setinitialstep(double d)
  {
    initial_step = d;
  }

  double getfmin() const
  {
    return fmin;
  }

  void *func_par;
  double (*vfunc)(const std::vector<double> &, std::vector<double> &, void *);
  nlopt::result getresult() const
  {
    return result;
  }

 private:
  std::vector<double> x;
  double xtol;
  double fmin;
  double ftol;
  unsigned int NMAX;
  unsigned int nfunc;
  unsigned int ndim;

  nlopt::result result;
  double initial_step;
};

#endif
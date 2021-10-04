
#ifndef DIDC_H_
#define DIDC_H_

#include <eigen3/Eigen/Core>

#ifndef STRUCTURE_HPP_
class Structure;
#endif

#ifndef SPHERICALHARMONICS_HPP_
class SphericalHarmonics;
#endif

struct B_opt_data {
  double *UD;
  double *UDt;
  double *Ident;
  double *work;
  double *sol;
  int NOR;
};

int calculate_lambda_B(Eigen::MatrixXd &, Eigen::MatrixXd &, Eigen::MatrixXd &);

int refine_B_Matrix(Eigen::MatrixXd &, Eigen::MatrixXd &, Eigen::MatrixXd &);

int analyse_B_ref(Eigen::MatrixXd &, SphericalHarmonics *);

int compute_B_Angles(Structure *);

void LambdaOpt(double *, double *, int, int, void *);

#endif

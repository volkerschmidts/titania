#ifndef MC_OUTPUT_TITANIA_
#define MC_OUTPUT_TITANIA_

#include <eigen3/Eigen/Core>

struct MC_Output {
  Eigen::MatrixXd p_mean;
  Eigen::MatrixXd p_trig_mean;
  Eigen::MatrixXd p_sigm;

  Eigen::VectorXd R_2_mean;
  Eigen::VectorXd R_mean;
  Eigen::VectorXd circular_sigma;

  Eigen::MatrixXd D_mean;
  Eigen::MatrixXd D_sigm;

  Eigen::MatrixXd D_calc_mean;
  Eigen::MatrixXd D_calc_sigm;

  Eigen::MatrixXd Aligns;
  Eigen::MatrixXd Aligns_sigm;

  Eigen::MatrixXd Saupe_tensor;
  Eigen::MatrixXd Saupe_tensor_sigm;

  Eigen::MatrixXd Evals;
  Eigen::MatrixXd Evals_sigm;

  Eigen::MatrixXd Evecs;
  Eigen::MatrixXd Evecs_sigm;

  Eigen::MatrixXd Euler;
  Eigen::MatrixXd Euler_sigm;

  unsigned int steps;
};
#endif
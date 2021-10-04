
#ifndef INDEPENDENCY_ANALYSIS_HPP_
#define INDEPENDENCY_ANALYSIS_HPP_

#include <eigen3/Eigen/Core>

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
enum struct StructureOptions : unsigned int;
#endif

#ifndef SECONDA_OUTPUT_TITANIA_
struct SECONDA_Output;
#endif

#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef STRUCTURE_HPP_
class Structure;
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

SECONDA_Output
SECONDA_analysis(Molecule &, Structure *, BasicInformation &, Flags &);

Eigen::MatrixXd
SECONDA_CovarianceMatrix(Eigen::MatrixXd &, Molecule &, Flags &);

Eigen::MatrixXd SECONDA_CovarianceMatrix_reduced(Eigen::MatrixXd &, Molecule &);

Eigen::MatrixXd
SECONDA_CovarianceMatrix_Hus(Eigen::MatrixXd &); // RDC_normalized

Eigen::VectorXd SECONDA_a_2(Eigen::MatrixXd &, Eigen::VectorXd &);

Eigen::VectorXd SECONDA_Kappa(Eigen::MatrixXd &); // Covariance eigenvectors


double SECONDA_rho(Eigen::MatrixXd &); // Covariance eigenvalues


double getConditionNumber(Eigen::MatrixXd &,
                          unsigned int &,
                          double,
                          unsigned int base = 0);


Eigen::VectorXd getSingularValues(Eigen::MatrixXd &);


Eigen::MatrixXd getTwoSetSingularValues(Eigen::MatrixXd &);


SECONDA_Output generate_Seconda_Output(Eigen::MatrixXd &RDCnorm,
                                       Eigen::VectorXd &Evals,
                                       Eigen::MatrixXd &Evecs,
                                       Eigen::VectorXd &KappaQ,
                                       Eigen::MatrixXd &R_2,
                                       Eigen::MatrixXd &m,
                                       Eigen::MatrixXd &b,
                                       Eigen::VectorXd &a_2,
                                       BasicInformation &baseInformation);


void full_Regression(Eigen::MatrixXd &, // RDC_normalized
                     Eigen::MatrixXd &, // Pearson R^2
                     Eigen::MatrixXd &, // slopes
                     Eigen::MatrixXd &  // intercepts
);


void rdc_vector_analysis(Structure *, StructureOptions);


void SECONDA_sensitivity(Molecule &, Structure *, BasicInformation &, Flags &);

#endif

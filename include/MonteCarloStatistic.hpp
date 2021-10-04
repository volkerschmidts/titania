
#ifndef MONTECARLOSTATISTIC_HPP_
#define MONTECARLOSTATISTIC_HPP_

#include <eigen3/Eigen/Core>
#include <omp.h>

#define MC_MINSTEPS_TITANIA 1000


#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef STRUCTURE_HPP_
class Structure;
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

#ifndef MC_RUN_INFORMATION_TITANIA_
struct MC_Run_Information;
#endif

#ifndef MC_OUTPUT_TITANIA_
#include <MC_Output.hpp>
#endif

class MonteCarloStep {
 private:
  unsigned int step;
  unsigned int NOS;
  unsigned int NOR;

  bool numericalDeviation;
  bool useMatrixFormalism;

  // RDC Matrizes
  Eigen::MatrixXd RDCmatrix;
  Eigen::MatrixXd RDCmatrixScaled;
  Eigen::MatrixXd DeltaRDC;
  Eigen::MatrixXd DeltaRDCScaled;
  Eigen::MatrixXd randMat;
  Eigen::MatrixXd RDCVar;
  Eigen::MatrixXd RDCVarScaled;
  Eigen::MatrixXd RDC_Calc;
  Eigen::MatrixXd Dmax;

  // Molecular Orientation Matrizes
  Eigen::MatrixXd A;
  Eigen::MatrixXd S;
  Eigen::MatrixXd EulerAngles;
  Eigen::MatrixXd EigenVals;
  Eigen::MatrixXd EigenVecs;
  Eigen::MatrixXcd F;
  Eigen::MatrixXcd F_i;

  // Structural Information matrizes
  Eigen::MatrixXd C;
  Eigen::MatrixXcd Y;
  Eigen::MatrixXd p;
  Eigen::MatrixXd direction;
  Eigen::MatrixXd p_trig;
  Eigen::MatrixXd weights;

  //      Eigen::MatrixXd polar2cartesian;

  MonteCarloStep *next;
  MonteCarloStep *prev;

 public:
  MonteCarloStep(Molecule &, Structure *, Flags &);
  MonteCarloStep(MonteCarloStep *);
  ~MonteCarloStep();

  void SetupMatrizes(Structure *);
  void VaryMatrix();
  static Eigen::MatrixXd MC_randMatrix(unsigned int &, unsigned int &);

  MonteCarloStep *addNextStep();
  MonteCarloStep *getNext()
  {
    return next;
  }
  MonteCarloStep *getPrev()
  {
    return prev;
  }
  void run(BasicInformation &, Eigen::MatrixXd &);

  Eigen::MatrixXd getA() const
  {
    return A;
  }
  Eigen::MatrixXd getS() const
  {
    return S;
  }
  Eigen::MatrixXd getEigenVecs() const
  {
    return EigenVecs;
  }
  Eigen::MatrixXd getD_calc() const
  {
    return RDC_Calc;
  }
  Eigen::MatrixXd getD() const
  {
    return RDCVar;
  }
  Eigen::MatrixXd getDeltaD() const
  {
    return DeltaRDC;
  }
  Eigen::MatrixXd getStatP();
  Eigen::MatrixXd getStatP_trig() const
  {
    return p_trig;
  }
  Eigen::MatrixXd getStat_direction() const
  {
    return direction;
  }
  Eigen::MatrixXd getWeights() const
  {
    return weights;
  }

  bool useMatrizes() const
  {
    return useMatrixFormalism;
  }
  void encrypt();
  void cleanup();
};

class MonteCarloStatistic {
 private:
  unsigned int mc_steps;
  unsigned int stan_steps;
  unsigned int NOS;
  unsigned int NOR;
  unsigned int NOR2;
  double MCSd;
  double rmsd;

  Eigen::MatrixXd MCPi;
  Eigen::MatrixXd A;
  Eigen::MatrixXd A_sigm;
  Eigen::MatrixXd S;
  Eigen::MatrixXd S_sigm;
  Eigen::MatrixXd D;
  Eigen::MatrixXd D_sigm;
  Eigen::MatrixXd D_calc;
  Eigen::MatrixXd D_calc_sigm;
  Eigen::MatrixXd p;
  Eigen::MatrixXd direction_mean;
  Eigen::MatrixXd direction_cov;
  Eigen::MatrixXd p_sigm;
  Eigen::MatrixXd p_trig;
  Eigen::VectorXd R_2_mean;
  Eigen::VectorXd R_mean;
  Eigen::VectorXd circular_sigma;
  Eigen::MatrixXd dD;

  omp_lock_t A_lock;
  omp_lock_t A_sigm_lock;
  omp_lock_t S_lock;
  omp_lock_t S_sigm_lock;
  omp_lock_t D_lock;
  omp_lock_t D_sigm_lock;
  omp_lock_t D_calc_lock;
  omp_lock_t D_calc_sigm_lock;
  omp_lock_t p_lock;
  omp_lock_t p_sigm_lock;

  MonteCarloStep *HeadStep;

 public:
  MonteCarloStatistic();
  ~MonteCarloStatistic();
  int startMonteCarlo(Molecule &,
                      Structure *,
                      BasicInformation &,
                      Flags &,
                      MC_Run_Information &);
  void fetchValues(MonteCarloStep *);
  void reduceStatistic();
  void decrypt();

  void addA_sigm(MonteCarloStep *);
  void addD_sigm(MonteCarloStep *);
  void addD_calc_sigm(MonteCarloStep *);
  void addp_sigm(MonteCarloStep *);

  void reduceSigmas();
  int checkConsistency(MC_Run_Information &);
  void plot_mc_angles(Molecule &, Structure *, BasicInformation &, Flags &);
  void print_mc_polar_angles(Structure *, BasicInformation &);

  MC_Output getStatistics();
  void cleanup();
};


double phiOutput(double p);

Eigen::MatrixXd getMonteCarloD(Molecule &);

void scaleMonteCarloD(Molecule &, Eigen::MatrixXd &, Eigen::MatrixXd &);

void getMonteCarloAngles(Eigen::MatrixXd &,
                         Eigen::MatrixXd &,
                         Eigen::MatrixXcd &,
                         unsigned int,
                         BasicInformation &,
                         bool);


Eigen::MatrixXd estimateSaupeErrors(
    Eigen::MatrixXd,   // vec(A)
    Eigen::MatrixXd &, // ten(A)
    Eigen::MatrixXd &, // eigenvectors
    Eigen::MatrixXd &, // eigenvalues
    Eigen::Vector3d &, // Euler angles
    Eigen::MatrixXd,   // sigm(vec(A))
    // see .cpp file.                      Eigen::MatrixXd&, // sigm(ten(A))
    Eigen::MatrixXd &, // sigm(eigenvectors)
    Eigen::MatrixXd &, // sigm(eigenvalues)
    Eigen::Vector3d &  // sigm(Euler angles)
);

#endif


#ifndef LINALG_HPP_
#define LINALG_HPP_

#include <complex>
#include <eigen3/Eigen/Core>

#ifndef ZERO_CUTOFF_
#define ZERO_CUTOFF_ 1e-9
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef SPHERICALHARMONICS_HPP_
class SphericalHarmonics;
#endif

#ifndef ATOM_HPP_
class Atom;
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
enum struct StructureOptions : unsigned int;
#endif

extern const double fac4pi;

class Tensor {
 private:
  Eigen::Matrix3d T;     /* Tensor */
  Eigen::Vector3d Evals; /* Eigen values */
  Eigen::Matrix3d Evecs; /* Eigen vectors */
  Eigen::MatrixXd Euler;

 public:
  /* Standard constructor */
  Tensor();
  /* Constructor whith in itialization of T & full determination of eigen system
   */
  Tensor(Eigen::Matrix3d &);

  /* Check if destructor is needed later */
  //~Tensor();

  /* Accessing T */
  void setTensor(Eigen::Matrix3d &);
  Eigen::Matrix3d getTensor() const
  {
    return T;
  };

  /* Just getter for Eves/Evals. Nothing bad can happen if only Tensor can
   * manipulate them */
  Eigen::Vector3d getEvals() const
  {
    return Evals;
  };
  Eigen::Matrix3d getEvecs() const
  {
    return Evecs;
  };


  /* Determine & access euler angles */
  // Eigen::MatrixXd getEulerAngles ();   /* TODO define a flag to say if
  // all/positive/one defined line is needed */ void determineEuler ();


  /* Determine eigen system */
  // void determineEigenSystem ();
};

int Eigen2Double(Eigen::MatrixXd &, double **);
int Eigen2Double(Eigen::MatrixXd &, double **, const int, const int);
int Double2Eigen(Eigen::MatrixXd &, double *);
int Double2Eigen(Eigen::MatrixXd &, double *, const int, const int);

double normRand(double, // mu
                double  // sigma
);


void A2goodA(double &, // theta
             double &  // phi
);


template<typename X>
Eigen::Matrix<X, -1, -1, 0, -1, -1>
MoorePenroseInverse(Eigen::Matrix<X, -1, -1, 0, -1, -1> &,
                    double cutoff   = ZERO_CUTOFF_,
                    bool enable_gpu = true);

#ifdef USE_CUDA
Eigen::MatrixXd GPU_MoorePenroseInverse(Eigen::MatrixXd &, double);
#endif

Eigen::MatrixXcd MoorePenroseInverse(Eigen::MatrixXcd &,
                                     double cutoff = ZERO_CUTOFF_);
Eigen::MatrixXcd get_svd_U(Eigen::MatrixXcd &, bool full = true);
Eigen::MatrixXd get_svd_U(Eigen::MatrixXd &, bool full = true);

/* Scales and normalizes the rdc matrix based on Dmax */
void ScaleRDCMatrix(Eigen::MatrixXd &, // RDC
                    Eigen::MatrixXd &, // RDC_scaled
                    Eigen::MatrixXd &, // RDC_normalized
                    Molecule &,
                    BasicInformation &,
                    StructureOptions);


/* Calculates spherical harmonics from plgndr */
std::complex<double> Ylm(const int,    // l
                         const int,    // m
                         const double, // theta
                         const double  // phi
);


/* Factorial = Product ( i ) */
double fac(const int a);


/* Calculates B row based on polar/euler angles theta, phi, alpha, beta & gamma
 */
Eigen::MatrixXd Sphere2B(const double, // theta
                         const double, // phi
                         const double alpha = .0,
                         const double beta  = .0,
                         const double gamma = .0);


/* Calculates small wigner d from l, m & beta */
double d2Mm(const int,   // l
            const int,   // m
            const double // beta
);

double dd2Mmdb(const int,   // l
               const int,   // m
               const double // beta
);


/* Calculates big wigner D from l, m, alpha, beta & gamma */
std::complex<double> D2Mm(const int,    // l
                          const int,    // m
                          const double, // alpha
                          const double, // beta
                          const double  // gamma
);

/* Derivative of D2Mm with respect to beta. Unused sice LM runs proper */
std::complex<double> DeltaD2MmDeltaBeta(const int,    // l
                                        const int,    // m
                                        const double, // alpha
                                        const double, // beta
                                        const double  // gamma
);


/* Calculates the rotated spherical harmonic based on m, theta, phi, alpha, beta
 * & gamma */
std::complex<double> D2MmY2m(const int,    // m
                             const double, // theta
                             const double, // phi
                             const double, // alpha
                             const double, // beta
                             const double  // gamma
);

void motionalAnalysisOfSphericalHarmonics(Eigen::MatrixXcd, // Y
                                          double,           // theta_av
                                          double,           //   phi_av
                                          double &,         // S^2 (rdc)
                                          double &,         // S^2 (ax)
                                          double &,         // eta (rdc)
                                          double &          // phi'(rdc)
);

double SrdcBySoverall(Eigen::MatrixXd);


/* Decomposition of B-matrix with calculation of A */
Eigen::MatrixXd BD2A(Eigen::MatrixXd &, // B
                     Eigen::MatrixXd &, // w
                     Eigen::MatrixXd &  // A
);

Eigen::MatrixXd BD2A_unref(Eigen::MatrixXd, // B
                           Eigen::MatrixXd, // w
                           Eigen::MatrixXd  // A
);


/* Transformes the Saupe vector to the Saupe tensor */
Eigen::MatrixXd Sv2St(Eigen::MatrixXd);

double S2R(double, double, double);
double dS2dR(double, double, double, double, double, double);

/* Perform the full decomposition */
void SaupeEigenSystems(Eigen::MatrixXd &, // B
                       Eigen::MatrixXd &, // w
                       Eigen::MatrixXd &, // D
                       Eigen::MatrixXd &, // A
                       Eigen::MatrixXd &, // eigenvalues
                       Eigen::MatrixXd &, // Euler angles
                       Eigen::MatrixXd &, // eigenvectors
                       bool               // use matrix formalism
);


/* Return only eigen values, eigen vectors & positive euler angles from Saupe
 * matrix */
Eigen::Vector3d Tensor2Euler(Eigen::MatrixXd &, // A
                             Eigen::MatrixXd &, // eigenvalues
                             Eigen::MatrixXd &, // eigenvectors
                             const int,         // permutation
                             const bool sort    = true,
                             const bool positiv = true);

void EigVecs2Quat(Eigen::Matrix3d &, double &, double &, double &, double &);


Eigen::MatrixXcd S2F(Eigen::MatrixXd &, // eigenvalues
                     Eigen::MatrixXd &  // Euler angles
);


Eigen::MatrixXcd back_calc_Y_ref(Eigen::MatrixXd &,  // weights
                                 Eigen::MatrixXd &,  // rdcs
                                 Eigen::MatrixXcd &, // Fmatrix
                                 bool                // use matrix formalism
);


/* Perform DIDC */
void TolmanApproach(Molecule &, BasicInformation &, Eigen::MatrixXd &);


double calculateDihedralAngle(Atom *, Atom *, Atom *, Atom *, StructureOptions);


Eigen::VectorXd getQfacs(Eigen::MatrixXd &, // B
                         Eigen::MatrixXd &, // D
                         Eigen::MatrixXd &, // D_calc
                         Eigen::MatrixXd,   // Dmax
                         Eigen::MatrixXd,   // w
                         Flags &);


Eigen::VectorXd DDcalc2Q(Eigen::MatrixXd &, // D
                         Eigen::MatrixXd &  // D_calc
);

Eigen::Vector3d Polar2Eigen(double, double);

Eigen::Vector3d Polar2Eigen(SphericalHarmonics *, StructureOptions);

void Eigen2Polar(Eigen::Vector3d &, double &, double &);


void joinVectors(Eigen::VectorXd &, Eigen::VectorXd &, Eigen::MatrixXd &);

void joinVectors(Eigen::MatrixXd &, Eigen::VectorXd &, Eigen::MatrixXd &);


double linear_Regression(Eigen::MatrixXd &, // x
                         Eigen::MatrixXd &, // y
                         double &,          // m
                         double &           // b
);


double mean_Value(Eigen::MatrixXd);

Eigen::VectorXd gradatan2(double, double);

Eigen::VectorXd dtheta_dxi(double, double, double);

Eigen::VectorXd dphi_dxi(double, double, double);

double Chi_square(Eigen::MatrixXd &, Eigen::MatrixXd &, double);

double *axb3d(double *, double *);

/*********************************************************
 *                                                       *
 *    Derivative of determinante (chiral volume) with    *
 *     respect to the individiual atoms spanning the     *
 *                  respective vectors.                  *
 *                                                       *
 *********************************************************/

Eigen::Matrix3i get_S3_permutation_matrix(int);

int get_S3_permuted_index(int, int);

int get_S3_permutation_order(int);

inline double
get_permutation_parity(int order)
{
  return ((order % 2) ? (1.0) : (-1.0));
}

double get_derived_Leibnitz_element(Eigen::Matrix3d &determinante,
                                    int permutation,
                                    int axis_to_derive,
                                    int coordinate_to_derive);

double get_derived_Lebinitz_determinante(Eigen::Matrix3d &determinante,
                                         int axis_to_derive,
                                         int coordinate_to_derive);

#endif

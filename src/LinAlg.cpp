
#include <Atom.hpp>
#include <Declarations.hpp>
#include <LinAlg.hpp>
#include <Molecule.hpp>
#include <RDCdata.hpp>
#include <RDCset.hpp>
#include <SphericalHarmonics.hpp>
#include <Structure.hpp>
#include <cmath>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include <random>

#ifdef USE_CUDA
#include <cuda_mpi_double.hpp>
#endif

/*
 * Factor for interconversion of normalized
 * order tensor and Saupe tensor
 */
const double Sfac = sqrt(3.0) / 2.0;

/*
 * Factor to build the derivative of
 * Wigner rotation element d(2)Mm with
 * respect to angle beta.
 */
const double dd2Mmdb_fac = sqrt(3.0 / 2.0);

/*
 * Standard scaling factors for
 * spherical harmonics.
 * - fac4pi: Needed for scaling
 *    of S^2.
 * - facSQ4pi: Needed for standard
 *    scaling.
 */
const double fac4pi    = (4.0 * PI_ / 5.0);
const double facSQR4pi = sqrt(fac4pi);

/*
 * Scaling factors for F matrix
 */
const double facSQR3o8 = sqrt(3.0 / 8.0);
const double fac2o3    = 2.0 / 3.0;

/*
 * Scaling factor for building B matrix.
 */
const double SQRhalf = sqrt(1.0 / 2.0);

Tensor::Tensor()
{
  T     = Eigen::Matrix3d::Zero(3, 3);
  Evals = Eigen::Vector3d::Zero(NUM_EIGENVALUES_);
  Evecs = Eigen::Matrix3d::Zero(EIGENVECTOR_ELEMENTS_, EIGENVECTOR_ELEMENTS_);
  Euler = Eigen::MatrixXd::Zero(NUM_EULER_PERMUTATIONS_, NUM_EULER_ANGLES_);
}

Tensor::Tensor(Eigen::Matrix3d &Tnew)
{
  setTensor(Tnew);
}

void
Tensor::setTensor(Eigen::Matrix3d &Tnew)
{
  int maxIndex, minIndex;
  maxIndex = minIndex = 0;

  T     = Tnew;
  Euler = Eigen::MatrixXd::Zero(NUM_EULER_PERMUTATIONS_, NUM_EULER_ANGLES_);
  Eigen::EigenSolver<Eigen::MatrixXd> es(T);

  Eigen::MatrixXd tmp_evec = es.eigenvectors().real();
  Eigen::VectorXd tmp_eval = es.eigenvalues().real();

  tmp_eval.cwiseAbs().maxCoeff(&maxIndex);
  tmp_eval.cwiseAbs().minCoeff(&minIndex);

  Evecs.col(2) = tmp_evec.col(maxIndex);
  Evecs.col(1) = tmp_evec.col(3 - minIndex - maxIndex);
  Evecs.col(0) = tmp_evec.col(minIndex);
  Evals(2, 0)  = tmp_eval(maxIndex);
  Evals(1, 0)  = tmp_eval(3 - minIndex - maxIndex);
  Evals(0, 0)  = tmp_eval(minIndex);

  if (Evecs.determinant() < 0.0)
  {
    Evecs.col(0) *= -1.0;
  }
}

int
Eigen2Double(Eigen::MatrixXd &E, double **D)
{
  int R, C;
  R = E.rows();
  C = E.cols();
  return Eigen2Double(E, D, R, C);
}

int
Eigen2Double(Eigen::MatrixXd &E, double **D, const int R, const int C)
{
  int r, c, d;
  for (r = d = 0; r < R; ++r)
  {
    for (c = 0; c < C; ++c, ++d)
    {
      D[0][d] = E(r, c);
    }
  }
  return d;
}

int
Double2Eigen(Eigen::MatrixXd &E, double *D)
{
  int R, C;
  R = E.rows();
  C = E.cols();
  return Double2Eigen(E, D, R, C);
}

int
Double2Eigen(Eigen::MatrixXd &E, double *D, const int R, const int C)
{
  int r, c, d;
  E = Eigen::MatrixXd::Zero(R, C);
  for (r = d = 0; r < R; ++r)
  {
    for (c = 0; c < C; ++c, ++d)
    {
      E(r, c) = D[d];
    }
  }
  return d;
}


double
normRand(double mu, double sig)
{
  static std::default_random_engine generator(time(NULL));
  static std::normal_distribution<double> distribution(0.0, 1.0);
  return ((distribution(generator) * sig) + mu);
}

/*
 * Transform arbitrary polar angles
 * to standard polar angle convention
 * (theta=[0,pi], phi=[-pi,pi]).
 */

void
A2goodA(double &t, double &p)
{
  if (t < .0)
  {
    t = fabs(t);
    p = PI_ + p;
  }
  while (t > (2.0 * PI_))
  {
    t = t - 2.0 * PI_;
  }
  while (t > PI_)
  {
    t = (2.0 * PI_ - t);
    p += PI_;
  }
  while (p > PI_)
  {
    p = -(2.0 * PI_ - p);
  }
  while (p < -PI_)
  {
    p = p + 2.0 * PI_;
  }
}

/*
 * Build the Moore-Penrose inverse
 * (= pseudo inverse) of a matrix.
 * This function is enables for
 * provided for "normal" and
 * complex matrizes.
 */

template<typename X>
Eigen::Matrix<X, -1, -1, 0, -1, -1>
large_MoorePenroseInverse(Eigen::Matrix<X, -1, -1, 0, -1, -1> &Mat,
                          double cutoff)
{
  Eigen::BDCSVD<Eigen::Matrix<X, Eigen::Dynamic, Eigen::Dynamic>> svd(
      Mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
  int i;
  int rows = svd.matrixU().cols();
  int cols = svd.matrixV().rows();
  Eigen::Matrix<X, Eigen::Dynamic, Eigen::Dynamic> SigmaT =
      Eigen::Matrix<X, Eigen::Dynamic, Eigen::Dynamic>::Zero(cols, rows);
  Eigen::Matrix<X, Eigen::Dynamic, 1> s = svd.singularValues();
  for (i = 0; i < rows && i < cols; ++i)
  {
    if (s(i) < cutoff)
      SigmaT(i, i) = .0;
    else
      SigmaT(i, i) = 1.0 / s(i);
  }
  Eigen::Matrix<X, Eigen::Dynamic, Eigen::Dynamic> MPI =
      svd.matrixV() * SigmaT * svd.matrixU().transpose();

  return MPI;
}

template<typename X>
Eigen::Matrix<X, -1, -1, 0, -1, -1>
small_MoorePenroseInverse(Eigen::Matrix<X, -1, -1, 0, -1, -1> &Mat,
                          double cutoff)
{
  Eigen::JacobiSVD<Eigen::Matrix<X, Eigen::Dynamic, Eigen::Dynamic>> svd(
      Mat, Eigen::ComputeFullU | Eigen::ComputeFullV);

  int i;
  int rows = svd.matrixU().cols();
  int cols = svd.matrixV().rows();
  Eigen::Matrix<X, Eigen::Dynamic, Eigen::Dynamic> SigmaT =
      Eigen::Matrix<X, Eigen::Dynamic, Eigen::Dynamic>::Zero(cols, rows);
  Eigen::Matrix<X, Eigen::Dynamic, 1> s = svd.singularValues();
  for (i = 0; i < rows && i < cols; ++i)
  {
    if (s(i) < cutoff)
      SigmaT(i, i) = .0;
    else
      SigmaT(i, i) = 1.0 / s(i);
  }
  Eigen::Matrix<X, Eigen::Dynamic, Eigen::Dynamic> MPI =
      svd.matrixV() * SigmaT * svd.matrixU().transpose();

  return MPI;
}

template<typename X>
Eigen::Matrix<X, -1, -1, 0, -1, -1>
MoorePenroseInverse(Eigen::Matrix<X, -1, -1, 0, -1, -1> &Mat,
                    double cutoff,
                    bool enable_gpu)
{
#ifdef USE_CUDA
  Eigen::MatrixXd inverse;
  if (enable_gpu && Mat.rows() > LARGE_MATRIX_LIMIT_)
  {
    inverse = GPU_MoorePenroseInverse(Mat, cutoff);
    if (inverse.cols() == 1 && inverse.rows() == 1)
      inverse = small_MoorePenroseInverse<X>(Mat, cutoff);
  }
  else
    inverse = small_MoorePenroseInverse<X>(Mat, cutoff);
  return inverse;
#else
  return small_MoorePenroseInverse<X>(Mat, cutoff);
#endif
}

template Eigen::Matrix<double, -1, -1, 0, -1, -1>
MoorePenroseInverse(Eigen::Matrix<double, -1, -1, 0, -1, -1> &,
                    double cutoff,
                    bool enable_gpu);
// template Eigen::Matrix<float, -1, -1, 0, -1, -1> MoorePenroseInverse (
// Eigen::Matrix<float, -1, -1, 0, -1, -1>&, double cutoff );

#ifdef USE_CUDA
Eigen::MatrixXd
GPU_MoorePenroseInverse(Eigen::MatrixXd &A, double cutoff)
{
  const double m   = A.rows();
  const double n   = A.cols();
  const double lda = m;
  double *MPI_A; // = (double**) malloc ( sizeof(double*) );
  Eigen::MatrixXd inverse = Eigen::MatrixXd::Identity(1, 1);

  int cuda_code = cuda_double_moore_penrose_inverse_host_memory(
      A.data(), &MPI_A, m, n, lda, cutoff);
  if (cuda_code)
  {
    std::cerr << "cuda_double_moore_penrose_inverse exited with " << cuda_code
              << std::endl;
    return inverse;
  }
  inverse = Eigen::Map<Eigen::MatrixXd>(MPI_A, n, m);
  free(MPI_A);
  return inverse;
}
#endif

Eigen::MatrixXcd
large_MoorePenroseInverse(Eigen::MatrixXcd &M, double cutoff)
{
  Eigen::BDCSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullU |
                                             Eigen::ComputeFullV);

  std::complex<double> const n(.0, .0);
  int i;
  int rows                = svd.matrixU().cols();
  int cols                = svd.matrixV().rows();
  Eigen::MatrixXcd SigmaT = Eigen::MatrixXcd::Zero(cols, rows);
  Eigen::VectorXcd s      = svd.singularValues();
  for (i = 0; i < rows && i < cols; ++i)
  {
    if ((s(i).imag() < cutoff) && (s(i).real() < cutoff))
      SigmaT(i, i) = n;
    else
      SigmaT(i, i) = 1.0 / s(i);
  }
  Eigen::MatrixXcd MPI = svd.matrixV() * SigmaT * svd.matrixU().adjoint();

  return MPI;
}

Eigen::MatrixXcd
small_MoorePenroseInverse(Eigen::MatrixXcd &M, double cutoff)
{
  Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullU |
                                                Eigen::ComputeFullV);

  std::complex<double> const n(.0, .0);
  int i;
  int rows                = svd.matrixU().cols();
  int cols                = svd.matrixV().rows();
  Eigen::MatrixXcd SigmaT = Eigen::MatrixXcd::Zero(cols, rows);
  Eigen::VectorXcd s      = svd.singularValues();
  for (i = 0; i < rows && i < cols; ++i)
  {
    if ((s(i).imag() < cutoff) && (s(i).real() < cutoff))
      SigmaT(i, i) = n;
    else
      SigmaT(i, i) = 1.0 / s(i);
  }
  Eigen::MatrixXcd MPI = svd.matrixV() * SigmaT * svd.matrixU().adjoint();

  return MPI;
}

Eigen::MatrixXcd
MoorePenroseInverse(Eigen::MatrixXcd &M, double cutoff)
{
  Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullU |
                                                Eigen::ComputeFullV);
  if (M.cols() > LARGE_MATRIX_LIMIT_)
    return large_MoorePenroseInverse(M, cutoff);
  else
    return small_MoorePenroseInverse(M, cutoff);
}

Eigen::MatrixXcd
get_svd_U(Eigen::MatrixXcd &M, bool full)
{
  if (full)
  {
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullU);
    return svd.matrixU();
  }
  else
  {
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeThinU);
    return svd.matrixU();
  }
}

Eigen::MatrixXd
get_svd_U(Eigen::MatrixXd &M, bool full)
{
  if (full)
  {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullU);
    return svd.matrixU();
  }
  else
  {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU |
                                                 Eigen::ComputeThinV);
    return svd.matrixU();
  }
}

/*
 * Scales and normalizes the provided
 * RDC matrix to Dmax and a normalization
 * factor.
 */

void
ScaleRDCMatrix(Eigen::MatrixXd &rdcMatrix,
               Eigen::MatrixXd &rdcMatrixScaled,
               Eigen::MatrixXd &rdcMatrixNorm,
               Molecule &CurrMol,
               struct BasicInformation &baseInformation,
               StructureOptions options)
{
  unsigned int i, j;
  RDCdata *CurrRDC;
  RDCset *CurrSet;
  if (rdcMatrixScaled.cols() != rdcMatrix.cols() ||
      rdcMatrixScaled.rows() != rdcMatrix.rows())
  {
    rdcMatrixScaled.resize(rdcMatrix.rows(), rdcMatrix.cols());
  }
  if (rdcMatrixNorm.cols() != rdcMatrix.cols() ||
      rdcMatrixNorm.rows() != rdcMatrix.rows())
  {
    rdcMatrixNorm.resize(rdcMatrix.rows(), rdcMatrix.cols());
  }
  for (i = 0, CurrSet = CurrMol.getHeadSet(); CurrSet;
       CurrSet = CurrSet->getNext(), ++i)
  {
    for (j = 0, CurrRDC = CurrSet->getHeadData(); CurrRDC;
         CurrRDC = CurrRDC->getNext(), ++j)
    {
      CurrRDC->setKappa();

      rdcMatrixScaled(j, i) =
          rdcMatrix(j, i) /
          (CurrRDC->getKappa() / pow(CurrRDC->getDistance(options), 3.0));
      rdcMatrixNorm(j, i) = rdcMatrixScaled(j, i) * baseInformation.normFactRDC;

      CurrRDC->setD_scaled(rdcMatrixScaled(j, i));
      CurrRDC->setD_norm(rdcMatrixNorm(j, i));
    }
  }
}

std::complex<double>
Ylm(const int l, const int m, const double theta, const double phi)
{
  std::complex<double> Y(.0, .0);
  std::complex<double> e(.0, .0);
  std::complex<double> P(.0, .0);
  std::complex<double> const i(.0, 1.0);

  double efac, prefac;
  int pm = abs(m);

  P.real(std::assoc_legendre(l, pm, cos(theta)) * (1 & m ? -1.0 : 1.0));

  efac = ((double) pm) * phi;

  e = std::exp(i * efac);

  prefac = ((double) (2 * l + 1) / 2.0) * (fac(l - pm) / fac(l + pm));
  prefac = sqrt(prefac) / sqrt(2.0 * PI_);

  Y = prefac * P * e;
  if (m < 0)
  {
    Y = pow(-1.0, pm) * std::conj(Y);
  }
  return Y;
}

/*
 * Builds the factorial of an int
 * number and returns the respective
 * double.
 */

double
fac(const int a)
{
  if (a < 0)
  {
    std::cerr << "ERROR:\tFactorial of Number < 0 was requested...\n\tResult 0 "
                 "will be returned...\n";
    return .0;
  }
  if (a <= 1)
    return 1.0;

  return fac(a - 1) * (double) a;
}

/*
 * Builds a row of the normalized
 * cosine matrix B based on spherical
 * angles theta and phi.
 *
 * If needed the cosine row can be
 * rotated based on euler angles alpha,
 * beta and gamma (ZYZ) convention.
 */

Eigen::MatrixXd
Sphere2B(const double theta,
         const double phi,
         const double alpha,
         const double beta,
         const double gamma)
{
  Eigen::MatrixXd output(1, COSINE_ELEMENTS_);
  output(0, 0) =
      (facSQR4pi * D2MmY2m(0, theta, phi, alpha, beta, gamma).real());
  output(0, 1) = (SQRhalf * facSQR4pi *
                  (D2MmY2m(2, theta, phi, alpha, beta, gamma).real() +
                   D2MmY2m(-2, theta, phi, alpha, beta, gamma).real()));
  output(0, 2) = -(SQRhalf * facSQR4pi *
                   (D2MmY2m(-2, theta, phi, alpha, beta, gamma).imag() -
                    D2MmY2m(2, theta, phi, alpha, beta, gamma).imag()));
  output(0, 3) = (SQRhalf * facSQR4pi *
                  (D2MmY2m(-1, theta, phi, alpha, beta, gamma).real() -
                   D2MmY2m(1, theta, phi, alpha, beta, gamma).real()));
  output(0, 4) = -(SQRhalf * facSQR4pi *
                   (D2MmY2m(1, theta, phi, alpha, beta, gamma).imag() +
                    D2MmY2m(-1, theta, phi, alpha, beta, gamma).imag()));

  return output;
}

/*
 * Builds the small d(2)Mm for Wigner
 * rotations. This part just includes
 * the the rotational part beta. For
 * more information see:
 *
 * M. Eden, Computer simulations in
 * solid-state NMR. I. Spin dynamics
 * theory, Concepts in Magnetic Resonance,
 * 17a (1), 2003, 117-154.
 * https://doi.org/10.1002/cmr.a.10061
 */

double
d2Mm(const int M, const int m, const double beta)
{
  switch (M)
  {
    case 2: {
      switch (m)
      {
        case 2:
          return (0.25 * pow((1.0 + cos(beta)), 2.0));
        case 1:
          return (-0.5 * sin(beta) * (1.0 + cos(beta)));
        case 0:
          return (facSQR3o8 * pow(sin(beta), 2.0));
        case -1:
          return (-0.5 * sin(beta) * (1.0 - cos(beta)));
        case -2:
          return (0.25 * pow((1.0 - cos(beta)), 2.0));
        default:
          std::cerr
              << "ERROR:\tThe small Wigner rotation elment d" << M << "," << m
              << "could not be build...\n"
              << "\t\tTITANIA does just support values for m = [-2,2]...\n";
          return .0;
      };
    };
    break;

    case 1: {
      switch (m)
      {
        case 2:
          return (0.5 * sin(beta) * (1.0 + cos(beta)));
        case 1:
          return (0.5 * (2.0 * pow(cos(beta), 2.0) + cos(beta) - 1.0));
        case 0:
          return (-facSQR3o8 * sin(2.0 * beta));
        case -1:
          return (-0.5 * (2.0 * pow(cos(beta), 2.0) - cos(beta) -
                          1.0)); /* Check for correctness */
        case -2:
          return (-0.5 * sin(beta) * (1.0 - cos(beta)));
        default:
          std::cerr
              << "ERROR:\tThe small Wigner rotation elment d" << M << "," << m
              << "could not be build...\n"
              << "\t\tTITANIA does just support values for m = [-2,2]...\n";
          return .0;
      };
    };
    break;

    case 0: {
      switch (m)
      {
        case 2:
          return (facSQR3o8 * pow(sin(beta), 2.0));
        case 1:
          return (facSQR3o8 * sin(2.0 * beta));
        case 0:
          return (0.5 * (3.0 * pow(cos(beta), 2.0) - 1.0));
        case -1:
          return (-facSQR3o8 * sin(2.0 * beta));
        case -2:
          return (facSQR3o8 * pow(sin(beta), 2.0));
        default:
          std::cerr
              << "ERROR:\tThe small Wigner rotation elment d" << M << "," << m
              << "could not be build...\n"
              << "\t\tTITANIA does just support values for m = [-2,2]...\n";
          return .0;
      };
    };
    break;

    case -1: {
      switch (m)
      {
        case 2:
          return (0.5 * sin(beta) * (1.0 - cos(beta)));
        case 1:
          return (-0.5 * (2.0 * pow(cos(beta), 2.0) - cos(beta) - 1.0));
        case 0:
          return (facSQR3o8 * sin(2.0 * beta));
        case -1:
          return (0.5 * (2.0 * pow(cos(beta), 2.0) + cos(beta) - 1.0));
        case -2:
          return (-0.5 * sin(beta) * (1.0 + cos(beta)));
        default:
          std::cerr
              << "ERROR:\tThe small Wigner rotation elment d" << M << "," << m
              << "could not be build...\n"
              << "\t\tTITANIA does just support values for m = [-2,2]...\n";
          return .0;
      };
    };
    break;


    case -2: {
      switch (m)
      {
        case 2:
          return (0.25 * pow((1.0 - cos(beta)), 2.0));
        case 1:
          return (0.5 * sin(beta) * (1.0 - cos(beta)));
        case 0:
          return (facSQR3o8 * pow(sin(beta), 2.0));
        case -1:
          return (0.5 * sin(beta) * (1.0 + cos(beta)));
        case -2:
          return (0.25 * pow((1.0 + cos(beta)), 2.0));
        default:
          std::cerr
              << "ERROR:\tThe small Wigner rotation elment d" << M << "," << m
              << "could not be build...\n"
              << "\t\tTITANIA does just support values for m = [-2,2]...\n";
          return .0;
      };
    };
    break;

    default:
      std::cerr << "ERROR:\tThe small Wigner rotation elment d" << M << "," << m
                << "could not be build...\n"
                << "\t\tTITANIA does just support values for M = [-2,2]...\n";
      return .0;
  }
}

/*
 * Builds the respective derivative
 * of the small Wigner elements d(2)Mm
 * with respect to the angle beta.
 */

double
dd2Mmdb(const int M, const int m, const double beta)
{
  switch (M)
  {
    case 2: {
      switch (m)
      {
        case 2:
          return (d2Mm(M, (m - 1), beta));
        case 1:
          return (-0.5 * (cos(beta) + cos(2.0 * beta)));
        case 0:
          return (dd2Mmdb_fac * sin(beta) * cos(beta));
        case -1:
          return (-0.5 * (cos(beta) - cos(2.0 * beta)));
        case -2:
          return (-d2Mm(M, (m + 1), beta));
        default:
          std::cerr
              << "ERROR:\tThe derivative small Wigner rotation elment d" << M
              << "," << m << "could not be build...\n"
              << "\t\tTITANIA does just support values for m = [-2,2]...\n";
          return .0;
      };
    };
    break;

    case 1: {
      switch (m)
      {
        case 2:
          return (-dd2Mmdb(m, M, beta));
        case 1:
          return (-0.5 * sin(beta) * (1.0 + 4.0 * cos(beta)));
        case 0:
          return (-dd2Mmdb_fac * cos(2.0 * beta));
        case -1:
          return (0.5 * sin(beta) * (1.0 - 4.0 * cos(beta)));
        case -2:
          return (dd2Mmdb(-m, -M, beta));
        default:
          std::cerr
              << "ERROR:\tThe derivative small Wigner rotation elment d" << M
              << "," << m << "could not be build...\n"
              << "\t\tTITANIA does just support values for m = [-2,2]...\n";
          return .0;
      };
    };
    break;

    case 0: {
      switch (m)
      {
        case 2:
          return (dd2Mmdb(m, M, beta));
        case 1:
          return (-dd2Mmdb(m, M, beta));
        case 0:
          return (-3.0 * sin(beta) * cos(beta));
        case -1:
          return (dd2Mmdb(-m, M, beta));
        case -2:
          return (dd2Mmdb(-m, M, beta));
        default:
          std::cerr
              << "ERROR:\tThe derivative small Wigner rotation elment d" << M
              << "," << m << "could not be build...\n"
              << "\t\tTITANIA does just support values for m = [-2,2]...\n";
          return .0;
      };
    };
    break;

    case -1: {
      switch (m)
      {
        case 2:
          return (-dd2Mmdb(m, M, beta));
        case 1:
          return (-dd2Mmdb(m, M, beta));
        case 0:
          return (-dd2Mmdb(m, M, beta));
        case -1:
          return (-dd2Mmdb(-m, -M, beta));
        case -2:
          return (dd2Mmdb(-m, -M, beta));
        default:
          std::cerr
              << "ERROR:\tThe derivative small Wigner rotation elment d" << M
              << "," << m << "could not be build...\n"
              << "\t\tTITANIA does just support values for m = [-2,2]...\n";
          return .0;
      };
    };
    break;


    case -2: {
      switch (m)
      {
        case 2:
          return (dd2Mmdb(m, M, beta));
        case 1:
          return (-dd2Mmdb(-M, -m, beta));
        case 0:
          return (dd2Mmdb(-M, m, beta));
        case -1:
          return (-dd2Mmdb(-M, -m, beta));
        case -2:
          return (dd2Mmdb(-M, -m, beta));
        default:
          std::cerr
              << "ERROR:\tThe derivative small Wigner rotation elment d" << M
              << "," << m << "could not be build...\n"
              << "\t\tTITANIA does just support values for m = [-2,2]...\n";
          return .0;
      };
    };
    break;

    default:
      std::cerr << "ERROR:\tThe derivative small Wigner rotation elment d" << M
                << "," << m << "could not be build...\n"
                << "\t\tTITANIA does just support values for M = [-2,2]...\n";
      return .0;
  }
}

/*
 * Builds the full Wigner rotation element
 * D2Mm for the sperical harmonics.
 */

std::complex<double>
D2Mm(const int M,
     const int m,
     const double alpha,
     const double beta,
     const double gamma)
{
  std::complex<double> result(.0, .0);
  const std::complex<double> i(.0, 1.0);
  result = (std::exp(-i * ((double) M) * alpha) * d2Mm(M, m, beta) *
            std::exp(-i * ((double) m) * gamma));
  return result;
}

/*
 * Builds the derivative of the full
 * Wigner rotation element D(2)Mm
 * with respect to beta.
 */

std::complex<double>
DeltaD2MmDeltaBeta(const int M,
                   const int m,
                   const double alpha,
                   const double beta,
                   const double gamma)
{
  const std::complex<double> i(.0, 1.0);
  const std::complex<double> efac =
      std::exp(-i * ((double) M) * alpha) * std::exp(-i * ((double) m) * gamma);

  if (M == (-2))
  {
    return (2.0 * sqrt(3.0 / 8.0) * cos(beta) * sin(beta) * efac);
  }
  if (M == (-1))
  {
    return (2.0 * sqrt(3.0 / 8.0) * cos(2.0 * beta) * efac);
  }
  if (M == 0)
  {
    return (-3.0 * cos(beta) * sin(beta) * efac);
  }
  if (M == 1)
  {
    return (-2.0 * sqrt(3.0 / 8.0) * cos(2.0 * beta) * efac);
  }
  if (M == 2)
  {
    return (2.0 * sqrt(3.0 / 8.0) * cos(beta) * sin(beta) * efac);
  }
  return .0;
}

/*
 * Builds the rotated Y2m for the
 * specific polar/Euler angles theta,
 * phi, alpha, beta and gamma.
 */

std::complex<double>
D2MmY2m(const int m,
        const double theta,
        const double phi,
        const double alpha,
        const double beta,
        const double gamma)
{
  std::complex<double> Y2m(.0, .0);
  int M;
  for (M = -2; M < 3; ++M)
  {
    Y2m += (D2Mm(M, m, alpha, beta, gamma) * Ylm(2, M, theta, phi));
  }

  return Y2m;
}

/*
 * Performs the full analyis for
 * refined spherical harmonics of
 * a structure. The values for
 * S^2(RDC), S^2(axial), eta
 * and phi(aniso) are computed.
 *
 * For more information check:
 *
 * J. Meiler et. al., JACS, 2001, 123, 6098-6107.
 * DOI: 10.1021/ja010002z
 *
 * W. Peti et. al., JACS, 2002, 124, 5822-5833.
 * DOI: 10.1021/ja011883c
 *
 *
 * S^2(axial), also lambda(rdc): see Meiler.
 */

void
motionalAnalysisOfSphericalHarmonics(Eigen::MatrixXcd Y,
                                     double thetaAv,
                                     double phiAv,
                                     double &sRDC,
                                     double &sAx,
                                     double &etaRDC,
                                     double &phiRDC)
{
  Eigen::MatrixXcd Y_VF = Eigen::MatrixXcd::Zero(1, HARMONIC_ELEMENTS_);
  std::complex<double> etaD, etaN, Sc;
  int m, i;
  for (m = -2; m < 3; ++m)
  {
    for (int M = -2; M < 3; ++M)
    {
      Y_VF(0, (m + 2)) += (D2Mm(M, m, phiAv, thetaAv, 0) * Y(0, (M + 2)));
    }
  }

  etaD = etaN = Sc = .0;
  etaN             = 2.0 * Y_VF(0, 0) * Y_VF(0, 4);

  for (i = 0; i < HARMONIC_ELEMENTS_; ++i)
  {
    Sc += (Y_VF(0, i) * std::conj(Y_VF(0, i)));
    etaD += (Y_VF(0, i) * Y_VF(0, 4 - i));
  }

  etaRDC = sqrt(etaN / etaD).real();
  phiRDC = atan2(Y_VF(0, 0).imag(), Y_VF(0, 0).real());
  sAx    = fac4pi * (Y_VF(0, 2) * Y_VF(0, 2)).real();
  sRDC   = fac4pi * Sc.real();
}

/*
 * Computes the S(overall) for a structure.
 * This is basically just a scaling of all
 * S^2(rdc) to obtain values <= 1.0
 *
 * For more information check:
 *
 * J. Meiler et. al., JACS, 2001, 123, 6098-6107.
 * DOI: 10.1021/ja010002z
 *
 * W. Peti et. al., JACS, 2002, 124, 5822-5833.
 * DOI: 10.1021/ja011883c
 *
 *
 * S^2(axial), also lambda(rdc): see Meiler.
 */

double
SrdcBySoverall(Eigen::MatrixXd Srdc)
{
  double maxSrdc = -1.0;
  for (int i = 0; i < Srdc.cols(); ++i)
  {
    if (Srdc(0, i) > maxSrdc)
      maxSrdc = Srdc(0, i);
  }
  maxSrdc = sqrt(1.0 / maxSrdc);
  return maxSrdc;
}

/*
 * Calculates the alignment tensor for a D matrix,
 * cosine matrix and free to choose weighting.
 */

Eigen::MatrixXd
BD2A(Eigen::MatrixXd &B, Eigen::MatrixXd &w, Eigen::MatrixXd &D)
{
  Eigen::MatrixXd C     = w * B;
  Eigen::MatrixXd MPI_B = MoorePenroseInverse<double>(C);

  /* TODO
   * Previously a Nullspace expression analoug to Tolman.
   * This lead to differences for matrix and vector formalism of the MFA
   * approach. Maybe this should be a flag? Eigen::MatrixXd Nullspace =  (
   * Eigen::MatrixXd::Identity(SAUPE_ELEMENTS_, SAUPE_ELEMENTS_) - MPI_B * B ) *
   * Eigen::MatrixXd::Random(SAUPE_ELEMENTS_,D.cols());
   */
  Eigen::MatrixXd A_bf = MPI_B * w * D; // + Nullspace;

  return MPI_B * w * D;
}

/*
 * At some points a version of BD2A witch uses a
 * call by value (not call by reference) is needed.
 * Therefor a second function has to be declared.
 * Overloading is not allowed due to the fact, that
 * the compiler is unable to determine which function
 * has to be used.
 */

Eigen::MatrixXd
BD2A_unref(Eigen::MatrixXd B, Eigen::MatrixXd w, Eigen::MatrixXd D)
{
  return BD2A(B, w, D);
}

/*
 * Transformas the normalized order tensor
 * to the Saupe tensor.
 */

Eigen::MatrixXd
Sv2St(Eigen::MatrixXd A)
{
  Eigen::MatrixXd S =
      Eigen::MatrixXd::Zero(SAUPE_TENSOR_DIM_, SAUPE_TENSOR_DIM_);
  const int i = 0;

  S(STANDARD_AXIS_X_, STANDARD_AXIS_X_) =
      (2.0 * Sfac * A(SAUPE_XX_YY_, i) - A(SAUPE_ZZ_, i)) / 2.0;
  S(STANDARD_AXIS_Y_, STANDARD_AXIS_Y_) =
      -S(STANDARD_AXIS_X_, STANDARD_AXIS_X_) - A(SAUPE_ZZ_, i);
  S(STANDARD_AXIS_Z_, STANDARD_AXIS_Z_) = A(SAUPE_ZZ_, i);
  S(STANDARD_AXIS_X_, STANDARD_AXIS_Y_) =
      S(STANDARD_AXIS_Y_, STANDARD_AXIS_X_) = Sfac * A(SAUPE_XY_, i);
  S(STANDARD_AXIS_X_, STANDARD_AXIS_Z_) =
      S(STANDARD_AXIS_Z_, STANDARD_AXIS_X_) = Sfac * A(SAUPE_XZ_, i);
  S(STANDARD_AXIS_Y_, STANDARD_AXIS_Z_) =
      S(STANDARD_AXIS_Z_, STANDARD_AXIS_Y_) = Sfac * A(SAUPE_YZ_, i);
  return S;
}

double
S2R(double Sxx, double Syy, double Szz)
{
  return (fac2o3 * (Sxx - Syy) / Szz);
}

double
dS2dR(double Sxx, double Syy, double Szz, double dSxx, double dSyy, double dSzz)
{
  // the signs do not matter due to square!
  double dRdSxx = fac2o3 / Szz;
  double dRdSyy = dRdSxx;
  double dRdSzz = S2R(Sxx, Syy, Szz) / (Szz);
  double dSdR   = pow(dRdSxx * dSxx, 2.0) + pow(dRdSyy * dSyy, 2.0) +
                pow(dRdSzz * dSzz, 2.0);
  return sqrt(dSdR);
}

/*
 * Solves the full Saupe problem.
 *
 * First the vector representation of all
 * order tensors is calculated.
 * For this the function uses the normalized
 * cosine tensor B, a weighting matrix w and
 * the scaled (to Dmax and the last known Szz)
 * rdcs D.
 *
 * The resulting order matrizes are transformed
 * to the standard Saupe tensor which is
 * decomposed using the eigenvalue decomposition.
 * All resulting eigenvalues/-vectors are sorted
 * to the standard scheme (|Szz| > |Syy| >= |Sxx|).
 */

void
SaupeEigenSystems(Eigen::MatrixXd &B,
                  Eigen::MatrixXd &weights,
                  Eigen::MatrixXd &D,
                  Eigen::MatrixXd &orderMatrix,
                  Eigen::MatrixXd &Eval,
                  Eigen::MatrixXd &Eang,
                  Eigen::MatrixXd &Evec,
                  bool useMatrixFormalism)
{
  unsigned int NOS = D.cols(), NOR = D.rows(), i, j;
  Eigen::MatrixXd A;
  Eigen::MatrixXd tmpEVal(NUM_EIGENVALUES_, 1);
  Eigen::MatrixXd tmpEVec(EIGENVECTOR_ELEMENTS_, EIGENVECTOR_ELEMENTS_);
  Eigen::Vector3d ea;
  Eigen::MatrixXd w;
  if (useMatrixFormalism)
    orderMatrix = BD2A(B, weights, D);
  else
  {
    orderMatrix = Eigen::MatrixXd::Zero(SAUPE_ELEMENTS_, NOS);
    w           = Eigen::MatrixXd::Identity(NOR, NOR);
  }
  for (i = 0; i < NOS; ++i)
  {
    if (!useMatrixFormalism)
    {
      Eigen::MatrixXd D_tmp = D.col(i);
      for (j = 0; j < NOR; ++j)
        w(j, j) = weights(j, i);
      orderMatrix.col(i) = BD2A(B, w, D_tmp);
    }

    A = Sv2St(orderMatrix.col(i));

    ea          = Tensor2Euler(A, tmpEVal, tmpEVec, 0, true, true);
    Eval.col(i) = tmpEVal;
    for (j = 0; j < NUM_EIGENVECTORS_; ++j)
      Evec(j, i) =
          tmpEVec((j % EIGENVECTOR_ELEMENTS_), (j / EIGENVECTOR_ELEMENTS_));
    for (j = 0; j < NUM_EULER_ANGLES_; ++j)
      Eang(j, i) = -ea(j);
  }
}

/*
 * Performes the eigenvalue decomposition
 * of an arbitrary tensor, sors the results
 * (if requested) and returns the (positiv)
 * euler angles.
 */

Eigen::Vector3d
Tensor2Euler(Eigen::MatrixXd &A,
             Eigen::MatrixXd &eigenVal,
             Eigen::MatrixXd &eigenVec,
             const int perm,
             const bool sort,
             const bool positiv)
{
  int minIndex, maxIndex;
  minIndex = maxIndex = 0;
  Eigen::EigenSolver<Eigen::MatrixXd> es(A);
  Eigen::VectorXd ev(NUM_EIGENVALUES_), tmpV(NUM_EIGENVALUES_);
  Eigen::Matrix3d evecs;
  Eigen::Vector3d ea;
  Eigen::MatrixXd tmpM(EIGENVECTOR_ELEMENTS_, EIGENVECTOR_ELEMENTS_);

  evecs = es.eigenvectors().real();
  ev    = es.eigenvalues().real();

  if (sort)
  {
    if (perm == 1)
    {
      evecs.col(STANDARD_AXIS_X_) = -evecs.col(STANDARD_AXIS_X_);
      evecs.col(STANDARD_AXIS_Y_) = -evecs.col(STANDARD_AXIS_Y_);
    }
    else if (perm == 2)
    {
      evecs.col(STANDARD_AXIS_X_) = -evecs.col(STANDARD_AXIS_X_);
      evecs.col(STANDARD_AXIS_Z_) = -evecs.col(STANDARD_AXIS_Z_);
    }
    else if (perm == 3)
    {
      evecs.col(STANDARD_AXIS_Y_) = -evecs.col(STANDARD_AXIS_Y_);
      evecs.col(STANDARD_AXIS_Z_) = -evecs.col(STANDARD_AXIS_Z_);
    }
    ev.cwiseAbs().maxCoeff(&maxIndex);
    ev.cwiseAbs().minCoeff(&minIndex);

    tmpM.col(STANDARD_AXIS_Z_) = evecs.col(maxIndex);
    tmpM.col(STANDARD_AXIS_Y_) =
        evecs.col(EIGENVECTOR_ELEMENTS_ - minIndex - maxIndex);
    tmpM.col(STANDARD_AXIS_X_)    = evecs.col(minIndex);
    eigenVal(STANDARD_AXIS_Z_, 0) = ev(maxIndex);
    eigenVal(STANDARD_AXIS_Y_, 0) = ev(NUM_EIGENVALUES_ - minIndex - maxIndex);
    eigenVal(STANDARD_AXIS_X_, 0) = ev(minIndex);
  }

  for (int i = 0; i < NUM_EIGENVECTORS_; ++i)
    if (fabs(tmpM((i / EIGENVECTOR_ELEMENTS_), (i % EIGENVECTOR_ELEMENTS_))) <
        1e-10)
      tmpM((i / EIGENVECTOR_ELEMENTS_), (i % EIGENVECTOR_ELEMENTS_)) = .0;

  if (tmpM.determinant() < 0.0)
  {
    tmpM.col(STANDARD_AXIS_X_) = -tmpM.col(STANDARD_AXIS_X_);
  }

  if (positiv)
    evecs = tmpM.transpose();
  else
    evecs = tmpM;

  ea = evecs.eulerAngles(STANDARD_AXIS_Z_, STANDARD_AXIS_Y_, STANDARD_AXIS_Z_);

  if (positiv &&
      (ea(EULER_ALPHA_) < .0 || ea(EULER_BETA_) < .0 || ea(EULER_GAMMA_) < .0))
  {
    return Tensor2Euler(A, eigenVal, eigenVec, perm + 1, sort, positiv);
  }

  eigenVec = tmpM;
  return ea;
}

void
EigVecs2Quat(Eigen::Matrix3d &EigVec,
             double &x,
             double &y,
             double &z,
             double &w)
{
  Eigen::Quaterniond q(EigVec);
  x = q.x();
  y = q.y();
  z = q.z();
  w = q.w();
}

/*
 * Builds the Wigner rotation matrix F
 * using the eigenvalues and Euler angles
 * of the Saupe tensors.
 */

Eigen::MatrixXcd
S2F(Eigen::MatrixXd &SaupeEigenValues, Eigen::MatrixXd &SaupeEulerAngles)
{
  const int NOS = SaupeEigenValues.cols();
  int i, j;
  double R;
  std::complex<double> tmp;

  Eigen::MatrixXcd Fmatrix = Eigen::MatrixXcd::Zero(SAUPE_ELEMENTS_, NOS);

  for (i = 0; i < NOS; ++i)
  {
    R = fac2o3 * ((SaupeEigenValues(STANDARD_AXIS_X_, i) -
                   SaupeEigenValues(STANDARD_AXIS_Y_, i)) /
                  SaupeEigenValues(STANDARD_AXIS_Z_, i));
    for (j = 0; j < SAUPE_ELEMENTS_; ++j)
    {
      Fmatrix(j, i) =
          facSQR4pi * (D2Mm(j - 2, 0, SaupeEulerAngles(EULER_GAMMA_, i),
                            SaupeEulerAngles(EULER_BETA_, i),
                            SaupeEulerAngles(EULER_ALPHA_, i)) +
                       facSQR3o8 * R *
                           (D2Mm((j - 2), 2, SaupeEulerAngles(EULER_GAMMA_, i),
                                 SaupeEulerAngles(EULER_BETA_, i),
                                 SaupeEulerAngles(EULER_ALPHA_, i)) +
                            D2Mm((j - 2), -2, SaupeEulerAngles(EULER_GAMMA_, i),
                                 SaupeEulerAngles(EULER_BETA_, i),
                                 SaupeEulerAngles(EULER_ALPHA_, i))));
    }
  }
  return Fmatrix;
}

Eigen::MatrixXcd
back_calc_Y_ref(Eigen::MatrixXd &weights,
                Eigen::MatrixXd &rdcs,
                Eigen::MatrixXcd &Fmatrix,
                bool useMatrixFormalism)
{
  const unsigned int NOR = rdcs.rows(), NOS = rdcs.cols();
  unsigned int r, s;
  Eigen::MatrixXcd Yref = Eigen::MatrixXcd::Zero(NOR, HARMONIC_ELEMENTS_);
  Eigen::MatrixXcd Fmatrix_i, Fweighted;
  Eigen::MatrixXd w, D;

  if (useMatrixFormalism)
  {
    Fweighted             = Fmatrix.transpose();
    Fmatrix_i             = MoorePenroseInverse(Fweighted);
    Eigen::MatrixXcd Yref = (Fmatrix_i * rdcs.transpose()).transpose();
    return Yref;
  }
  else
  {
    w = Eigen::MatrixXd::Zero(NOS, NOS);
    for (r = 0; r < NOR; ++r)
    {
      for (s = 0; s < NOS; ++s)
        w(s, s) = weights(r, s);
      Fweighted   = w * Fmatrix.transpose();
      D           = rdcs.row(r).transpose();
      Fmatrix_i   = MoorePenroseInverse(Fweighted);
      Yref.row(r) = (Fmatrix_i * w * D).transpose();
    }
  }
  return Yref;
}


/*
 * Calculates the sign sensitive dihedral angle
 * of four atoms.
 */


double
calculateDihedralAngle(Atom *A1,
                       Atom *A2,
                       Atom *A3,
                       Atom *A4,
                       enum StructureOptions options)
{
  double angle = .0;

  // Load coordinates
  // All labels according to Bakken, JChemPhys, 2002

  Eigen::Vector3d m = A1->Coordinates2Eigen(options);
  Eigen::Vector3d o = A2->Coordinates2Eigen(options);
  Eigen::Vector3d p = A3->Coordinates2Eigen(options);
  Eigen::Vector3d n = A4->Coordinates2Eigen(options);

  // Calculate bond vectors

  Eigen::Vector3d u = m - o;
  u                 = u / u.norm();
  Eigen::Vector3d v = n - p;
  v                 = v / v.norm();
  Eigen::Vector3d w = p - o;
  w                 = w / w.norm();

  // Calculate | Torsion angle |
  Eigen::Vector3d uXw = u.cross(w);
  uXw                 = uXw / uXw.norm();
  Eigen::Vector3d vXw = v.cross(w);
  vXw                 = vXw / vXw.norm();

  angle = acos(uXw.dot(vXw));

  // Calculate vector which codes for the sign
  Eigen::Vector3d V = vXw.cross(uXw);
  V                 = V / V.norm();
  if (w.dot(V) > .0)
    return -angle;
  else
    return angle;
}

/*
 * Caclulates all Q-factors using D and C matrix.
 */

Eigen::VectorXd
getQfacs(Eigen::MatrixXd &C,
         Eigen::MatrixXd &D,
         Eigen::MatrixXd &D_calc,
         Eigen::MatrixXd Dmax,
         Eigen::MatrixXd w,
         Flags &flags)
{
  const unsigned int NOS = D.cols(), NOR = D.rows();
  unsigned int s, r, i;
  Eigen::MatrixXd A           = Eigen::MatrixXd::Zero(SAUPE_ELEMENTS_, NOS);
  Eigen::MatrixXd EulerAngles = Eigen::MatrixXd::Zero(NUM_EULER_ANGLES_, NOS);
  Eigen::MatrixXd EigenVals   = Eigen::MatrixXd::Zero(NUM_EIGENVALUES_, NOS);
  Eigen::MatrixXd EigenVecs   = Eigen::MatrixXd::Zero(NUM_EIGENVECTORS_, NOS);
  Eigen::MatrixXd Dscaled     = D;
  for (i = 0; i < Dmax.cols(); ++i)
  {
    Dscaled.row(i) /= Dmax(i, i);
  }
  SaupeEigenSystems(C, w, Dscaled, A, EigenVals, EulerAngles, EigenVecs,
                    flags.calculateFullMatrix);
  D_calc = Dcalc(Dmax, C, A);

  if (!flags.calculateFullMatrix)
  {
    for (r = 0; r < NOR; ++r)
    {
      for (s = 0; s < NOS; ++s)
      {
        if (w(r, s) == .0)
          D_calc(r, s) = .0;
      }
    }
  }
  return DDcalc2Q(D, D_calc);
}

/*
 * Calculates all Q-factors using D and Dcalc matrix.
 */

Eigen::VectorXd
DDcalc2Q(Eigen::MatrixXd &D, Eigen::MatrixXd &D_calc)
{
  const unsigned int NOR = D.rows(), NOS = D.cols();
  unsigned int r, s;
  double rezRDCs      = (1.0 / ((double) NOR));
  Eigen::VectorXd Q   = Eigen::VectorXd::Zero(NOS);
  Eigen::VectorXd sum = Q;
  for (r = 0; r < NOR; ++r)
  {
    for (s = 0; s < NOS; ++s)
    {
      Q(s) += (pow((D_calc(r, s) - D(r, s)), 2.0));
      sum(s) += (pow(D(r, s), 2.0));
    }
  }
  Q *= rezRDCs;
  sum *= rezRDCs;
  for (s = 0; s < NOS; ++s)
  {
    Q(s) = (sqrt(Q(s) / sum(s)));
  }
  return Q;
}

/*
 * Builds a unit vector from polar angles.
 */


Eigen::Vector3d
Polar2Eigen(double t, double p)
{
  Eigen::Vector3d V;
  V(STANDARD_AXIS_X_) = sin(t) * cos(p);
  V(STANDARD_AXIS_Y_) = sin(t) * sin(p);
  V(STANDARD_AXIS_Z_) = cos(t);
  return V;
}

Eigen::Vector3d
Polar2Eigen(SphericalHarmonics *Y, StructureOptions opt)
{
  return Polar2Eigen(Y->getTheta(opt), Y->getPhi(opt));
}

void
Eigen2Polar(Eigen::Vector3d &vec, double &theta, double &phi)
{
  Eigen::Vector3d V = vec / vec.norm();
  theta             = acos(V(STANDARD_AXIS_Z_));
  phi               = atan2(V(STANDARD_AXIS_Y_), V(STANDARD_AXIS_X_));
}

void
joinVectors(Eigen::VectorXd &v1, Eigen::VectorXd &v2, Eigen::MatrixXd &out)
{
  unsigned int i;
  Eigen::VectorXd j1 = v1;
  if (v1.cols() != 1)
    j1 = v1.transpose();
  Eigen::VectorXd j2 = v2;
  if (v2.cols() != 1)
    j2 = v2.transpose();

  out = Eigen::MatrixXd::Zero(2, v1.rows());
  for (i = 0; i < j1.rows(); ++i)
  {
    out(0, i) = j1(i);
    out(1, i) = j2(i);
  }
}

void
joinVectors(Eigen::MatrixXd &m1, Eigen::VectorXd &v2, Eigen::MatrixXd &out)
{
  unsigned int j, i;
  Eigen::VectorXd j2 = v2;
  if (v2.cols() != 1)
    j2 = v2.transpose();
  Eigen::MatrixXd j1 = m1;
  if (m1.rows() != j2.rows())
    m1 = j1.transpose();
  out = Eigen::MatrixXd::Zero(j1.cols(), j2.rows());
  for (i = 0; i < j1.rows(); ++i)
  {
    for (j = 0; j < j1.cols(); ++j)
      out(j, i) = j1(j, i);
    out(j, i) = j2(i);
  }
}

double
linear_Regression(Eigen::MatrixXd &x, Eigen::MatrixXd &y, double &m, double &b)
{
  const unsigned int r = x.rows(), c = x.cols();
  unsigned int i, j;

  Eigen::MatrixXd x_xB  = Eigen::MatrixXd::Zero(r, c);
  Eigen::MatrixXd y_yB  = Eigen::MatrixXd::Zero(r, c);
  Eigen::MatrixXd yR_yB = Eigen::MatrixXd::Zero(r, c);

  double xB, yB, xxB_yyB, xxB_xxB, yyB_yyB, yRyB_yRyB, R;

  xB      = mean_Value(x);
  yB      = mean_Value(y);
  xxB_yyB = xxB_xxB = yRyB_yRyB = yyB_yyB = .0;
  for (i = 0; i < r; ++i)
  {
    for (j = 0; j < c; ++j)
    {
      x_xB(i, j) = x(i, j) - xB;
      y_yB(i, j) = y(i, j) - yB;
      xxB_yyB += (x_xB(i, j) * y_yB(i, j));
      xxB_xxB += (x_xB(i, j) * x_xB(i, j));
      yyB_yyB += (y_yB(i, j) * y_yB(i, j));
    }
  }
  m = (xxB_yyB / xxB_xxB);
  b = yB - m * xB;

  for (i = 0; i < r; ++i)
  {
    for (j = 0; j < c; ++j)
    {
      yR_yB(i, j) = (m * x(i, j) + b - yB);
      yRyB_yRyB += (yR_yB(i, j) * yR_yB(i, j));
    }
  }
  R = (yRyB_yRyB / yyB_yyB);
  return R;
}

double
mean_Value(Eigen::MatrixXd M)
{
  double M_bar = .0, N = .0;
  unsigned int i, j;
  for (i = 0; i < M.cols(); ++i)
  {
    for (j = 0; j < M.rows(); ++j)
    {
      M_bar += M(j, i);
      N += 1.0;
    }
  }
  return (M_bar / N);
}

/***************************************
 *            gradatan2                *
 *                                     *
 * grad(0) = d(atan2)/dx               *
 * grad(1) = d(atan2)/dy               *
 *                                     *
 ***************************************/

Eigen::VectorXd
gradatan2(double y, double x)
{
  Eigen::VectorXd grad = Eigen::VectorXd::Zero(2);

  double denom = x * x + y * y;

  grad(0) = -y / denom;
  grad(1) = x / denom;

  return grad;
}

Eigen::VectorXd
dtheta_dxi(double x, double y, double z)
{
  double R             = sqrt(x * x + y * y + z * z);
  double R3            = pow(R, 3.0);
  Eigen::VectorXd dtdx = Eigen::VectorXd::Zero(NUMBER_OF_AXIS_);
  // d(acos(u))/d(u(x_i))
  double dadu            = -(1.0 / sqrt(1.0 - pow((z / R), 2.0)));
  dtdx(STANDARD_AXIS_X_) = dadu * (-z * x / R3);
  dtdx(STANDARD_AXIS_Y_) = dadu * (-z * y / R3);
  dtdx(STANDARD_AXIS_Z_) = dadu * ((R - z * z / R) / pow(R, 2.0));
  return dtdx;
}

Eigen::VectorXd
dphi_dxi(double x, double y, double z)
{
  double R             = sqrt(x * x + y * y + z * z);
  double R2            = pow(R, 2.0);
  double R3            = pow(R, 3.0);
  Eigen::VectorXd dpdx = Eigen::VectorXd::Zero(NUMBER_OF_AXIS_);
  // dp/dx = gradatan2(0)
  // dp/dy = gradatan2(1)
  Eigen::VectorXd datan2 = gradatan2((y / R), (x / R));
  dpdx(STANDARD_AXIS_X_) =
      datan2(1) * (-y * x / R3) + datan2(0) * ((R - x * x / R) / R2);
  dpdx(STANDARD_AXIS_Y_) =
      datan2(1) * ((R - y * y / R) / R2) + datan2(0) * (-x * y / R3);
  dpdx(STANDARD_AXIS_Z_) =
      datan2(1) * (-y * z / R3) + datan2(0) * (-x * z / R3);
  return dpdx;
}

double
Chi_square(Eigen::MatrixXd &x1, Eigen::MatrixXd &x2, double sigma_square)
{
  double chi_2 = .0;
  int r, c, R, C;
  R = x1.rows();
  C = x1.cols();
  for (c = 0; c < C; ++c)
  {
    for (r = 0; r < R; ++r)
    {
      chi_2 += pow((x1(r, c) - x2(r, c)), 2.0);
    }
  }
  return (chi_2 / sigma_square);
}

double *
axb3d(double *a, double *b)
{
  double *result = (double *) malloc(3 * sizeof(double));

  result[STANDARD_AXIS_X_] = (a[STANDARD_AXIS_Y_] * b[STANDARD_AXIS_Z_] -
                              a[STANDARD_AXIS_Z_] * b[STANDARD_AXIS_Y_]);
  result[STANDARD_AXIS_Y_] = (a[STANDARD_AXIS_Z_] * b[STANDARD_AXIS_X_] -
                              a[STANDARD_AXIS_X_] * b[STANDARD_AXIS_Z_]);
  result[STANDARD_AXIS_Z_] = (a[STANDARD_AXIS_X_] * b[STANDARD_AXIS_Y_] -
                              a[STANDARD_AXIS_Y_] * b[STANDARD_AXIS_X_]);
  return result;
}

Eigen::Matrix3i
get_S3_permutation_matrix(int permutation)
{
  Eigen::Matrix3i S3_permutation_matrix;
  switch (permutation)
  {
    case 1:
      S3_permutation_matrix << 1, 0, 0, 0, 1, 0, 0, 0, 1;
      break;
    case 2:
      S3_permutation_matrix << 0, 1, 0, 0, 0, 1, 1, 0, 0;
      break;
    case 3:
      S3_permutation_matrix << 0, 0, 1, 1, 0, 0, 0, 1, 0;
      break;
    case 4:
      S3_permutation_matrix << 1, 0, 0, 0, 0, 1, 0, 1, 0;
      break;
    case 5:
      S3_permutation_matrix << 0, 0, 1, 0, 1, 0, 1, 0, 0;
      break;
    case 6:
      S3_permutation_matrix << 0, 1, 0, 1, 0, 0, 0, 0, 1;
      break;
    default:
      S3_permutation_matrix << 1, 0, 0, 0, 1, 0, 0, 0, 1;
      std::cerr
          << "ERROR:\tRequested invalid S3 permutation matrix (permutation: "
          << permutation << ")\n";
      break;
  }
  return S3_permutation_matrix;
}

int
get_S3_permuted_index(int permutation, int base_index)
{
  Eigen::Vector3i base_indizes(0, 1, 2);
  Eigen::Matrix3i permutation_matrix = get_S3_permutation_matrix(permutation);
  Eigen::Vector3i permuted_indizes   = permutation_matrix * base_indizes;
  return permuted_indizes(base_index);
}

int
get_S3_permutation_order(int permutation)
{
  switch (permutation)
  {
    case 1:
      return 1;
    case 2:
      return 3;
    case 3:
      return 3;
    case 4:
      return 2;
    case 5:
      return 2;
    case 6:
      return 2;
    default:
      return 0;
  }
}

double
get_derived_Leibnitz_element(Eigen::Matrix3d &determinante,
                             int permutation,
                             int axis_to_derive,
                             int coordinate_to_derive)
{
  /*
   * Expression for determinante:
   *
   * SUM_permutation [ sgn(permutation) * PRODUCT_i { Determinante(i, tau(i)) }
   * ]
   *
   * sgn(permutation) => parity of permutation (see: get_permutation_parity)
   */

  double parity = get_permutation_parity(get_S3_permutation_order(permutation));
  double element = (coordinate_to_derive == CENTRAL_ATOM_ ? (-parity) : parity);
  int axis, permuted_index;
  bool non_zero_element = false;

  for (axis = 0; axis < NUMBER_OF_AXIS_; ++axis)
  {
    if (coordinate_to_derive == CENTRAL_ATOM_)
    {
      /*
       * Central atom is present in every element.
       * To get the correct derivative use everything
       * but the element that is present in the
       * current derivative.
       */

      if (axis != axis_to_derive)
      {
        element *= determinante(axis, get_S3_permuted_index(permutation, axis));
        non_zero_element = true;
      }
    }
    /*
     * Dependending on the parity (permutation 1-3 vs. 4-6) the
     * selection of the elementes needed for the derivates
     * is different.
     *
     * This explantion will contain some permutation...
     *
     * All elements occure twice in the Leibnitz form.
     * One time in the positve parity part (permutation
     * 1-3) and one time in the negative parity
     * (permutation 4-6). To check which permutation
     * in the respective blocks contain the current
     * element one can make use of the permutation!
     * Each block will explain their choise individually.
     */
    else if (parity > .0)
    {
      /*
       * When differing axis "i" in the respective vectors
       * of the determinante (so a[i,1], a[i,2]& a[i,3])
       * the choise of the needed permutation (permutation
       * of Leibnitz summand!) follows the odd ordered
       * permutations!
       *
       *                [      Vector       ]
       * Axis 0(x) needs P(1), P(2) and P(3) for vectors 1, 2 and 3.
       * (Permutation 1) Axis 1(y) needs P(3), P(1) and P(2) for vectors 1, 2
       * and 3. (Permutation 3) Axis 2(z) needs P(2), P(3) and P(1) for vectors
       * 1, 2 and 3. (Permutation 2)
       *
       * So the needed permutation permutes respectively to the
       * current axis!
       * Due to the implementation (P(1) = 1) and the zero base
       * of the axis (X-axis = 0) we have to shift the values.
       *
       * The permutations you need (defined in the braces) permutes
       * according to P(4) = [ 1 3 2 ].
       * So if you want to check if your current permutation is
       * needed, you first need the permuted index of
       * (axis_to_derive+1) in permutation 4.
       * The result is an index (Range 0-2)! This index will now be
       * shifted according to the permutation (coordinate_to_derive+1).
       *
       * This will now result in an index, which encodes the needed
       * permutation in the positve parity (index + 1).
       *
       */
      if (permutation == (get_S3_permuted_index(
                              (coordinate_to_derive + 1),
                              get_S3_permuted_index(
                                  // Check the permutation for the current axis
                                  4, axis_to_derive)) +
                          1))
      {
        permuted_index = get_S3_permuted_index(permutation, axis);
        // Now extract the element we are deriving and build the product of the
        // others.
        if (!((coordinate_to_derive == permuted_index) &&
              (axis_to_derive == axis)))
        {
          non_zero_element = true;
          element *= determinante(axis, permuted_index);
        }
      }
    }
    else
    {
      /*
       * So, this is just like above, but since we need the negative
       * parity, we will look for the permutations 4-6. This means
       * we need a "(+3)-shift".
       *
       *                [         Vector          ]
       * Axis 0(x) needs P(4/1), P(6/3) and P(5/2) for vectors 1, 2 and 3.
       * (Permutation 4) Axis 1(y) needs P(6/3), P(5/2) and P(4/1) for vectors
       * 1, 2 and 3. (Permutation 5) Axis 2(z) needs P(5/2), P(4/1) and P(6/3)
       * for vectors 1, 2 and 3. (Permutation 6)
       *
       * Since the first step is just the permuation 1 we can just skip the
       * first encoding step.
       *
       * Just shift the axis_to_derive by 4 (axis 0 -> 4, axis 1 -> 5....)
       * and permute the coordinate_to_derive index by this value.
       * Now shift back by 4 to end up with the correct permutation and go!
       */
      if (permutation ==
          (get_S3_permuted_index((axis_to_derive + 4), coordinate_to_derive) +
           4))
      {
        permuted_index = get_S3_permuted_index(permutation, axis);
        if (!((coordinate_to_derive == permuted_index) &&
              (axis_to_derive == axis)))
        {
          non_zero_element = true;
          element *= determinante(axis, permuted_index);
        }
      }
    }
  }

  return (non_zero_element ? element : .0);
}

double
get_derived_Lebinitz_determinante(Eigen::Matrix3d &determinante,
                                  int axis_to_derive,
                                  int coordinate_to_derive)
{
  double derivative = .0;
  int permutation;
  for (permutation = 1; permutation <= NUMBER_OF_S3_PERMUTATIONS_;
       ++permutation)
  {
    derivative += get_derived_Leibnitz_element(
        determinante, permutation, axis_to_derive, coordinate_to_derive);
  }
  return derivative;
}

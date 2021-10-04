#include <Atom.hpp>
#include <DIDC.hpp>
#include <Declarations.hpp>
#include <LinAlg.hpp>
#include <SphericalHarmonics.hpp>
#include <SphericalHarmonicsMatrix.hpp>
#include <Structure.hpp>
#include <cblas.h>
#include <eigen3/Eigen/Cholesky>
#include <iomanip>
#include <iostream>
#include <phobos.h>

/*
 * Calculate the B-matrix according to equation 15 in Joel R. Tolman, Journal of
the American Chemical Society 2002, 124, 12020–12030, DOI 10.1021/ja0261123. The
producet Lambda*Lambda(tr) is computed using a Levenberg-Marquardt algorithm
followed by a Cholesky decomposition to obtain Lambda.
 */

int
calculate_lambda_B(Eigen::MatrixXd &U_D,
                   Eigen::MatrixXd &B_ref,
                   Eigen::MatrixXd &lambda)
{
  B_opt_data Bdata;
  int NOR, I;
  NOR         = U_D.rows();
  Bdata.UD    = (double *) malloc((NOR * COSINE_ELEMENTS_) * sizeof(double));
  Bdata.UDt   = (double *) malloc((NOR * COSINE_ELEMENTS_) * sizeof(double));
  Bdata.Ident = (double *) malloc((NOR * NOR) * sizeof(double));
  Bdata.work  = (double *) malloc((NOR * COSINE_ELEMENTS_) * sizeof(double));
  Bdata.sol   = (double *) malloc((NOR * NOR) * sizeof(double));
  Bdata.NOR   = NOR;
  Eigen::MatrixXd UD_red = Eigen::MatrixXd::Zero(NOR, 5);
  for (int i = 0; i < NOR; ++i)
  {
    for (int j = 0; j < COSINE_ELEMENTS_; ++j)
      UD_red(i, j) = U_D(i, j);
  }
  Eigen2Double(U_D, &(Bdata.UD), NOR, COSINE_ELEMENTS_);

  double *p, *x;
  p = (double *) malloc(25 * sizeof(double));
  x = (double *) malloc(NOR * sizeof(double));

  for (int i = I = 0; i < NOR; ++i)
  {
    for (int j = 0; j < NOR; ++j, ++I)
    {
      Bdata.Ident[I] = (i == j ? 1.0 : .0);
      if (I < 25)
        p[I] = .0;
    }
  }

  phobos(LambdaOpt, NULL, p, x, 25, 25, 1000, NULL, NULL, NULL, NULL,
         (void *) &Bdata);

  LambdaOpt(p, x, 25, 25, (void *) &Bdata);

  Double2Eigen(lambda, p, COSINE_ELEMENTS_, COSINE_ELEMENTS_);
  Eigen::MatrixXd tmp = UD_red * lambda * UD_red.transpose() -
                        Eigen::MatrixXd::Identity(NOR, NOR);

  Eigen::LLT<Eigen::MatrixXd> lltOfA(
      lambda); // compute the Cholesky decomposition of A
  Eigen::MatrixXd L =
      lltOfA.matrixL(); // retrieve factor L  in the decomposition

  double frob = .0;
  for (int i = 0; i < NOR; ++i)
    frob += (tmp(i, i) * tmp(i, i));

  B_ref = UD_red * L;

  return 0;
}

/*
 * Calculate the refinded B-matrix according to equation 17 in Joel R. Tolman,
Journal of the American Chemical Society 2002, 124, 12020–12030,
DOI 10.1021/ja0261123.
 */

int
refine_B_Matrix(Eigen::MatrixXd &B_md,
                Eigen::MatrixXd &D,
                Eigen::MatrixXd &B_ref)
{
  Eigen::MatrixXd MPI_B = MoorePenroseInverse(B_md);
  Eigen::MatrixXd MPI_D = MoorePenroseInverse(D);

  B_ref = D * MPI_D * B_md + B_md - B_md * MPI_B * D * MPI_D * B_md;

  return 0;
}

int
analyse_B_ref(Eigen::MatrixXd &B_ref, SphericalHarmonics *CurrY)
{
  Eigen::MatrixXd B;
  Eigen::MatrixXd tmpEVal(NUM_EIGENVALUES_, 1);
  Eigen::MatrixXd tmpEVec(EIGENVECTOR_ELEMENTS_, EIGENVECTOR_ELEMENTS_);
  Eigen::Vector3d ea;

  for (int i = 0; i < B_ref.rows(); ++i)
  {
    B  = Sv2St(B_ref.row(i).transpose());
    ea = Tensor2Euler(B, tmpEVal, tmpEVec, 0, true, false);

    CurrY->setOptimizedAngles(ea(1), ea(0));
    CurrY = CurrY->getNext();
  }
  return 0;
}

int
compute_B_Angles(Structure *CurrStruc)
{
  Eigen::MatrixXd B_ref;
  Eigen::MatrixXd lambda;
  Eigen::MatrixXd B_md =
      CurrStruc->getCosineMatrix(StructureOptions::Optimized);

  Eigen::MatrixXd DS = CurrStruc->getRDCmatrix(rdcMatrixOptions::Scaled);
  refine_B_Matrix(B_md, DS, B_ref);

  analyse_B_ref(B_ref, CurrStruc->getHeadYmatrix()->getHeadHarmonic());
  return 0;
}

void
LambdaOpt(double *p, double *x, int par, int mea, void *data)
{
  if (par != mea)
    std::cerr << "ERROR:\tStrange parameters in function LambdaOpt "
                 "(double*,double*,int,int,void*)...\n\t\tContact your trusted "
                 "nerd...\n";
  struct B_opt_data *dptr;
  dptr = (struct B_opt_data *) data;

  // Calculate estimate for hessian matrix
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dptr->NOR, 5, 5, 1.0,
              dptr->UD, 5, p, 5, .0, dptr->work, 5);


  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, dptr->NOR, dptr->NOR, 5,
              1.0, dptr->work, 5, dptr->UD, 5, .0, dptr->sol, dptr->NOR);


  for (int i = 0; i < (dptr->NOR); ++i)
    x[0] += ((dptr->sol[(i * (dptr->NOR + 1))] - 1.0) *
             (dptr->sol[(i * (dptr->NOR + 1))] - 1.0));
  x[0] = sqrt(x[0]);
  for (int i = 1; i < mea; ++i)
    x[i] = x[0];

  dptr = NULL;
  delete dptr;
}

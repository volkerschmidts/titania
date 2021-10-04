#include <cblas.h>
#include <lapacke.h>
#include <malloc.h>
#include <math.h>
#include <phobos.h>
#include <stdio.h>
#include <stdlib.h>


#define min(x, y)       (((x) < (y)) ? (x) : (y))
#define max(x, y)       (((x) > (y)) ? (x) : (y))
#define KRONDELTA(i, j) (((i) == (j)) ? (1.0) : (0.0))

void
copyPointer(double *initial, double *destination, int elements)
{
  int i;
  for (i = 0; i < elements; ++i)
    destination[i] = initial[i];
}

double
secondNorm(double *vec, int elements)
{
  int i;
  double mag = .0;
  for (i = 0; i < elements; ++i)
    mag += (vec[i] * vec[i]);
  return sqrt(mag);
}

double
infiniteNorm(double *vec, int elements)
{
  int i;
  double res = 0;
  for (i = 0; i < elements; ++i)
    res = max(res, fabs(vec[i]));
  return res;
}

double *
allocateWork(int worksize)
{
  int i;
  double *work = (double *) malloc((worksize) * sizeof(double));
  for (i = 0; i < (worksize); ++i)
    work[i] = .0;
  return work;
}

double
computeRho(
    void (*func)(double *parVec, double *meaVec, int par, int mea, void *data),
    double *parNew,
    double *delta,
    double *g,
    double *epsilon,
    double *epsilonNew,
    double *meaVec,
    double mu,
    int par,
    int mea,
    void *data)
{
  int i;
  double rho = .0, eMag, eNewMag, denom;
  func(parNew, epsilonNew, par, mea, data);
  for (i = 0; i < mea; ++i)
    epsilonNew[i] = (meaVec[i] - epsilonNew[i]);
  eMag    = secondNorm(epsilon, par);
  eNewMag = secondNorm(epsilonNew, par);
  denom   = .0;
  for (i = 0; i < par; ++i)
    denom += (delta[i] * (mu * delta[i] + g[i]));
  rho = (eMag * eMag - eNewMag * eNewMag) / denom;
  return rho;
}

void
numericalGradient(
    void (*func)(double *parVec, double *meaVec, int par, int mea, void *data),
    double *parVec,
    double *jac,
    double *work,
    int par,
    int mea,
    void *data)
{
  int i, j;
  double h;
  h = 1e-4;
  func(parVec, work, par, mea, data);
  double *tmp_mea = work + mea;
  for (i = 0; i < par; ++i)
  {
    parVec[i] += h;
    func(parVec, tmp_mea, par, mea, data);
    for (j = 0; j < mea; ++j)
      jac[j * mea + i] = ((tmp_mea[j] - work[j]) / h);
    parVec[i] -= h;
  }
  tmp_mea = NULL;
}

extern int
phobos(
    void (*func)(double *parVec, double *meaVec, int par, int mea, void *data),
    void (*jac)(double *parVec, double *jacobi, int par, int mea, void *data),
    double *parVec,
    double *meaVec,
    int par,
    int mea,
    int itmax,
    double *opts,
    double *info,
    double *work,
    double *covar,
    void *data)
{
  int i, j, index, freeWork, freeCovar, freeInfo, measXpar, parSq, meaSq,
      worksize, iter, incx, stop, numGrad;
  int *ipiv = (int *) malloc(par * sizeof(int));

  double alpha, beta, mu, muNew, nu, rho, tau, e1, e2, e3, maxHess,
      oneOverThree, compVal;
  double *jacobi, *hessian, *nMatrix, *currMeasure, *epsilon, *tmpMeasure, *g,
      *delta, *pNew, *numWork;

  numWork = NULL;
  if (jac == NULL)
  {
    numGrad = 1;
    numWork = (double *) malloc(2 * mea * sizeof(double));
  }
  else
    numGrad = 0;

  if (opts)
  {
    tau = opts[0];
    e1  = opts[1];
    e2  = opts[2];
    e3  = opts[3];
  }
  else
  {
    tau = STAN_TAU;
    e1  = STAN_EPSILON;
    e2  = STAN_EPSILON;
    e3  = STAN_EPSILON;
  }

  if (info)
    freeInfo = 0;
  else
  {
    info     = (double *) malloc(PHOBOS_INFOSIZE * sizeof(double));
    freeInfo = 1;
  }
  for (i = 0; i < PHOBOS_INFOSIZE; ++i)
    info[i] = .0;

  freeWork  = 0;
  freeCovar = 0;
  stop      = 0;
  iter      = 0;
  incx      = 1;

  measXpar = mea * par;
  parSq    = par * par;
  meaSq    = mea * mea;
  worksize = PHOBOS_WORKSIZE(par, mea);

  // cblas_dgemm implementation: C = alpha*A*B + beta*C
  alpha        = 1.0;
  beta         = .0;
  rho          = .0;
  nu           = 2.0;
  maxHess      = .0;
  oneOverThree = 1.0 / 3.0;
  compVal      = .0;

  // Handle the workspace
  if (work == NULL && covar == NULL)
  {
    work     = allocateWork((worksize + parSq));
    covar    = work + worksize;
    freeWork = 1;
  }
  else if (work == NULL || covar == NULL)
  {
    if (work == NULL)
    {
      work     = allocateWork(worksize);
      freeWork = 1;
    }
    else
    {
      covar     = allocateWork(parSq);
      freeCovar = 1;
    }
  }
  // Clean memory
  // Care with worksize: covar matrix is included
  // in standard implementation of PHOBOS_WORKSIZE.
  // This means we might run in trouble when just
  // running i < workspace ans user defined
  // workspace and covar seperatly.
  for (i = 0; i < (worksize - parSq); ++i)
    work[i] = .0;
  for (i = 0; i < parSq; ++i)
    covar[i] = .0;
  for (i = 0; i < mea; ++i)
    meaVec[i] = .0;

  // manage memory
  // for more information see phobos.h
  jacobi      = work;
  hessian     = (work + measXpar);
  nMatrix     = hessian + parSq;
  currMeasure = (nMatrix + parSq);
  epsilon     = currMeasure + mea;
  tmpMeasure  = epsilon + mea;
  g           = tmpMeasure + mea;
  delta       = g + par;
  pNew        = delta + par;

  // Calculate initial jacobi matrix
  if (numGrad)
    numericalGradient(func, parVec, jacobi, numWork, par, mea, data);
  else
    jac(parVec, jacobi, par, mea, data);

  // Calculate estimate for hessian matrix
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, par, par, mea, alpha,
              jacobi, par, jacobi, par, beta, hessian, par);

  // find maxHess
  for (i = 0; i < meaSq; i += (mea + 1))
    maxHess = max(maxHess, hessian[i]);

  // calculate mu = tau * maxHess
  mu = tau * maxHess;

  // Calculate current measurement vector
  func(parVec, currMeasure, par, mea, data);

  // Shift the arrays to fit the cblas implementation
  copyPointer(meaVec, epsilon, mea);
  copyPointer(currMeasure, tmpMeasure, mea);

  // calculate epsilon
  cblas_daxpy(mea, -alpha, tmpMeasure, incx, epsilon, incx);
  info[1] = secondNorm(epsilon, mea);

  // Calculate g = J(t) * epsilon
  cblas_dgemv(CblasRowMajor, CblasTrans, mea, par, alpha, jacobi, par, epsilon,
              incx, beta, g, incx);

  info[2] = infiniteNorm(g, par);
  if (info[2] <= e1)
  {
    info[7] = 5.0;
    info[8] = e1;
    stop    = 1;
  }
  while (!stop && (iter < itmax))
  {
    ++iter;
    info[0] += 1.0;
    do
    {
      //         printf("%d", ++it );
      for (i = 0; i < par; ++i)
      {
        for (j = 0; j < par; ++j)
        {
          index          = i * par + j;
          nMatrix[index] = hessian[index] + KRONDELTA(i, j) * mu;
        }
        delta[i] = g[i];
      }

      LAPACKE_dgesv(LAPACK_ROW_MAJOR, par, 1, nMatrix, par, ipiv, delta, 1);
      info[5] = secondNorm(delta, par);
      compVal = (e2 * secondNorm(parVec, par));
      if (info[5] <= compVal)
      {
        info[7] = 1.0;
        info[8] = compVal;
        stop    = 1;
      }
      else
      {
        for (j = 0; j < par; ++j)
          pNew[j] = parVec[j] + delta[j];
        rho = computeRho(func, pNew, delta, g, epsilon, tmpMeasure, meaVec, mu,
                         par, mea, data);
        //            printf( "  %f", rho);
        if (rho > .0)
        {
          copyPointer(pNew, parVec, par);

          if (numGrad)
            numericalGradient(func, parVec, jacobi, numWork, par, mea, data);
          else
            jac(parVec, jacobi, par, mea, data);

          cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, par, par, mea,
                      alpha, jacobi, par, jacobi, par, beta, hessian, par);

          // Calculate current measurement vector
          func(parVec, currMeasure, par, mea, data);
          copyPointer(meaVec, epsilon, par);
          copyPointer(currMeasure, tmpMeasure, par);

          // calculate epsilon
          cblas_daxpy(par, -alpha, tmpMeasure, incx, epsilon, incx);

          // Calculate g = J(t) * epsilon
          cblas_dgemv(CblasRowMajor, CblasTrans, mea, par, alpha, jacobi, par,
                      epsilon, incx, beta, g, incx);

          // calculate mu = tau * maxHess
          muNew = (1.0 - pow((2.0 * rho - 1), 3.0));
          mu *= (max(oneOverThree, muNew));
          nu      = 2.0;
          info[3] = secondNorm(epsilon, mea);
          muNew *= info[3];
          info[4] = infiniteNorm(g, par);
          if (((info[4]) <= e1))
          {
            stop    = 1;
            info[7] = 2.0;
            info[8] = e1;
          }
          else if ((info[3] * info[3] < e3))
          {
            stop    = 1;
            info[7] = 3.0;
            info[8] = e3;
          }
        }
        else
        {
          mu *= nu;
          nu *= 2.0;
        }
      }
    } while (!(rho < 0 || stop));
  }
  if (iter == itmax && info[7] == .0)
  {
    info[7] = 4.0;
    info[8] = (double) itmax;
  }
  for (i = 0; i < meaSq; i += (mea + 1))
    info[6] = max(info[6], hessian[i]);
  /*   info[1] = mu;
     info[2] = nu;
     info[3] = infiniteNorm(g, par);
     info[4] = rho;
  */
  // Clean the workspace
  free(ipiv);
  if (freeWork)
    free(work);
  if (freeCovar)
    free(covar);
  if (numGrad)
    free(numWork);
  if (freeInfo)
    free(info);
  return iter;
}

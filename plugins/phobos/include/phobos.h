
#ifndef INCLUDE_PHOBOS_H_
#define INCLUDE_PHOBOS_H_

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************
 *                                                     *
 *                     Memory map                      *
 *                                                     *
 *                        work                         *
 * jacobiMatrix       [ N(measure)   x N(parameter) ]  *
 * hessianMatrix      [ N(parameter) x N(parameter) ]  *
 * N-Matrix           [ N(parameter) x N(parameter) ]  *
 * covariance Matrix  [ N(parameter) x N(parameter) ]  *
 * x = f(p)           [ N(measure)   x      1       ]  *
 * eps. = x(i)-x(i-1) [ N(measure)   x      1       ]  *
 * x(tmp)             [ N(measure)   x      1       ]  *
 * g                  [ N(parameter) x      1       ]  *
 * delta              [ N(parameter) x      1       ]  *
 * p(new)             [ N(parameter) x      1       ]  *
 *                                                     *
 * Additional workspace:                               *
 * - Allocate the "dialog" p-/x-vectors for the user.  *
 *   (p[npar x 1], x[nmea x 1])                        *
 *                                                     *
 *                                                     *
 * By this the pointers are shifted by the respective  *
 * sizes of the matrizes.                              *
 * For "calculation" see phobos.c                      *
 *                                                     *
 *******************************************************/

#define PHOBOS_WORKSIZE(npar, nmeas) \
  (1 * (nmeas * npar) + 3 * (npar * npar) + 3 * nmeas + 3 * npar)

/******************************************************
 *                     Info map                       *
 *                                                    *
 * info[0]: number of iterations.                     *
 * info[1]: ||epsilon||(2) of p(initial).             *
 * info[2]: ||g||(inf) of p(initial).                 *
 * info[3]: ||epsilon||(2) of p(final).               *
 * info[4]: ||g||(inf) of p(final).                   *
 * info[5]: ||delta(p)||(2) at p(final).              *
 * info[6]: max(diag[hessian]) of p(final).           *
 * info[7]: reason for termination.                   *
 *          - 1: low delta(p)                         *
 *          - 2: low g(p)                             *
 *          - 3: low epsilon^2(p)                     *
 *          - 4: max ierations                        *
 *          - 5: low initial g(p)                     *
 * info[8]: respective criterion for info[7].         *
 *                                                    *
 *                                                    *
 *                                                    *
 ******************************************************/


#define PHOBOS_INFOSIZE 9

#define STAN_THRESH     1E-17
#define STAN_DIFF_DELTA 1E-06
#define STAN_EPSILON    1E-15
#define STAN_TAU        1E-03


int
phobos(void (*func)(double *parVec, double *meaVec, int m, int n, void *data),
       void (*jac)(double *parVec, double *jacobi, int m, int n, void *data),
       double *parVec,
       double *meaVec,
       int m,
       int n,
       int itmax,
       double *opts,
       double *info,
       double *work,
       double *covar,
       void *adata);

#ifdef __cplusplus
}
#endif

double computeRho(
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
    void *data);

void numericalGradient(
    void (*func)(double *parVec, double *meaVec, int par, int mea, void *data),
    double *p,
    double *jac,
    double *work,
    int par,
    int mea,
    void *data);

double *allocateWork(int worksize);

void copyPointer(double *initial, double *destination, int elements);

double secondNorm(double *vec, int elements);

double infiniteNorm(double *vec, int elements);

#endif

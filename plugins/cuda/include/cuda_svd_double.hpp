
#ifndef CUDA_SVD_DOUBLE_
#define CUDA_SVD_DOUBLE_

#include <cusolverDn.h>

int
cuda_double_svd(double *, double *, double *, double *, int, int, int, bool);
int cuda_double_svd_device_memory(double *,
                                  double *,
                                  double *,
                                  double *,
                                  int,
                                  int,
                                  int);


// int CUDA_DOUBLE_SVD_GPU_MEMORY( double*, double*, double*, double*, int, int,
// int, int, int* );

#endif // CUDA_SVD_DOUBLE_

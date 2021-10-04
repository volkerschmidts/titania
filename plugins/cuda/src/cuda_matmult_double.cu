
#include "../../cuda/include/cuda_matmult_double.hpp"

#include <cuda_runtime.h>
#include <cublas.h>

void cuda_dgemm(
             char transa,
             char transb,
             const int &m,
             const int &n,
             const int &k,
             const double alpha,
             const double *A,
             const int &lda,
             const double *B,
             const int &ldb,
             const double beta,
             double *C,
             const int &ldc)
{
   cublasDgemm( transa,
                transb,
                m,
                n,
                k,
                alpha,
                A,
                lda,
                B,
                ldb,
                beta,
                C,
                ldc);
}

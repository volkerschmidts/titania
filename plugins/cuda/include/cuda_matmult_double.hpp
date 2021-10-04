
#ifndef CUDA_MATMULT_DOUBLE_
#define CUDA_MATMULT_DOUBLE_


void cuda_dgemm(char transa,
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
                const int &ldc);

#endif

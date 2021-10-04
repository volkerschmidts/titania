
#include "../../cuda/include/cuda_svd_double.hpp"

#include <iostream>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include "../include/cuda_svd.hpp"
#include "../include/cuda_helper.hpp"


int cuda_double_svd( double *A, double *S, double *U, double *V, int m, int n, int lda, bool reset )
{
   cudaError_t cudaStat1 = cudaSuccess;
   cudaError_t cudaStat2 = cudaSuccess;
   cudaError_t cudaStat3 = cudaSuccess;
   cudaError_t cudaStat4 = cudaSuccess;
   cudaError_t cudaStat5 = cudaSuccess;

   // Device pointers
   double *d_A, *d_S, *d_U, *d_V;
   constexpr size_t size_d = sizeof(double);
   d_A = d_S = d_U = d_V = NULL;

   const long unsigned int S_size   = get_size_(  1, min(m,n), size_d);
   const long unsigned int U_size   = get_size_(  m, min(m,n), size_d);
   const long unsigned int V_size   = get_size_(  n, min(m,n), size_d);
   const long unsigned int A_size   = get_size_(lda,       n , size_d);

   // Allocate memory on GPU
   cudaStat1 = cudaMalloc ((void**)&d_A   , A_size);
   cudaStat2 = cudaMalloc ((void**)&d_S   , S_size);
   cudaStat3 = cudaMalloc ((void**)&d_U   , U_size);
   cudaStat4 = cudaMalloc ((void**)&d_V   , V_size);

   if (cudaSuccess != cudaStat1) return ERROR_CUDA_MALLOC;
   if (cudaSuccess != cudaStat2) return ERROR_CUDA_MALLOC;
   if (cudaSuccess != cudaStat3) return ERROR_CUDA_MALLOC;
   if (cudaSuccess != cudaStat4) return ERROR_CUDA_MALLOC;

   // Copy initial matrix on GPU
   cudaStat1 = cudaMemcpy(d_A, A, A_size, cudaMemcpyHostToDevice);
   
   if (cudaSuccess != cudaStat1)  return ERROR_CUDA_COPY;

   cuda_double_svd_device_memory( d_A, d_S, d_U, d_V, m, n, lda );

   // Copy results
   cudaStat1 = cudaMemcpy(U, d_U, U_size, cudaMemcpyDeviceToHost);
   cudaStat2 = cudaMemcpy(V, d_V, V_size, cudaMemcpyDeviceToHost);
   cudaStat3 = cudaMemcpy(S, d_S, S_size, cudaMemcpyDeviceToHost);
   cudaStat5 = cudaDeviceSynchronize();

   if (cudaSuccess != cudaStat1) return ERROR_CUDA_COPY;
   if (cudaSuccess != cudaStat2) return ERROR_CUDA_COPY;
   if (cudaSuccess != cudaStat3) return ERROR_CUDA_COPY;
   if (cudaSuccess != cudaStat4) return ERROR_CUDA_COPY;
   if (cudaSuccess != cudaStat5) return ERROR_CUDA_COPY;

// Free device memory
   if (d_A    ) cudaFree(d_A);
   if (d_S    ) cudaFree(d_S);
   if (d_U    ) cudaFree(d_U);
   if (d_V    ) cudaFree(d_V);

   if ( reset ) cudaDeviceReset();
   return 0;
}

int cuda_double_svd_device_memory( double *d_A, double *d_S, double *d_U, double *d_V, int m, int n, int lda )
{
   cusolverDnHandle_t cusolverH = NULL;
   cudaStream_t stream = NULL;

   cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
   cudaError_t cudaStat = cudaSuccess;

   // Device pointers
   double *d_work;
   d_work = NULL;

   int info, return_value;
   int *d_info = NULL;  /* error info */
   int lwork = 0;       /* size of workspace */

   return_value = 0;

   // Create cusolver handle
   status = cusolverDnCreate(&cusolverH);
   if (CUSOLVER_STATUS_SUCCESS != status ) return ERROR_CUSOLVER_DN_CREATE;

   cudaStat = cudaMalloc ((void**)&d_info, sizeof(int));
   if (cudaSuccess != cudaStat) { return_value = ERROR_CUDA_MALLOC; goto cleanup; }

   // Querry buffer
   status = cusolverDnDgesvd_bufferSize( cusolverH, m, n, &lwork );
   if ( CUSOLVER_STATUS_SUCCESS != status ) { return_value = ERROR_CUDA_BUFFER_SIZE; goto cleanup; }

   // Allocate workspace
   cudaStat = cudaMalloc((void**)&d_work , sizeof(double)*lwork);
   if (cudaSuccess != cudaStat ) { return_value = ERROR_CUDA_MALLOC; goto cleanup; }

   // Compute SVD
   status = cusolverDnDgesvd (
                              cusolverH,
                              'S',        // compute vectors in the range of A
                              'S',        // compute vectors in the range of A
                              m,
                              n,
                              d_A,
                              lda,
                              d_S,
                              d_U,
                              m,
                              d_V,
                              n,          // The lead dimension is smaller due to the fact that we calculate reduced vectors
                              d_work,
                              lwork,
                              NULL,
                              d_info
                             );

   // wait for the calculation
   cudaStat = cudaDeviceSynchronize();
   if ( cudaSuccess != cudaStat ) { std::cout << "cudastat = " << cudaStat << std::endl; return_value = ERROR_CUDA_DEVICE_SYNCHRONIZE; goto cleanup; }
   if ( CUSOLVER_STATUS_SUCCESS != status )
   {
      cudaStat = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
      if (cudaSuccess == cudaStat)
      {
         if ( info ) std::cerr << "ERROR:\tcudasolverDnDgesvd error code: " << info << std::endl;
      }
      else
      {
         std::cerr << "ERROR:\tcudasolverDnDgesvd failed. TITANIA was not able to get addional information on the\n"
                   << "\tproblem via cudasolverDnDgesvd info...\n";
      }
      return_value = ERROR_CUDA_SVD_S;
   }

   cleanup:
      if (d_work ) cudaFree(d_work);
      if (cusolverH) cusolverDnDestroy(cusolverH);
      if (stream   ) cudaStreamDestroy(stream);
      return return_value;
}

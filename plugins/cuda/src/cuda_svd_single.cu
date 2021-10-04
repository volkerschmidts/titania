
#include "../../cuda/include/cuda_svd_single.hpp"

#include <iostream>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include "../../cuda/include/cuda_svd.hpp"

/*constexpr float mb_fac = 1024.0*1024.0;

void check_memory()
{
   size_t free_byte ;
   size_t total_byte ;
   cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
   if ( cudaSuccess != cuda_status )
   {
      printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
      exit(1);
   }

   float free_db = (float)free_byte ;
   float total_db = (float)total_byte ;
   float used_db = total_db - free_db ;

   printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n", 
           used_db/mb_fac, 
           free_db/mb_fac, 
           total_db/mb_fac);
}*/

int cuda_single_svd( float *A, float *S, float *U, float *V, int m, int n, int lda )
{
   cusolverDnHandle_t cusolverH = NULL;
   cudaStream_t stream = NULL;
   gesvdjInfo_t gesvdj_params = NULL;

   cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
   cudaError_t cudaStat1 = cudaSuccess;
   cudaError_t cudaStat2 = cudaSuccess;
   cudaError_t cudaStat3 = cudaSuccess;
   cudaError_t cudaStat4 = cudaSuccess;
   cudaError_t cudaStat5 = cudaSuccess;

   // Device pointers
   float *d_work, *d_A, *d_S, *d_U, *d_V;
   d_work = d_A = d_S = d_U = d_V = NULL;

   int *d_info = NULL;  /* error info */
   int lwork = 0;       /* size of workspace */

   // Configuration of gesvdj
   const float tol = 1.e-7;
   const int econ = 0; /* econ = 1 for economy size */

   // Numerical results of gesvdj
//   double residual = 0;
//   int executed_sweeps = 0;

   // Create cusolver handle
   status = cusolverDnCreate(&cusolverH);
   if (CUSOLVER_STATUS_SUCCESS != status ) return ERROR_CUSOLVER_DN_CREATE;

   cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
   if ( cudaSuccess != cudaStat1 ) return ERROR_STREAM_CREATE_WITH_FLAGS;
   
   // Bind a stream
   status = cusolverDnSetStream(cusolverH, stream);
   if ( CUSOLVER_STATUS_SUCCESS != status ) return ERROR_SET_STREAM;

   // Configure gesvdj 
   status = cusolverDnCreateGesvdjInfo(&gesvdj_params);
   if ( CUSOLVER_STATUS_SUCCESS != status ) return ERROR_CREATE_SVD_INFO;

   // Set tolerance
   status = cusolverDnXgesvdjSetTolerance( gesvdj_params, tol);
   if ( CUSOLVER_STATUS_SUCCESS != status ) return ERROR_SET_SVD_PARAMETER;

   // Set sweeps
//   status = cusolverDnXgesvdjSetMaxSweeps( gesvdj_params, max_sweeps);
//   if ( CUSOLVER_STATUS_SUCCESS != status ) return ERROR_SET_SVD_PARAMETER;

   // Allocate memory on GPU
   cudaStat1 = cudaMalloc ((void**)&d_A   , sizeof(float)*lda*n);
   cudaStat2 = cudaMalloc ((void**)&d_S   , sizeof(float)*n);
   cudaStat3 = cudaMalloc ((void**)&d_U   , sizeof(float)*lda*m);
   cudaStat4 = cudaMalloc ((void**)&d_V   , sizeof(float)*lda*n);
   cudaStat5 = cudaMalloc ((void**)&d_info, sizeof(int));
   
   if (cudaSuccess != cudaStat1) return ERROR_CUDA_MALLOC;
   if (cudaSuccess != cudaStat2) return ERROR_CUDA_MALLOC;
   if (cudaSuccess != cudaStat3) return ERROR_CUDA_MALLOC;
   if (cudaSuccess != cudaStat4) return ERROR_CUDA_MALLOC;
   if (cudaSuccess != cudaStat5) return ERROR_CUDA_MALLOC;

   // Copy initial matrix on GPU
   cudaStat1 = cudaMemcpy(d_A, A, sizeof(float)*lda*n, cudaMemcpyHostToDevice);
   
   if (cudaSuccess != cudaStat1)  return ERROR_CUDA_COPY;
 
   // Querry buffer
   status = cusolverDnSgesvdj_bufferSize(
                                         cusolverH,
                                         CUSOLVER_EIG_MODE_VECTOR,
                                         econ,
                                         m,
                                         n,
                                         d_A,
                                         lda,
                                         d_S,
                                         d_U,
                                         lda,
                                         d_V,
                                         lda,
                                         &lwork,
                                         gesvdj_params);
   if ( CUSOLVER_STATUS_SUCCESS != status ) return ERROR_CUDA_BUFFER_SIZE;

   // Allocate workspace
   cudaStat1 = cudaMalloc((void**)&d_work , sizeof(float)*lwork);
   if (cudaSuccess != cudaStat1 ) return ERROR_CUDA_MALLOC;

   // Compute SVD
   status = cusolverDnSgesvdj(
                              cusolverH,
                              CUSOLVER_EIG_MODE_VECTOR,
                              econ,
                              m,
                              n,
                              d_A,
                              lda,
                              d_S,
                              d_U,
                              lda,
                              d_V,
                              lda,
                              d_work,
                              lwork,
                              d_info,
                              gesvdj_params);
   // wait for the calculation
   cudaStat1 = cudaDeviceSynchronize();
   if ( CUSOLVER_STATUS_SUCCESS != status ) return ERROR_CUDA_SVD_S;
   if ( cudaSuccess != cudaStat1 ) return ERROR_CUDA_DEVICE_SYNCHRONIZE;

   // Copy results
   cudaStat1 = cudaMemcpy(U, d_U, sizeof(float)*lda*m, cudaMemcpyDeviceToHost);
   cudaStat2 = cudaMemcpy(V, d_V, sizeof(float)*lda*n, cudaMemcpyDeviceToHost);
   cudaStat3 = cudaMemcpy(S, d_S, sizeof(float)*n    , cudaMemcpyDeviceToHost);
   cudaStat5 = cudaDeviceSynchronize();

   if (cudaSuccess != cudaStat1) return ERROR_CUDA_COPY;
   if (cudaSuccess != cudaStat2) return ERROR_CUDA_COPY;
   if (cudaSuccess != cudaStat3) return ERROR_CUDA_COPY;
   if (cudaSuccess != cudaStat4) return ERROR_CUDA_COPY;
   if (cudaSuccess != cudaStat5) return ERROR_CUDA_COPY;

//   status = cusolverDnXgesvdjGetSweeps( cusolverH, gesvdj_params, &executed_sweeps);
//   if ( CUSOLVER_STATUS_SUCCESS != status ) return ERROR_READ_INFORMATION;

//   status = cusolverDnXgesvdjGetResidual( cusolverH, gesvdj_params, &residual);
//   if (CUSOLVER_STATUS_SUCCESS != status) return ERROR_READ_INFORMATION;

// Free device memory
   if (d_A    ) cudaFree(d_A);
   if (d_S    ) cudaFree(d_S);
   if (d_U    ) cudaFree(d_U);
   if (d_V    ) cudaFree(d_V);
   if (d_work ) cudaFree(d_work);

   if (cusolverH) cusolverDnDestroy(cusolverH);
   if (stream      ) cudaStreamDestroy(stream);
   if (gesvdj_params) cusolverDnDestroyGesvdjInfo(gesvdj_params);

   cudaDeviceReset();
   return 0;
}


#include "../include/cuda_mpi_double.hpp"

#include <cuda_runtime.h>
#include <cusolverDn.h>

#include <iostream>

#include "../include/cuda_matmult_double.hpp"
#include "../include/cuda_svd.hpp"
#include "../include/cuda_svd_double.hpp"
#include "../include/cuda_get_information.hpp"
#include "../include/cuda_helper.hpp"

#define THREADS_PER_BLOCK 512

   constexpr size_t size_d = sizeof(double);

__global__ void cuda_invert_singular_values_device_memory( double *d_S, double *d_S_inv, const double cutoff, const int num_of_sigma, const int num_of_el )
{
   int index = threadIdx.x + blockIdx.x * blockDim.x;
   if ( index < num_of_el )
   {
      int sigma_index = ( index / num_of_sigma );
      if ( sigma_index >= num_of_sigma ) d_S_inv[index] = .0;
      else
      {
         bool is_sigma = (( index - sigma_index ) % num_of_sigma ) == 0;
         if ( is_sigma && d_S[sigma_index] > cutoff ) d_S_inv[index] = ( 1.0 / d_S[sigma_index] );
         else d_S_inv[index] = .0;
      }
   }
}

void printDeviceMatrix ( double* M, int m, int n, int lda )
{
   double *tmp = (double*) malloc( m*n*sizeof(double));
   cudaMemcpy(tmp, M, lda*n*sizeof(double),cudaMemcpyDeviceToHost);
   for ( int r = 0; r < m; ++r )
   {
      for ( int c = 0; c < n; ++c ) std::cout << tmp[r+c*lda] << "  ";
      std::cout << std::endl;
   }
   free(tmp);
}

int check_needed_memory ( int lda, int m, int n, size_t type_size )
{
   // Calculate the needed space on dRAM if all matrizes are kept on it,
   // but are deleted as soon as they are not needed anymore.
   long unsigned int requested_matrizes = type_size * n * ( m + 3*n );

   // Calculate the needed space on dRAM if S_inv and U share the memory
   // since they are the largest matrizes. U will be buffered on main RAM.
   long unsigned int reducable_size = type_size * n * ( m + n +1 );

   // Additional worksize of cusolverDNDgesvd
   cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
   cusolverDnHandle_t cusolverH = NULL;

   int work = 0;
   long unsigned int free_i;
   size_t free_s, total;

   // Create cusolver handle
   status = cusolverDnCreate(&cusolverH);
   if (CUSOLVER_STATUS_SUCCESS != status ) return ERROR_CUSOLVER_DN_CREATE;

   // Querry buffer
   status = cusolverDnDgesvd_bufferSize( cusolverH, m, n, &work );
   if ( CUSOLVER_STATUS_SUCCESS != status ) return ERROR_CUDA_BUFFER_SIZE;

   // Run cudaMemGetInfor here since cuslverDnCreate allocates
   // space for the cusolverH
   cudaMemGetInfo(&free_s,&total);

   if (cusolverH) cusolverDnDestroy(cusolverH);

   work *= type_size;
   requested_matrizes += ((long unsigned int) work);
   reducable_size += ((long unsigned int) work);

   free_i = (long unsigned int) free_s;

   if ( free_i > requested_matrizes && (free_i - requested_matrizes) > CUDA_SAVE_MEMORY_BUFFER_ ) return 0;
   else if ( free_i > reducable_size && (free_i - reducable_size) > CUDA_SAVE_MEMORY_BUFFER_ ) return 1;

   return -1;
}

int cuda_double_moore_penrose_inverse_host_memory( double *A, double **A_inv, int m, int n, int lda, double cutoff )
{
   cudaError_t cudaStat = cudaSuccess;

   int return_value = 0;
   const long unsigned int A_size =   get_size_(lda,       n , size_d);
   const long unsigned int A_i_size = get_size_(  m,       n , size_d);
   double *d_A;      // Initial matrix A
   double *d_A_inv;      // Initial matrix A

   cudaStat = cudaMalloc ((void**)&d_A, A_size); if (cudaSuccess != cudaStat) { return_value = ERROR_CUDA_MALLOC; goto cleanup; }
   cudaStat = cudaMalloc ((void**)&d_A_inv, A_i_size); if (cudaSuccess != cudaStat) { return_value = ERROR_CUDA_MALLOC; goto cleanup; }

   return_value = cuda_double_moore_penrose_inverse_device_memory( d_A, d_A_inv, m, n, lda, cutoff );

   A_inv[0] = (double*) malloc(A_size);
   cudaStat = cudaMemcpy(A_inv[0], d_A_inv, A_size, cudaMemcpyDeviceToHost);

   cleanup:
      if (d_A)     cudaFree(d_A);
      if (d_A_inv) cudaFree(d_A_inv);
      cudaDeviceReset();
      return return_value;
}

int cuda_double_moore_penrose_inverse_device_memory( double *d_A, double *d_A_inv, int m, int n, int lda, double cutoff )
{
   cudaError_t cudaStat = cudaSuccess;

   // Start with checking the expected size of the problem.
   int memory_information = check_needed_memory (lda, m, n, sizeof(double));
   if ( memory_information == -1 )
   {
      std::cerr << "ERROR:\tNot enough GPU memory...\n";
      return ERROR_CUDA_MALLOC;
   }

   // Some values for the following calculations
   int return_value = 0;
   const int S_elements = min(m,n)*n;

   const long unsigned int S_size   = get_size_(  1, min(m,n), size_d);
   const long unsigned int S_i_size = get_size_(  n, min(m,n), size_d);
   const long unsigned int U_size   = get_size_(  m, min(m,n), size_d);
   const long unsigned int V_size   = get_size_(  n, min(m,n), size_d);
   const long unsigned int A_size   = get_size_(lda,       n , size_d);

   // Device pointers
//   double *d_A;      // Initial matrix A
   double *d_S;      // Vector of singular values
   double *d_S_inv;  // diag[ 1 / sigma_i ]
   double *d_U;      // SVD matrix U
   double *d_V;      // SVD matrix V
   double *d_VS;     // V * diag[1/sigma_i]
   double *U;
   d_S = d_S_inv = d_U = d_V = d_VS = U = NULL;

   // Allocate svd memory on GPU
   cudaStat = cudaMalloc ((void**)&d_U, U_size); if (cudaSuccess != cudaStat) { return_value = ERROR_CUDA_MALLOC; goto cleanup; }
//   cudaStat = cudaMalloc ((void**)&d_A, A_size); if (cudaSuccess != cudaStat) { return_value = ERROR_CUDA_MALLOC; goto cleanup; }
   cudaStat = cudaMalloc ((void**)&d_V, V_size); if (cudaSuccess != cudaStat) { return_value = ERROR_CUDA_MALLOC; goto cleanup; }
   cudaStat = cudaMalloc ((void**)&d_S, S_size); if (cudaSuccess != cudaStat) { return_value = ERROR_CUDA_MALLOC; goto cleanup; }

   if ( !memory_information )
   {
      cudaStat = cudaMalloc ((void**)&d_S_inv, S_i_size);
      if (cudaSuccess != cudaStat)
      {
         return_value = ERROR_CUDA_MALLOC;
         goto cleanup;
      }
   }

   // Copy initial matrix on GPU
//   cudaStat = cudaMemcpy(d_A, A, A_size, cudaMemcpyHostToDevice); if (cudaSuccess != cudaStat) { return_value = ERROR_CUDA_COPY; goto cleanup; }

   // Compute the svd
   return_value = cuda_double_svd_device_memory( d_A, d_S, d_U, d_V, m, n, lda );
   if ( return_value )
   {
      std::cout << "svd_error = " << return_value << std::endl;
      goto cleanup;
   }

   if ( memory_information )
   {
      U = (double*) malloc (U_size);
      cudaStat = cudaMemcpy(U, d_U, U_size, cudaMemcpyDeviceToHost); if (cudaSuccess != cudaStat) { return_value = ERROR_CUDA_COPY; goto cleanup; }
//      if (d_A) { cudaFree(d_A); d_A = NULL; }
      if (d_U) { cudaFree(d_U); d_U = NULL; }
      cudaStat = cudaMalloc ((void**)&d_S_inv, S_i_size); if (cudaSuccess != cudaStat) { return_value = ERROR_CUDA_MALLOC; goto cleanup; }
   }

   // Invert the singular values and write them on the diagonal of S_inv
   if ( n < THREADS_PER_BLOCK )
      cuda_invert_singular_values_device_memory<<<m,n>>> ( d_S, d_S_inv, cutoff, n, S_elements );
   else
   {
      int Number_of_blocks = (S_elements/THREADS_PER_BLOCK+1);  // +1 -> else last sigma might be skiped!
      cuda_invert_singular_values_device_memory<<<Number_of_blocks,THREADS_PER_BLOCK>>> ( d_S, d_S_inv, cutoff, n, S_elements );
   }

   // Free the vector S
   if (d_S) { cudaFree(d_S); d_S = NULL; }

   // Allocate the memory for the product V*S_inv
   cudaStat = cudaMalloc ((void**)&d_VS, V_size); if (cudaSuccess != cudaStat) { return_value = ERROR_CUDA_MALLOC; goto cleanup; }

   // Compute the product
   // If "small svd is used we have to choose 't' for matrix d_V!
   cuda_dgemm('t', 'n', n, n, n, 1.0, d_V,  n, d_S_inv,   n, 0.0,    d_VS, n );

   // Free V since it is not needed anymore
   if (d_V) { cudaFree(d_V); d_V = NULL; }
   if (d_S_inv) { cudaFree(d_S_inv); d_S_inv = NULL; }

   // If U was shifted to main RAM write the respective values back
   // in initial dRAM and remove S_inv
   if ( memory_information)
   {
      cudaStat = cudaMalloc ((void**)&d_U, U_size); if (cudaSuccess != cudaStat) { return_value = ERROR_CUDA_MALLOC; goto cleanup; }
//      cudaStat = cudaMalloc ((void**)&d_A, U_size); if (cudaSuccess != cudaStat) return ERROR_CUDA_MALLOC;
      cudaStat = cudaMemcpy(d_U, U, U_size, cudaMemcpyHostToDevice); if (cudaSuccess != cudaStat) { return_value = ERROR_CUDA_COPY; goto cleanup; }
      free(U); U = NULL;
   }

   // Allocate the memory for the final matrix A
   //cudaStat = cudaMalloc ((void**)&d_A_inv, A_size); if (cudaSuccess != cudaStat) { return_value = ERROR_CUDA_MALLOC; goto cleanup; }

   // And perform the final product
   cuda_dgemm('n', 't', n, m, n, 1.0, d_VS,   n,     d_U, lda, 0.0, d_A_inv, n );

   // Copy A_inverse back on host memory
//   A_inv[0] = (double*) malloc( A_size);
//   cudaStat = cudaMemcpy(A_inv[0], d_A, A_size, cudaMemcpyDeviceToHost);

// Free device memory and make sure that everyhting is freed!
   cleanup:
//      if (d_A)     cudaFree(d_A);
      if (d_S)     cudaFree(d_S);
      if (d_S_inv) cudaFree(d_S_inv);
      if (d_U)     cudaFree(d_U);
      if (d_V)     cudaFree(d_V);
      if (d_VS)    cudaFree(d_VS);
      if (U)       free(U);

//   cudaDeviceReset();
   return return_value;
}

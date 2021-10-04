
#ifndef CUDA_MPI_DOUBLE_
#define CUDA_MPI_DOUBLE_

#ifndef ZERO_CUTOFF_
#define ZERO_CUTOFF_ 1e-9
#endif

int cuda_double_moore_penrose_inverse_host_memory(double *,
                                                  double **,
                                                  int,
                                                  int,
                                                  int,
                                                  double);

int cuda_double_moore_penrose_inverse_device_memory(double *,
                                                    double *,
                                                    int,
                                                    int,
                                                    int,
                                                    double);

#endif

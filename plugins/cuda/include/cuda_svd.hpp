
#ifndef CUDA_SVD_
#define CUDA_SVD_

#define ERROR_CUSOLVER_DN_CREATE       1
#define ERROR_STREAM_CREATE_WITH_FLAGS 2
#define ERROR_SET_STREAM               3
#define ERROR_CREATE_SVD_INFO          4
#define ERROR_SET_SVD_PARAMETER        5
#define ERROR_CUDA_MALLOC              6
#define ERROR_CUDA_COPY                7
#define ERROR_CUDA_BUFFER_SIZE         8
#define ERROR_CUDA_DEVICE_SYNCHRONIZE  9
#define ERROR_CUDA_SVD_S               10
#define ERROR_READ_INFORMATION         11

constexpr auto CUDA_SAVE_MEMORY_BUFFER_ = 1024 * 1024 * 100;

#endif // CUDA_SVD

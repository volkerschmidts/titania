
#include <cuda_helper.hpp>
#include <iostream>
std::string get_gpu_device_name ( int count )
{
   cudaDeviceProp deviceProp;
   cudaGetDeviceProperties (&deviceProp, count);
   cudaDeviceReset();
   return deviceProp.name;
}

double get_gpu_memory ( int count )
{
   cudaDeviceProp deviceProp;
   cudaGetDeviceProperties (&deviceProp, count);
   cudaDeviceReset();
   return (deviceProp.totalGlobalMem / (1024.0*1024.0));
}

int get_number_of_gpu_devices()
{
   int count;
   cudaError_t error_code = cudaGetDeviceCount(&count);
   if ( error_code ) return 0;
   cudaDeviceReset();
   return count;
}

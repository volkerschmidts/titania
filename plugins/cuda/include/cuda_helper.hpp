
#ifndef CUDA_HELPER_HPP_
#define CUDA_HELPER_HPP_

#include <string>

inline long unsigned int
get_size_(int m, int n, size_t type_size)
{
  long unsigned int r = ((long unsigned int) n) * ((long unsigned int) m);
  return (r * type_size);
}

std::string get_gpu_device_name(int count = 0);

double get_gpu_memory(int count = 0);

int get_number_of_gpu_devices();

#endif /* CUDA_HELPER_HPP_ */

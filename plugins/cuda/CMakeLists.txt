project(titania_cuda)
enable_language(CUDA)

include(GNUInstallDirs)

add_library(
  ${PROJECT_NAME}
  src/cuda_get_information.cu
  src/cuda_helper.cu
  src/cuda_matmult_double.cu
  src/cuda_mpi_double.cu
  src/cuda_svd_double.cu
  src/cuda_svd_single.cu)

target_include_directories(${PROJECT_NAME} PRIVATE include)

set_target_properties(
  ${PROJECT_NAME}
  PROPERTIES CUDA_SEPARABLE_COMPILATION ON # something changed in recent CMake versions
                                           # ...
             CUDA_ARCHITECTURES OFF)

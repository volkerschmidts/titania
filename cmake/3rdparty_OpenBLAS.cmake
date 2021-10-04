include(ExternalProject)
ExternalProject_Add(
  project_lapacke
  GIT_REPOSITORY "https://github.com/xianyi/OpenBLAS.git"
  GIT_TAG "v0.3.17"
  GIT_SHALLOW ON
  PREFIX ${3RDPARTY_DIR}
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DBUILD_SHARED_LIBS=off
             -DUSE_THREAD=ON)
add_library(LAPACKE::LAPACKE STATIC IMPORTED)
set_target_properties(
  LAPACKE::LAPACKE
  PROPERTIES INTERFACE_COMPILE_OPTIONS "-pthread"
             INTERFACE_INCLUDE_DIRECTORIES "${3RDPARTY_DIR}/include/openblas"
             INTERFACE_LINK_LIBRARIES "gfortran;pthread"
             IMPORTED_LOCATION "${3RDPARTY_DIR}/lib/libopenblas.a")

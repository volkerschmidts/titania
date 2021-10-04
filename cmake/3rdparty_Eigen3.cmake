set(3RDPARTY_DIR ${CMAKE_BINARY_DIR}/3rdparty)

include(ExternalProject)
ExternalProject_Add(
  project_eigen3
  GIT_REPOSITORY "https://gitlab.com/libeigen/eigen.git"
  GIT_TAG "3.3.9"
  GIT_SHALLOW ON
  PREFIX ${3RDPARTY_DIR}
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>)
add_library(Eigen3::Eigen INTERFACE IMPORTED)
set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                               "${3RDPARTY_DIR}/include/eigen3")

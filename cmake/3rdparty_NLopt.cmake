set(3RDPARTY_DIR ${CMAKE_BINARY_DIR}/3rdparty)

include(ExternalProject)
ExternalProject_Add(
  project_nlopt
  GIT_REPOSITORY "https://github.com/stevengj/nlopt.git"
  GIT_TAG "v2.7.0"
  GIT_SHALLOW ON
  PREFIX ${3RDPARTY_DIR}
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
             -DBUILD_SHARED_LIBS=off
             -DNLOPT_CXX=on
             -DNLOPT_PYTHON=off
             -DNLOPT_OCTAVE=off
             -DNLOPT_MATLAB=off
             -DNLOPT_GUILE=off
             -DNLOPT_SWIG=off)
add_library(NLopt::nlopt STATIC IMPORTED)
set_target_properties(
  NLopt::nlopt
  PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${3RDPARTY_DIR}/include"
             IMPORTED_LINK_INTERFACE_LANGUAGES "C;CXX"
             INTERFACE_LINK_LIBRARIES "m"
             IMPORTED_LOCATION "${3RDPARTY_DIR}/lib/libnlopt.a")

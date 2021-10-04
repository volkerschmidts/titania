set(3RDPARTY_DIR ${CMAKE_BINARY_DIR}/3rdparty)

set(BLA_STATIC ON)
include(ExternalProject)
ExternalProject_Add(
  project_cppunit
  GIT_REPOSITORY "https://anongit.freedesktop.org/git/libreoffice/cppunit.git"
  GIT_TAG "cppunit-1.15.1"
  GIT_SHALLOW ON
  PREFIX ${3RDPARTY_DIR}
  INSTALL_DIR ${3RDPARTY_DIR}
  BUILD_IN_SOURCE ON
  CONFIGURE_COMMAND ./autogen.sh
  COMMAND ./configure --disable-doxygen --prefix=<INSTALL_DIR>)
add_library(cppunit STATIC IMPORTED)
set(CPPUNIT_INCLUDE_DIRS "${3RDPARTY_DIR}/include")
set(CPPUNIT_LIBRARIES "${3RDPARTY_DIR}/lib/libcppunit.a")

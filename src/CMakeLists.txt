add_executable(${PROJECT_NAME} Main.cpp)

add_library(titania_lib)
target_sources(
  titania_lib
  PRIVATE Atom.cpp
          Bond.cpp
          DIDC.cpp
          Eckart.cpp
          Fragment.cpp
          Helper.cpp
          IndependencyAnalysis.cpp
          LinAlg.cpp
          Molecule.cpp
          MonteCarloStatistic.cpp
          Output.cpp
          Potential.cpp
          Properties.cpp
          RDCdata.cpp
          RDCset.cpp
          RedundantInternals.cpp
          SCRM.cpp
          Simplex.cpp
          ShoemakeEulerAngles.cpp
          SphericalHarmonics.cpp
          SphericalHarmonicsMatrix.cpp
          Structure.cpp
          StructureSimulator.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE titania_lib)

target_include_directories(
  ${PROJECT_NAME}
  PRIVATE $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
          $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)

target_include_directories(
  titania_lib
  PRIVATE $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
          $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
target_link_libraries(${PROJECT_NAME} PRIVATE project_options project_warnings)
target_link_libraries(titania_lib PRIVATE project_options project_warnings)

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
  target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(titania_lib PRIVATE OpenMP::OpenMP_CXX)
endif()

if(${PROJECT_NAME}_BUILD_3RDPARTY_EIGEN3)
  # HACK to support interface libraries pointing to not-yet-existing paths see
  # https://gitlab.kitware.com/cmake/cmake/-/issues/15052
  file(MAKE_DIRECTORY ${3RDPARTY_DIR}/include/eigen3)
  add_dependencies(titania_lib project_eigen3)
  target_include_directories(${PROJECT_NAME}
                             PUBLIC $<BUILD_INTERFACE:${3RDPARTY_DIR}/include>)
else()
  find_package(
    Eigen3
    3.3
    NO_MODULE
    REQUIRED)
endif()
target_link_libraries(titania_lib PRIVATE Eigen3::Eigen)

if(${PROJECT_NAME}_BUILD_3RDPARTY_LAPACKE)
  # HACK to support interface libraries pointing to not-yet-existing paths see
  # https://gitlab.kitware.com/cmake/cmake/-/issues/15052
  file(MAKE_DIRECTORY ${3RDPARTY_DIR}/include/openblas)
  add_dependencies(titania_lib project_lapacke)
  target_include_directories(
    titania_lib PRIVATE $<BUILD_INTERFACE:${3RDPARTY_DIR}/include/openblas>)
  target_link_libraries(titania_lib PRIVATE LAPACKE::LAPACKE)
else()
  find_package(BLAS REQUIRED)
  target_compile_options(titania_lib PRIVATE ${BLAS_LINKER_FLAGS})
  target_link_libraries(titania_lib PRIVATE ${BLAS_LIBRARIES})

  find_package(LAPACK REQUIRED)
  target_compile_options(titania_lib PRIVATE ${LAPACK_LINKER_FLAGS})
  # TODO this only brings in the default lapack libraries, but not LAPACKE, so we add it
  # manually
  target_link_libraries(titania_lib PRIVATE lapacke)
  target_link_libraries(titania_lib PRIVATE ${LAPACK_LIBRARIES})
endif()
if(BLA_STATIC)
  target_link_libraries(titania_lib PRIVATE gfortran)
endif()

if(${PROJECT_NAME}_BUILD_3RDPARTY_NLOPT)
  # HACK to support interface libraries pointing to not-yet-existing paths see
  # https://gitlab.kitware.com/cmake/cmake/-/issues/15052
  file(MAKE_DIRECTORY ${3RDPARTY_DIR}/include)
  add_dependencies(titania_lib project_nlopt)
  target_include_directories(titania_lib
                             PRIVATE $<BUILD_INTERFACE:${3RDPARTY_DIR}/include>)
else()
  find_package(NLopt REQUIRED)
endif()
target_link_libraries(titania_lib PRIVATE NLopt::nlopt)

option(${PROJECT_NAME}_BUILD_LIBPHOBOS "Enable compilation of libphobos plugin" ON)
if(${PROJECT_NAME}_BUILD_LIBPHOBOS)
  add_subdirectory(${CMAKE_SOURCE_DIR}/plugins/phobos phobos)
  target_include_directories(
    ${PROJECT_NAME}
    PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/plugins/phobos/include>)
  target_include_directories(
    titania_lib PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/plugins/phobos/include>)
  if(${PROJECT_NAME}_BUILD_3RDPARTY_LAPACKE)
    target_include_directories(
      phobos PRIVATE $<BUILD_INTERFACE:${3RDPARTY_DIR}/include/openblas>)
  endif()
  target_link_libraries(${PROJECT_NAME} PRIVATE phobos)
endif()

option(${PROJECT_NAME}_BUILD_PARSER "Enable compilation of Parser plugin" ON)
if(${PROJECT_NAME}_BUILD_PARSER)
  add_subdirectory(${CMAKE_SOURCE_DIR}/Parser Parser)
  target_link_libraries(${PROJECT_NAME} PRIVATE Parser)
  if(${PROJECT_NAME}_BUILD_3RDPARTY_EIGEN3)
    target_include_directories(Parser
                               PRIVATE $<BUILD_INTERFACE:${3RDPARTY_DIR}/include>)
  endif()
endif()

option(${PROJECT_NAME}_USE_CUDA "Use Nvidia's CUDA libraries" OFF)
if(${PROJECT_NAME}_USE_CUDA)
  add_subdirectory(${CMAKE_SOURCE_DIR}/plugins/cuda titania_cuda)
  target_compile_definitions(${PROJECT_NAME} PUBLIC USE_CUDA)
  target_include_directories(
    ${PROJECT_NAME} PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/plugins/cuda/include>)
  target_link_libraries(
    ${PROJECT_NAME}
    PRIVATE titania_cuda
            cudart
            cusolver
            cublas)
endif()

target_link_libraries(${PROJECT_NAME} PRIVATE m)

include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME} RUNTIME)

include(Packaging)

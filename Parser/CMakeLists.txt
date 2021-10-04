cmake_minimum_required(VERSION 3.15)

project(Parser CXX)

include(GNUInstallDirs)

add_library(
  ${PROJECT_NAME}
  KeyParser.cpp
  OwnHash.cpp
  Parser.cpp
  RDCParser.cpp
  StructureParser.cpp)

target_include_directories(
  ${PROJECT_NAME} PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
# target_include_directories(${PROJECT_NAME} PRIVATE src)

# install(TARGETS phobos LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
# PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

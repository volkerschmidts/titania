cmake_minimum_required(VERSION 3.15)

project(phobos C)

include(GNUInstallDirs)

add_library(${PROJECT_NAME} src/phobos.c)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER include/phobos.h)

target_include_directories(${PROJECT_NAME} PRIVATE include)
# target_include_directories(${PROJECT_NAME} PRIVATE src)

# install(TARGETS phobos LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER
# DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

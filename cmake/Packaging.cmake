include(InstallRequiredSystemLibraries)

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "${PROJECT_NAME} model-free interpretation of Residual Dipolar Couplings")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
set(CPACK_PACKAGE_VENDOR "Thielelab")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")
# set(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake
# ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
set(CPACK_STRIP_FILES "${CMAKE_INSTALL_BIN_DIR}/${PROJECT_NAME}")
set(CPACK_SOURCE_STRIP_FILES "")
set(CPACK_PACKAGE_EXECUTABLES ${PROJECT_NAME})
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY ON)

# CPack unfortunately packages the whole directory with **everything** in it and there
# is no clean way to only package the files tracked by git. There is
# https://github.com/RWTH-HPC/CMake-gitpack/ which filters files excluded by
# .gitattributes, but not all the other files potentially generated locally.
#
# include(CPack)

# This adds a custom (very much simplified) version of the build target which CPack
# would have created.
set(${PROJECT_NAME}_WITH_VERSION
    "${PROJECT_NAME}-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
)
set(PACKAGE_SOURCE_ARCHIVE_FILE
    ${CMAKE_BINARY_DIR}/${${PROJECT_NAME}_WITH_VERSION}-src.tar.gz)
add_custom_target(
  package_source
  COMMAND git archive --format=tar.gz --output ${PACKAGE_SOURCE_ARCHIVE_FILE}
          --prefix="${${PROJECT_NAME}_WITH_VERSION}/" HEAD
  COMMENT "-- Packaging Source Distribution to ${PACKAGE_SOURCE_ARCHIVE_FILE}"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

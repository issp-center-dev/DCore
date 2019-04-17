# Version number of the application
set (DCORE_VERSION "2.0")
set (DCORE_RELEASE "2.0.0")

# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# start configuration 
cmake_minimum_required(VERSION 2.8)
project(DCore NONE)
set(CMAKE_BUILD_TYPE Release)
enable_testing()

find_package(PythonInterp REQUIRED)
message(STATUS "Found python executable: ${PYTHON_EXECUTABLE}")

set(PYTHON_INTERPRETER ${PYTHON_EXECUTABLE})

#
# TO DO
#   * Git hash
#

message(STATUS "Installing DCore into ${CMAKE_INSTALL_PREFIX}...")

# Check version of TRIQS
execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cmake/show_triqs_version.py OUTPUT_VARIABLE TRIQS_VERSION)
if (TRIQS_VERSION STREQUAL "1.4" OR TRIQS_VERSION STREQUAL "2.1")
    message(STATUS "Found TRIQS Version ${TRIQS_VERSION}")
else()
    message(WARNING "Found TRIQS Version ${TRIQS_VERSION}, but may be imcompatible with DCore.")
endif()

# Find where the pytriqs is installed
execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cmake/show_triqs_site_packages.py OUTPUT_VARIABLE TRIQS_SITE_PACKAGES)

set(DCORE_SITE_PACKAGES "${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages")
set(DCORE_DEP_PYTHONPATH "${DCORE_SITE_PACKAGES}:${TRIQS_SITE_PACKAGES}")

add_subdirectory(python)
add_subdirectory(shells)
add_subdirectory(test)
add_subdirectory(tools)

option(BUILD_DOC "Build documentation" OFF)
if(${BUILD_DOC})
  add_subdirectory(doc)
endif(${BUILD_DOC})

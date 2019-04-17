# Version number of the application
set (DCORE_VERSION "2.0")
set (DCORE_RELEASE "2.0.0")

# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# start configuration 
cmake_minimum_required(VERSION 2.8)
#project(DCore CXX NONE)
project(DCore NONE)
set(CMAKE_BUILD_TYPE Release)
enable_testing()

find_package(PythonInterp REQUIRED)
message(STATUS "Found python executable: ${PYTHON_EXECUTABLE}")

set(PYTHON_INTERPRETER ${PYTHON_EXECUTABLE})

#
# TO DO
#   * Check version of TRIQS 
#   * Check version of Python
#   * Git hash
#   * Install version.py
#

message(STATUS "Installing DCore into ${CMAKE_INSTALL_PREFIX}...")

set(DCORE_SITE_PACKAGES "${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages")
#set(DCORE_DEP_PYTHONPATH "${DCORE_SITE_PACKAGES}:${TRIQS_PATH}/${TRIQS_PYTHON_LIB_DEST_ROOT}:${TRIQS_PATH}/libexec/venv/${TRIQS_PYTHON_LIB_DEST_ROOT}")
set(DCORE_DEP_PYTHONPATH "${DCORE_SITE_PACKAGES}")

add_subdirectory(python)
add_subdirectory(shells)
add_subdirectory(test)
add_subdirectory(tools)

option(BUILD_DOC "Build documentation" OFF)
if(${BUILD_DOC})
  add_subdirectory(doc)
endif(${BUILD_DOC})

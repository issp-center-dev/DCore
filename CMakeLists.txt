# Version number of the application
set (DCORE_VERSION "1.0")
set (DCORE_RELEASE "1.0.0")

# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

# start configuration 
cmake_minimum_required(VERSION 2.8)
project(DCore CXX NONE)
set(CMAKE_BUILD_TYPE Release)
enable_testing()

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

if (NOT ${TRIQS_WITH_PYTHON_SUPPORT})
 MESSAGE(FATAL_ERROR "dft_tools require Python support in TRIQS")
endif()

# Get hash
triqs_get_git_hash(${CMAKE_SOURCE_DIR} "DCORE")
if(${GIT_RESULT} EQUAL 0)
  message(STATUS "Hash: ${DCORE_GIT_HASH}")
endif(${GIT_RESULT} EQUAL 0)

# DCore will be installed in the TRIQS tree by default.
if(NOT CMAKE_INSTALL_PREFIX)
    message("Using default install prefix: ${TRIQS_PATH}")
    set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})
endif()

add_subdirectory(python)
add_subdirectory(shells)
add_subdirectory(test)
add_subdirectory(tools)
option(BUILD_DOC "Build documentation" OFF)
if(${BUILD_DOC})
  if(NOT TRIQS_WITH_DOCUMENTATION)
    message("Error: TRIQS library has not been compiled with its documentation")
  endif()
  add_subdirectory(doc)
endif(${BUILD_DOC})

# Version number of the application
set (DCORE_VERSION "2.0")
set (DCORE_RELEASE "2.0.0")

# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

option(BUILD_DCORE "Build DCore" ON)

# start configuration 
cmake_minimum_required(VERSION 2.8)
project(DCore NONE)
set(CMAKE_BUILD_TYPE Release)
enable_testing()

find_package(PythonInterp REQUIRED)
message(STATUS "Found python executable: ${PYTHON_EXECUTABLE}")
if (${PYTHON_VERSION_MAJOR} EQUAL 3)
    message(FATAL_ERROR "Python 3.x is not supported!")
endif()

set(PYTHON_INTERPRETER ${PYTHON_EXECUTABLE})

if (BUILD_DCORE)
    #
    # TO DO
    #   * Git hash
    #

    message(STATUS "Installing DCore into ${CMAKE_INSTALL_PREFIX}...")

    # Add pytriqs to PYTHONPATH
    if (NOT TRIQS_PATH)
        message(FATAL_ERROR "Set TRIQS_PATH to the install directory of TRIQS (e.g. /opt/triqs).")
    endif()
    set(TRIQS_SITE_PACKAGES ${TRIQS_PATH}/lib/python2.7/site-packages)
    set(ENV{PYTHONPATH} ${TRIQS_SITE_PACKAGES}:$ENV{PYTHONPATH})

    # Check version of TRIQS
    execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cmake/show_triqs_version.py OUTPUT_VARIABLE TRIQS_VERSION)
    if (TRIQS_VERSION STREQUAL "NOT_FOUND")
        message(FATAL_ERROR "TRIQS not found! Set TRIQS_PATH properly.")
    elseif (TRIQS_VERSION STREQUAL "1.4" OR TRIQS_VERSION MATCHES "^2.1.")
        message(STATUS "Found TRIQS Version ${TRIQS_VERSION}")
    else()
        message(WARNING "TRIQS Version ${TRIQS_VERSION} may be imcompatible with DCore.")
    endif()

    set(DCORE_SITE_PACKAGES "${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages")
    set(DCORE_DEP_PYTHONPATH "${DCORE_SITE_PACKAGES}:${TRIQS_SITE_PACKAGES}")

    # BSE
    if (BSE_PATH)
        set(ENV{PYTHONPATH} ${BSE_PATH}/python:$ENV{PYTHONPATH})
        set(DCORE_DEP_PYTHONPATH "${DCORE_DEP_PYTHONPATH}:${BSE_PATH}/python")
    endif()

    add_subdirectory(python)
    add_subdirectory(shells)
    add_subdirectory(test)
    add_subdirectory(tools)
endif()

option(BUILD_DOC "Build documentation" OFF)
if(${BUILD_DOC})
  add_subdirectory(doc)
endif(${BUILD_DOC})

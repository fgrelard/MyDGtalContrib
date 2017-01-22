cmake_minimum_required(VERSION 2.8.1)

project(MyDGtalContrib)
find_package(DGtal REQUIRED)
include_directories(${DGTAL_INCLUDE_DIRS})
link_directories(${DGTAL_LIBRARY_DIRS})

set(TEST_OUTPUT_PATH ${CMAKE_BINARY_DIR}/tests)
set(EXAMPLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/examples)
#set(CMAKE_BUILD_TYPE Release)

FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
endif( Boost_FOUND )


include_directories(${CMAKE_SOURCE_DIR}/src/)
option(BUILD_TESTS "Build test programs" OFF)
option(BUILD_EXAMPLES "Build examples programs" ON)
option(WITH_OPENMP "Parallelization with OpenMP" OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if (WITH_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif(WITH_OPENMP)


if (NOT WITH_EIGEN)
  message(FATAL_ERROR "Please compile DGtal with Eigen")
endif(NOT WITH_EIGEN)

if (NOT WITH_VISU3D_QGLVIEWER)
  message(FATAL_ERROR "Please compile DGtal with QGLViewer")
endif(NOT WITH_VISU3D_QGLVIEWER)

if (BUILD_TESTS)
  add_subdirectory(tests)
endif(BUILD_TESTS)

if (BUILD_EXAMPLES)
  add_subdirectory(examples)
endif(BUILD_EXAMPLES)




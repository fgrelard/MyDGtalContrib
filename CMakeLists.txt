cmake_minimum_required(VERSION 2.8.1)

project(MyDGtalContrib)
find_package(DGtal REQUIRED)
include_directories(${DGTAL_INCLUDE_DIRS})
link_directories(${DGTAL_LIBRARY_DIRS})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(TEST_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin/tests)
#set(CMAKE_BUILD_TYPE Release)

FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
endif( Boost_FOUND )


include_directories(${CMAKE_SOURCE_DIR}/src/)
option(BUILD_TESTS "Build test programs" OFF)
option(BUILD_EXAMPLES "Build examples programs" ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


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
  #add_subdirectory(examples)
endif(BUILD_EXAMPLES)




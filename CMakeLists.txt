cmake_minimum_required(VERSION 2.8.1)

project(MyDGtalContrib)
find_package(DGtal REQUIRED)
include_directories(${DGTAL_INCLUDE_DIRS})
link_directories(${DGTAL_LIBRARY_DIRS})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(TEST_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin/tests)
set(CMAKE_BUILD_TYPE Release)

FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
endif( Boost_FOUND )


include_directories(${CMAKE_SOURCE_DIR}/src/utils)
option(BUILD_TESTS "Build test programs" OFF)
option(WITH_CXX11 "Use C++11 standard" ON)

if (NOT WITH_CXX11)
   message(FATAL_ERRROR "Please compile DGtal using C++11")
endif(NOT WITH_CXX11)

if (NOT WITH_EIGEN)
  message(FATAL_ERROR "Please compile DGtal with Eigen")
endif(NOT WITH_EIGEN)

if (BUILD_TESTS)
  add_subdirectory(tests)
endif(BUILD_TESTS)




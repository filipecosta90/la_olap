#Minimum version of CMake to build this project
cmake_minimum_required (VERSION 2.8)

#Declaration of the project
project (olap)

#Compilation flags
# first set "normal flags"
set(CMAKE_CXX_FLAGS "-O3 -Wall -W -pedantic")
# next check for C++11 support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


  add_subdirectory(src/)
  #find_package(CUDA)

#Defines subdirectory based on cuda support
#if (CUDA_FOUND)
#  message("CUDA found, using gpu code!")
#  LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)

#else()
#  message("CUDA not found, not implemented yet!! Call filipe!")
#endif()




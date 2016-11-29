#Declaration of the project
project (olap_engine)

#Minimum version of CMake to build this project
cmake_minimum_required (VERSION 3.0)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(GLIB2 REQUIRED)

include_directories(${GLIB2_INCLUDE_DIRS})
link_directories (${GLIB2_LIBRARY_DIRS})

#Compilation flags
# first set "normal flags"
set(CMAKE_CXX_FLAGS "-O3 -Wall -W -pedantic")
# next check for C++11 support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

if(COMPILER_SUPPORTS_CXX11)
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};--std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++0x")
else()
  message(STATUS "The compiler ${CMAKE_CUDA_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

##set(CUDA_HOST_COMPILER clang++) # shouldn't fail anymore

#Defines subdirectory based on cuda support
find_package(CUDA)
if (CUDA_FOUND)
  message("CUDA found, using gpu code!")
  set(CUDA_SEPARABLE_COMPILATION ON)
  set(BUILD_SHARED_LIBS OFF)
  set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -std=c++11  -lineinfo -use_fast_math -Xptxas -dlcm=cg)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
  add_subdirectory(src/)
else()
  message("CUDA not found, not implemented yet!! Call filipe!")
endif()


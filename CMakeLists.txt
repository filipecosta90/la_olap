#Minimum version of CMake to build this project
cmake_minimum_required (VERSION 2.8)

#Declaration of the project
project (olap)

#Compilation flags
set(CMAKE_CXX_FLAGS "-O3 -Wall -W -pedantic")

#Defines subdirectory
add_subdirectory(src/)

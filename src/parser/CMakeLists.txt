include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

file(GLOB HH_FILES *.hh)
file(GLOB SRC_FILES *.cc)
file(GLOB CUDA_SRC_FILES *.cu )
file(GLOB CUDA_HH_FILES *.cuh )
file(GLOB HXX_FILES *.hxx)

find_package(BISON)
bison_target(Parser olap_parser.yy ${CMAKE_CURRENT_BINARY_DIR}/olap_parser.cc)

find_package(FLEX)
flex_target(Scanner olap_scanner.ll ${CMAKE_CURRENT_BINARY_DIR}/olap_scanner.cc)

add_flex_bison_dependency(Scanner Parser)

SET_SOURCE_FILES_PROPERTIES(olap_parser.hh GENERATED)


add_library(parser
  STATIC
  ${FLEX_Scanner_OUTPUTS}
  ${BISON_Parser_OUTPUTS}
  ${HH_FILES}
  ${SRC_FILES}
  )

CUDA_ADD_LIBRARY(cudalib
  ${CUDA_HH_FILES}
  ${CUDA_SRC_FILES}
  )


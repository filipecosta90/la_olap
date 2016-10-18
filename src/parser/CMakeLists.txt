include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

file(GLOB HH_FILES *.hh)
file(GLOB SRC_FILES *.cc)
file(GLOB CUDA_SRC_FILES *.cu )
file(GLOB HXX_FILES *.hxx)

find_package(BISON)
find_package(FLEX)

bison_target(Parser olap_parser.yy ${CMAKE_CURRENT_SOURCE_DIR}/olap_parser.cc
DEFINES_FILE ${CMAKE_CURRENT_SOURCE_DIR}/olap_parser.hh)
flex_target(Scanner olap_scanner.ll ${CMAKE_CURRENT_SOURCE_DIR}/olap_scanner.cc)

add_flex_bison_dependency(Scanner Parser)

SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/olap_parser.hh GENERATED)

CUDA_COMPILE(CUDA_O ${CUDA_SRC_FILES})

cuda_add_library(parser
  STATIC
 ${HH_FILES}
  ${SRC_FILES}
  ${BISON_Parser_OUTPUTS}
  ${BISON_Parser_OUTPUT_HEADER}
   ${FLEX_Scanner_OUTPUTS}
  ${HXX_FILES}
  ${CUDA_O}
  )


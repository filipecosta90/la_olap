# Add subdirectories
add_subdirectory(parser)

# Add the binary and sources

#CUDA_COMPILE(CUDA_O  parser/${CUDA_SRC_FILES})

CUDA_ADD_EXECUTABLE(
  ../../olap_engine
  main.cc
  ${CUDA_O}
  )

target_link_libraries(../../olap_engine parser)

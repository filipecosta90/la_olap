# Add subdirectories
add_subdirectory(parser)

# Add the binary and sources
CUDA_ADD_EXECUTABLE(
  ../../olap_engine
  main.cu
  )

target_link_libraries(../../olap_engine parser)

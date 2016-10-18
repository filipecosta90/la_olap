# Add subdirectories
add_subdirectory(parser)

# Add the binary and sources
ADD_EXECUTABLE(
  ../../olap_engine
  main.cc
  )

target_link_libraries(../../olap_engine parser)

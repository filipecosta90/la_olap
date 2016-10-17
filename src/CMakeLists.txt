# Add subdirectories

add_subdirectory(parser)

# Add the binary and sources
add_executable(
  ../../olap_engine
  main.cpp
  )

target_link_libraries(../../olap_engine parser)

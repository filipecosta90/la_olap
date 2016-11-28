# Add subdirectories
add_subdirectory(parser)

# Add the binary and sources

add_executable(
  olap_engine
  main.cc
  )

target_link_libraries(olap_engine parser)
target_link_libraries(olap_engine cudalib)
target_link_libraries (olap_engine ${GLIB2_LIBRARIES})

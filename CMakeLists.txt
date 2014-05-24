project(evsel)

cmake_minimum_required(VERSION 2.8)

add_subdirectory(evgraph)

INSTALL(FILES bshm/evsel.bshm DESTINATION share/evsel)

cmake_minimum_required(VERSION 3.10)
project(SynCore)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -g -Wall -Wno-unused-variable -pthread")
set(CMAKE_CXX_FLAGS_RELEASE   "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -O3 -pthread -fopenmp")

#add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

add_executable(MLCS_baseline MLCS_baseline.cpp
        utilities/multilayer_graph.cpp
        utilities/multilayer_graph.h
        utilities/projected_graph.cpp
        utilities/projected_graph.h)

add_executable(MLCS_prune MLCS_prune.cpp
        utilities/multilayer_graph.cpp
        utilities/multilayer_graph.h
        utilities/projected_graph.cpp
        utilities/projected_graph.h
        utilities/search_space.cpp
        utilities/search_space.h)

add_executable(MLCS_decomposition MLCS_decomposition.cpp
        utilities/multilayer_graph.cpp
        utilities/multilayer_graph.h
        utilities/projected_graph.cpp
        utilities/projected_graph.h
        utilities/search_space.cpp
        utilities/search_space.h)

add_executable(MLCS_index MLCS_index.cpp
        utilities/multilayer_graph.cpp
        utilities/multilayer_graph.h
        utilities/projected_graph.cpp
        utilities/projected_graph.h
        utilities/search_space.cpp
        utilities/search_space.h
        utilities/trie.cpp
        utilities/trie.h)

add_executable(MLCS_maintenance MLCS_maintenance.cpp
        utilities/multilayer_graph.cpp
        utilities/multilayer_graph.h
        utilities/projected_graph.cpp
        utilities/projected_graph.h
        utilities/search_space.cpp
        utilities/search_space.h
        utilities/trie.cpp
        utilities/trie.h)


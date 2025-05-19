cmake_minimum_required(VERSION 3.16)
project(Benchmark_researching)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -O0 -g -ggdb3")
# set(CMAKE_CXX_FLAGS "-Wall -O3")

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(test pgm_index_test.cpp)


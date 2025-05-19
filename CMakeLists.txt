cmake_minimum_required(VERSION 3.10)
project(Benchmark_researching)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -O0 -g -ggdb3")
# set(CMAKE_CXX_FLAGS "-Wall -O3")

include_directories(${PROJECT_SOURCE_DIR}/include)

# 添加OpenMP支持
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(test pgm_index_test.cpp)
add_executable(test2 pgm_index_test_2.cpp)
add_executable(ft fitting_tree_test.cpp)
target_include_directories(ft PRIVATE include)

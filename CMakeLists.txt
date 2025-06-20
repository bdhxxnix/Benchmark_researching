cmake_minimum_required(VERSION 3.10)
project(Benchmark_researching)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -O0 -g -ggdb3")

# Add include directories to each target
set(PROJECT_INCLUDES
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/stx-btree-0.8.6/include)

# OpenMP Support
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(pgmtest pgm_index_test.cpp)
target_include_directories(pgmtest PRIVATE ${PROJECT_INCLUDES})

add_executable(ltest linear_model_test.cpp)
target_include_directories(ltest PRIVATE ${PROJECT_INCLUDES})

add_executable(ftest fitting_tree_test.cpp)
target_include_directories(ftest PRIVATE ${PROJECT_INCLUDES})

add_executable(threads ltestForThreads.cpp)
target_include_directories(threads PRIVATE ${PROJECT_INCLUDES})





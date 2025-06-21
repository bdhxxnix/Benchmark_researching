#!/bin/bash

# build.sh - Build script for PLA Benchmark Project
# Author: Li Hua
# Usage: Run this script from the root of your repo: ./build.sh

set -e  # Exit immediately if a command fails

# --------- Configuration ---------
BUILD_DIR="build"
CXX_FLAGS="-O0 -Wall -fopenmp"  # Enable OpenMP
NUM_THREADS=$(nproc)            # Use all cores

# --------- Build Process ---------
echo " Creating build directory: $BUILD_DIR"
mkdir -p $BUILD_DIR
cd $BUILD_DIR

echo " Configuring project with CMake..."
cmake .. -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS="$CXX_FLAGS"

echo " Building project with $NUM_THREADS threads..."
make -j$NUM_THREADS

echo " Build complete. Binaries are in: $BUILD_DIR"

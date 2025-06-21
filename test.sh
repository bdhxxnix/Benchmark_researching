#!/bin/bash

# test.sh - Interface to run PLA benchmark on a given dataset
# Usage:
#   ./test.sh ./data/osm_keys.txt

set -e

# -------- Argument Parsing --------
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <dataset_path>"
    exit 1
fi

DATASET_PATH="$1"

if [ ! -f "$DATASET_PATH" ]; then
    echo "Dataset file not found: $DATASET_PATH"
    exit 1
fi

# -------- Configuration --------
BIN_DIR="./build"
TARGETS=("ltest" "threads" "pgmtest" "ftest")

# -------- Check Binaries --------
for target in "${TARGETS[@]}"; do
    if [ ! -x "$BIN_DIR/$target" ]; then
        echo "Binary not found: $BIN_DIR/$target"
        echo "Build first using ./build.sh"
        exit 1
    fi
done

# -------- Run Benchmarks --------
echo "Running tests on dataset: $DATASET_PATH"

for target in "${TARGETS[@]}"; do
    echo "Running: $target"
    "$BIN_DIR/$target" "$DATASET_PATH"
    echo "Finished: $target"
    echo "------------------------"
done

echo "All tests completed."

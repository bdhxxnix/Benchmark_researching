# PLA Benchmark

This repository benchmarks various Piecewise Linear Approximation (PLA) algorithms like `GreedyPLA`, `OptimalPLA`, `SwingFilter` used in learned indexes (e.g. PGM-Index, FITing-Tree).

## Introduction

This project aims to systematically compare different ε-PLA algorithms under unified benchmarking settings. These algorithms are used for segmenting key-index mappings with an error bound ε, essential for learned index structures like:
- **PGM-Index**: hierarchical index using PLA at each level
- **FITing-Tree**: learned B+ Tree with PLA-based leaf segments

## Getting started
Clone the repo and install dependencies:
```bash
git clone https://github.com/bdhxxnix/PLABench.git
```

## Prepare datasets
We write the preparation scripts based on SOSD repository: [SOSD](https://github.com/learnedsystems/SOSD)
- First, please install Python and use pip to download necessary dependencies `numpy` and `scipy`.
- Second, you need to install the zstd software to decompress datasets:
```bash
sudo apt install zstd
```
- Third, run the scrip and wait for everything to be ready:
```bash
bash prepare_data.sh
```

## Compile

After preparing your dataset, build the project by running:

```bash
./build.sh
```
This will generate all test runners inside the `build/` directory.


# Run
Before running the tests, ensure your dataset is ready and placed in the correct directory.
Then run the test by:
```bash
./test.sh ./your_dataset/data
```
This script will automatically execute the following tests in order:
- Linear Test
- Linear Test with Different Threads
- PGM Test
- FIT Test

You can also navigate to the `build/` directory and execute individual test runners manually.
For example, to run only the PGM test:
```bash
cd build/
./pgmtest <dataset_path>
```
To capture all the outputs into a certain log file:
```bash
./test.sh ./your_dataset/data > output.log 2>&1
```
This will save everything printed to the terminal into output.log.

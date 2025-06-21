# PLA Benchmark

This repository benchmarks various Piecewise Linear Approximation (PLA) algorithms like `GreedyPLA`, `OptimalPLA`, `SwingFilter` used in learned indexes (e.g. PGM-Index, FITing-Tree).

# Features
- Parallel FRS and Greedy implementations
- Easy-to-run benchmarking scripts

# Getting started
## Clone the repo
Clone the repo and install dependencies:
```bash
git clone https://github.com/bdhxxnix/PLABench.git
```

## Prepare datasets
We write the preparation scripts based on SOSD repository: [text](https://github.com/learnedsystems/SOSD)
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

## Run


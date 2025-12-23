# Fast Matrix Multiplication Benchmark

## Overview

In this project I implemented and benchmarked three matrix multiplication algorithms:

- naive matrix multiplication
- Strassen algorithm
- AlphaEvolve (48-multiplication 4x4 kernel)

The goal was to compare their real performance and the real number of arithmetic operations, not only asymptotic complexity.

All experiments were run on Linux (Debian).  
Results are exported to CSV and visualized with Python.

Don't forget to look at the plots after reading this report.

---

## Command Line Interface

The benchmark is a single executable called `bench`.

Basic usage:

```./bench --algo <name> --sizes <list>```


Parameters:

- `--algo`  
  Algorithm to run. Possible values:
    - `naive`
    - `strassen`
    - `alphaevolve`

- `--sizes`  
  Comma-separated list of matrix sizes.  
  Example: `--sizes 4,128,256,512`


Optional parameters:

- `--runs`  
  Number of timed runs per size (default: 5)

- `--warmup`  
  Number of warmup runs before timing (default: 1)

The program prints results in CSV format to stdout.

---

## Algorithms

### Naive

The naive algorithm uses three nested loops.

Properties:
- Multiplications: n³
- Additions: n³
- Very simple memory access pattern
- Good cache locality

This algorithm is used as a baseline.

---

### Strassen

Strassen reduces the number of scalar multiplications by recursively splitting matrices.

Properties:
- Multiplications grow as n^(log₂7)
- Many more additions than naive
- Large temporary matrices
- Worse cache locality
- High memory traffic

In this project Strassen is implemented with a recursive algorithm and exact operation counting.

---

### AlphaEvolve

AlphaEvolve uses a fixed 4x4 kernel with 48 multiplications.

Properties:
- Exactly 48 multiplications per 4x4 block
- Many additions inside the kernel
- Implemented as a blocked algorithm for larger matrices
- Padding is used for sizes not divisible by 4

This implementation uses the exact formulas from the provided alpha48.c file.

---

## Correctness

Before running benchmarks, the program checks correctness using a fixed 4x4 test case.

The result of each algorithm is compared against the naive implementation.

If the test fails, the program stops.

---

## Results Summary

### Execution Time

- The naive algorithm is the fastest for all tested sizes.
- Strassen is slower than naive, even though it uses fewer multiplications.
- AlphaEvolve is the slowest.

This happens because:
- Strassen and AlphaEvolve use many more additions.
- They allocate and move much more data.
- Memory access dominates execution time.
- Modern CPUs handle simple loops very efficiently.

---

### Number of Operations

- Naive has the largest number of multiplications and additions.
- Strassen significantly reduces multiplications.
- AlphaEvolve reduces multiplications even further.

However:
- Fewer multiplications does not mean faster execution.
- Additions and memory access matter a lot.

---

## Conclusion

This project shows that asymptotic operation counts do not directly predict real performance.

On modern hardware:
- Memory access and cache behavior are critical.
- Simple algorithms can outperform theoretically faster ones.
- Reducing multiplications alone is not enough.

Strassen and AlphaEvolve are interesting from a theoretical point of view, but naive matrix multiplication is still the fastest in practice for the tested sizes.


# Chudnovsky + Binary Splitting: research-guided optimization patch

This package keeps the same core algorithm as your current fastest version and applies only lower-risk implementation optimizations.

## What changed

1. **Leaf-block backward summation**
   - Instead of recursing all the way down to single terms everywhere, small ranges are summed iteratively in one block.
   - This reduces recursion overhead and the number of `mpz_init/mpz_clear` hot-path calls.

2. **Cheaper merge step**
   - The merge `T = P1*T2 + T1*Q2` now uses `mpz_addmul`, which avoids one extra temporary in the hot path.

3. **OpenMP task granularity control**
   - Added `task_min_size` so that tiny subtrees do not become OpenMP tasks.
   - This usually improves task scheduling efficiency on laptop CPUs.

4. **Same output shape as your current code**
   - This keeps timing comparisons straightforward.

## Files

- `bs_single/pi_bs_single_research_opt.c`
- `bs_openmp/pi_bs_openmp_research_opt.c`

## Build

Single-thread:

```bash
gcc -O3 -march=native -std=c11 bs_single/pi_bs_single_research_opt.c -lmpfr -lgmp -o pi_bs_single_research_opt
```

OpenMP:

```bash
gcc -O3 -march=native -fopenmp -std=c11 bs_openmp/pi_bs_openmp_research_opt.c -lmpfr -lgmp -o pi_bs_openmp_research_opt
```

## Run

Single-thread:

```bash
./pi_bs_single_research_opt 1000000 64
```

Arguments:
- `digits`
- `leaf_size`

OpenMP:

```bash
export OMP_NUM_THREADS=8
./pi_bs_openmp_research_opt 1000000 4 64 4096
```

Arguments:
- `digits`
- `task_depth`
- `leaf_size`
- `task_min_size`

## Suggested tuning sweep on a laptop

Single-thread:
- `leaf_size`: 32, 64, 96, 128

OpenMP:
- `task_depth`: 3, 4, 5
- `leaf_size`: 64, 96, 128
- `task_min_size`: 1024, 2048, 4096, 8192

## What is intentionally NOT included here

The literature and y-cruncher internals suggest two larger optimizations that can still help, but they require a more invasive rewrite:

1. **Fully factored / GCD-factorized binary splitting**
2. **Top-level truncation / PiFast-style space reduction**

Those are promising, but they change the implementation much more substantially than this patch.

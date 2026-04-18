# y-cruncher-inspired Chudnovsky exact binary splitting

This package contains:

- `pi_bs_baseline_from_main`: a copy of the uploaded baseline program.
- `pi_bs_yc_inspired`: a conservative engineering rewrite inspired by y-cruncher-style ideas.

## What changed in `pi_bs_yc_inspired`

1. Keep the same core algorithm: Chudnovsky + exact binary splitting + OpenMP tasks.
2. Reuse per-thread leaf scratch (`lhs`, `rhs`) instead of `mpz_init/mpz_clear` in every leaf call.
3. Keep single-child task spawning for recursion.
4. Add *merge-task thresholds*: only higher tree levels parallelize the 4 expensive merge multiplications.
5. Keep a laptop-friendly execution model: exact arithmetic, bounded task creation, no risky truncation/GCD rewrites.

## Build

```bash
make
```

## Run baseline

```bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=close
export OMP_PLACES=cores
./pi_bs_baseline_from_main 10000000 5 128 4096
```

## Run y-cruncher-inspired version

```bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=close
export OMP_PLACES=cores
./pi_bs_yc_inspired 10000000 5 128 4096 2 65536
```

Arguments:

- `digits`
- `task_depth`
- `leaf_size`
- `task_min_size`
- `merge_task_depth_limit`
- `merge_task_min_size`

## Suggested starting points

10 million digits:

```bash
./pi_bs_yc_inspired 10000000 5 128 4096 2 65536
./pi_bs_yc_inspired 10000000 5 128 4096 1 131072
```

100 million digits:

```bash
./pi_bs_yc_inspired 100000000 6 256 4096 2 131072
./pi_bs_yc_inspired 100000000 6 256 8192 1 262144
```

## Notes

This is **not** y-cruncher source code.
It is a small C/OpenMP program that borrows a few practical ideas:

- multiplication-heavy hotspots matter most
- scratch reuse matters
- task decomposition matters more than simply adding more algorithmic variants

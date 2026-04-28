# pi_yc_learned_opt

This package contains:

- `pi_bs_baseline_from_user.c`: your uploaded 23s-class baseline.
- `pi_bs_optimized_yc_learned.c`: an optimized rewrite that keeps the same
  Chudnovsky + exact binary splitting + OpenMP family, but folds in several
  lessons inspired by how y-cruncher attacks this class of problem.
- `pi_bs_fixedpoint_yc_adaptive.c`: the recommended new version in this folder.
  It keeps the exact binary splitting core, but replaces the slow MPFR finish
  with a GMP-only fixed-point finish that is closer to how high-performance
  constant engines avoid general-purpose floating-point overhead.
- `pi_bs_gmp_kernel_tuned.c`: an experimental lower-level GMP-oriented variant.
  It adds root/right-subtree `QT` pruning and tries to reduce `mpz`
  reallocation overhead. On my current tests it did **not** beat the main
  fixed-point version, so it is kept as a research branch rather than the
  default recommendation.

## What changed in the optimized version

1. **Removed the mathematically-unused root `P` result**
   - The root only needs `Q` and `T` to evaluate `pi = 426880*sqrt(10005)*Q/T`.
   - The optimized version introduces a `QT` recursion variant that does **not**
     compute the aggregate `P` on the root/right-spine where it is never used.
   - This removes at least one full-size `P = P1 * P2` multiply and additional
     right-spine `P` work.

2. **Rewrote the leaf recurrence**
   - Old form:
     - `lhs = p*T`
     - `rhs = p*(A+Bk)*Q`
     - `T = lhs + rhs`
   - New exact form:
     - `T = p * (T +/- (A+Bk)*Q)`
   - This lets the leaf use `mpz_addmul_ui()` / `mpz_submul_ui()` and removes
     the `lhs/rhs` temporaries entirely.

3. **Parallel merge only where it is worth it**
   - The baseline always parallelized the 4 heavy merge multiplies on parallel
     nodes.
   - The optimized version adds two knobs:
     - `merge_task_depth`
     - `merge_task_min_size`
   - Above those thresholds, merge multiplies are parallelized; below them,
     they are done serially to avoid task overhead.

4. **Kept the good parts of the current baseline**
   - one-sided child tasking
   - exact binary splitting
   - OpenMP parallel recursion
   - balanced split

## What the new fixed-point version changes

1. **Leaf recurrence is rewritten in a tighter exact form**
   - Old baseline leaf:
     - builds explicit `lhs = p*T`
     - builds explicit `rhs = p*(A+Bk)*Q`
   - New fixed-point version:
     - updates `T = p * (T +/- (A+Bk) * Q)`
   - This removes per-leaf scratch temporaries and reduces big-integer traffic.

2. **Final evaluation no longer depends on MPFR**
   - Baseline ending:
     - convert giant `Q` and `T` into MPFR
     - compute `sqrt(10005)`
     - do floating-point multiply/divide at full precision
   - New ending:
     - compute `floor(sqrt(10005) * 10^(digits+guard))` with `mpz_sqrt`
     - compute `floor(426880 * Q * sqrt_scaled / T)` with integer division
   - This borrows a key lesson from y-cruncher:
     high-end pi code prefers exact/fixed-point pipelines and only converts to
     decimal output at the edge.
   - The program currently prints a **truncated** 50-digit preview from the
     fixed-point integer result. That preview may differ by 1 ulp from the
     MPFR-rounded preview printed by your original program, but the main target
     is lower endgame cost.

3. **Task fanout is adaptive instead of unconditional**
   - top-level merge multiplies can still run in parallel
   - deeper/smaller nodes fall back to serial merge to avoid task overhead

4. **Default task depth can auto-scale with core count**
   - passing `0` for `task_depth` lets the program pick a depth from
     `OMP_NUM_THREADS`

## What y-cruncher is doing that your GMP/MPFR code still cannot match

The logs and files inside `y-cruncher v0.8.7.9547-static/` show several
advantages that are outside the reach of a small `mpz_t/mpfr_t` rewrite:

- architecture-specific binaries and tuning profiles
  - your package contains separate binaries such as `14-BDW ~ Kurumi`,
    `22-ZN4 ~ Kizuna`, `24-ZN5 ~ Komari`
- a custom threaded runtime
  - the log shows `Threading Mode: Push Pool`
- a custom allocator
  - the log shows `Allocator : "mmap"` with multiple allocator threads
- specialized large multiply / divide / inverse-square-root kernels
  - the event log explicitly separates `Large Division`, `InvSqrt(10005)`,
    `Large Multiply`
- a different internal Chudnovsky representation
  - the config uses a hypergeometric form and the benchmark log reports
    `Series CommonP2B3`

So the realistic takeaway is:

- you can borrow algebraic pruning and fixed-point finishing
- you cannot expect GMP + MPFR + OpenMP C code to reach y-cruncher speed
  without dropping below the `mpz/mpfr` abstraction and writing/tuning
  low-level big-integer kernels

## Build on Ubuntu

Install dependencies:

```bash
sudo apt-get update
sudo apt-get install -y build-essential libgmp-dev libmpfr-dev
```

Build:

```bash
make
```

## Run

Recommended environment:

```bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=close
export OMP_PLACES=cores
```

Run the baseline:

```bash
./pi_bs_baseline_from_user 100000000
```

Run the older optimized version:

```bash
./pi_bs_optimized_yc_learned 100000000 12 256 4096 2 65536
```

Run the recommended fixed-point version:

```bash
./pi_bs_fixedpoint_yc_adaptive 100000000 7 256 4096 2 65536 32 1 0 1 0 0
```

Arguments for the fixed-point version:

1. `digits`
2. `task_depth`
   - use `0` for auto
3. `leaf_size`
4. `task_min_size`
5. `merge_task_depth`
6. `merge_task_min_size`
7. `guard_digits`
8. `print_pi`
9. `skew_num`
10. `skew_den`
11. `skew_depth`
12. `skew_min_size`

The recommended starting point above means:

- no skewed split
- `task_depth=7`
- `leaf_size=256`
- `task_min_size=4096`
- `merge_task_depth=2`
- `merge_task_min_size=65536`

Those were the best large-scale settings among the configurations I tested here.

## Experimental low-level GMP branch

Build and run:

```bash
./pi_bs_gmp_kernel_tuned 100000000 7 256 4096 2 65536 32 0 0 1 0 0
```

This version is for experimentation only. It tries to learn from the idea of
"specialized bigint kernels", but because GMP already uses `mpn`-level
optimized multiplication internally, these extra structural changes were not
enough to beat the simpler fixed-point version in my tests.

## Compare both

```bash
./run_compare.sh 100000000 0 256 4096 2 65536 32
```

## Tune on your laptop first at lower scale

```bash
./autotune.sh 10000000
```

This script now tunes `pi_bs_fixedpoint_yc_adaptive` rather than the older
`pi_bs_optimized_yc_learned`.

Then take the best command and rerun it at `100000000` digits.

# pi_yc_learned_opt

This package contains:

- `pi_bs_baseline_from_user.c`: your uploaded 23s-class baseline.
- `pi_bs_optimized_yc_learned.c`: an optimized rewrite that keeps the same
  Chudnovsky + exact binary splitting + OpenMP family, but folds in several
  lessons inspired by how y-cruncher attacks this class of problem.

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

Run the optimized version:

```bash
./pi_bs_optimized_yc_learned 100000000 12 256 4096 2 65536
```

Arguments for the optimized version:

1. `digits`
2. `task_depth`
3. `leaf_size`
4. `task_min_size`
5. `merge_task_depth`
6. `merge_task_min_size`

## Compare both

```bash
./run_compare.sh 100000000 12 256 4096 2 65536
```

## Tune on your laptop first at lower scale

```bash
./autotune.sh 10000000
```

Then take the best command and rerun it at `100000000` digits.

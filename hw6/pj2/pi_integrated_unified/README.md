# Unified Chudnovsky + Exact Binary Splitting Workspace

This package keeps the fast exact OpenMP binary-splitting path as the default, and integrates the *safe* optimizations that can coexist in one executable:

- OpenMP task parallelism
- leaf blocking with reverse accumulation
- task cutoff (`task_min_size`)
- optional top-level skewed splitting
- optional Newton reciprocal endgame

## Why not include naive GCD reduction here?
A local `gcd(Q, T)` reduction inside a normal exact BS tree is **not** algebraically safe unless the entire representation is redesigned around factorized forms. The earlier experimental version that did naive reductions can therefore break correctness. A true factorized-BS implementation is a larger rewrite.

## Build on Ubuntu

```bash
sudo apt-get update
sudo apt-get install -y build-essential libgmp-dev libmpfr-dev
make
```

## Run baseline-equivalent exact mode

```bash
export OMP_NUM_THREADS=8
./pi_bs_unified_exact 100000000 6 256 4096 0 1 0 0 0
```

Arguments:

1. `digits`
2. `task_depth`
3. `leaf_size`
4. `task_min_size`
5. `skew_num`
6. `skew_den`
7. `skew_depth`
8. `skew_min_size`
9. `endgame_mode` (`0` = MPFR direct divide, `1` = Newton reciprocal)

## Try light top-skew only near the root

```bash
./pi_bs_unified_exact 100000000 6 256 4096 2 3 1 131072 0
./pi_bs_unified_exact 100000000 6 256 4096 2 3 2 131072 0
```

## Try Newton endgame

```bash
./pi_bs_unified_exact 100000000 6 256 4096 0 1 0 0 1
```

## Autotune on your laptop

```bash
export OMP_NUM_THREADS=8
./autotune.sh 5000000
```

Then rerun the winning command at your target digit count.

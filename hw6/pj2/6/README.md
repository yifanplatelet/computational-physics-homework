# Strengthened Pi Program Workspace

This package contains two OpenMP programs for Chudnovsky + exact binary splitting:

- `pi_bs_openmp_baseline`: the simple exact-BS baseline.
- `pi_bs_openmp_factorized`: an enhanced experimental version inspired by `ilmp`.

## What is added in the factorized version

Compared with the baseline, the factorized version keeps:

- exact binary splitting
- OpenMP tasks
- backward block leaves

and additionally tries a safe version of the `ilmp` idea:

- maintain sparse prime-factor lists for subtree `P` and `Q`
- at selected internal levels, cancel common factors between **left `P`** and **right `Q`** before merging
- then merge the reduced subtrees exactly

This is **not** a full port of `ilmp`'s custom bigint engine. It still uses GMP/MPFR.

## Build on Ubuntu

```bash
sudo apt-get update
sudo apt-get install -y build-essential libgmp-dev libmpfr-dev
make
```

## Run

Set threads first:

```bash
export OMP_NUM_THREADS=8
```

### Baseline

```bash
./pi_bs_openmp_baseline 10000000 5 128 4096
```

Arguments:

1. `digits`
2. `task_depth`
3. `leaf_size`
4. `task_min_size`

### Factorized

```bash
./pi_bs_openmp_factorized 10000000 5 128 4096 4 8192
```

Arguments:

1. `digits`
2. `task_depth`
3. `leaf_size`
4. `task_min_size`
5. `cancel_level`
6. `cancel_min_size`

You can disable factorized cancellation with:

```bash
./pi_bs_openmp_factorized 10000000 5 128 4096 999 999999999
```

## Suggested experiments

Start with the baseline you already know is good, then try:

```bash
./pi_bs_openmp_factorized 10000000 5 128 4096 4 8192
./pi_bs_openmp_factorized 10000000 5 128 4096 5 16384
./pi_bs_openmp_factorized 10000000 6 256 4096 4 32768
```

For 1e8 digits, try the same shape around your current best settings.

## Notes

- This factorized version is a research-style strengthened prototype.
- It is meant to test whether `ilmp`-style cross-cancellation helps on your machine.
- Because the environment used to generate this package did not have MPFR headers available, the code was prepared carefully but could not be fully compiled here. If your compiler reports any issue, share the error output and it can be patched quickly.

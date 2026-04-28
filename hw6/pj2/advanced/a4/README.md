# advanced/a4

这一版只替换了末端求值：

- binary splitting 沿用最快的一条 OpenMP + GMP 主线
- `terms` 改成更准确的 Chudnovsky 项数估计
- final evaluation 从 MPFR 改为 GMP `mpf_t`

## 编译

```bash
make
```

## 运行

```bash
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_NUM_THREADS=8
./pi_bs_extreme_mpf 100000000 12 128 4096
```

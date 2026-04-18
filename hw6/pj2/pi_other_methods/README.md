# pi_other_methods

这一包包含 5 个可直接在 Ubuntu 上对比的实现：

- `pi_chudnovsky_openmp_baseline`：你当前主线的 Chudnovsky + OpenMP + exact binary splitting 基线
- `pi_ramanujan_openmp_bs`：Ramanujan 级数 + OpenMP + exact binary splitting
- `pi_gauss_legendre`：Gauss-Legendre / Brent-Salamin AGM 算法
- `pi_borwein_quartic`：Borwein 四次收敛算法
- `pi_machin_mpfr`：Machin 公式 + MPFR arctan 级数

## Ubuntu 依赖

```bash
sudo apt-get update
sudo apt-get install -y build-essential libgmp-dev libmpfr-dev
```

## 编译

```bash
make
```

## 运行

### 1) Chudnovsky 基线

```bash
export OMP_NUM_THREADS=8
./pi_chudnovsky_openmp_baseline 10000000 5 128 4096
```

参数：
- `digits`
- `task_depth`
- `leaf_size`
- `task_min_size`

### 2) Ramanujan + exact BS

```bash
export OMP_NUM_THREADS=8
./pi_ramanujan_openmp_bs 10000000 5 128 4096
```

参数同上。

### 3) Gauss-Legendre

```bash
./pi_gauss_legendre 10000000
./pi_gauss_legendre 10000000 28
```

参数：
- `digits`
- `iters`，省略或传负数表示自动估计

### 4) Borwein quartic

```bash
./pi_borwein_quartic 10000000
./pi_borwein_quartic 10000000 14
```

参数：
- `digits`
- `iters`，省略或传负数表示自动估计

### 5) Machin 公式

```bash
./pi_machin_mpfr 1000000
```

Machin 公式在高位数下一般会明显慢于前几种，仅适合当对照组。

## 建议测试顺序

先跑：

```bash
export OMP_NUM_THREADS=8
./pi_chudnovsky_openmp_baseline 10000000 5 128 4096
./pi_ramanujan_openmp_bs 10000000 5 128 4096
./pi_gauss_legendre 10000000
./pi_borwein_quartic 10000000
```

如果你要上 1e8，建议先只测：
- Chudnovsky 基线
- Ramanujan + exact BS
- Gauss-Legendre
- Borwein quartic

Machin 不建议在 1e8 上跑。

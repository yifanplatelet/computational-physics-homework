# pi_bs_openmp_rewrite

这是一份保守重写版：不再继续叠加那些已经被实测证明更慢的变体，而是把仍然“有赢面”的东西收敛到一条主线上：

- Chudnovsky
- exact binary splitting
- OpenMP tasks
- 叶子块反向累加
- 单侧 task 生成（只 spawn 一个子树，另一个在当前线程直接递归）
- 叶子层尽量改用 `mpz_mul_ui` 链，减少构造中间 `p/q/t` 大整数的开销

## Ubuntu 编译

```bash
sudo apt-get update
sudo apt-get install -y build-essential libgmp-dev libmpfr-dev
make
```

## 运行

建议先设线程绑定：

```bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=close
export OMP_PLACES=cores
```

然后运行：

```bash
./pi_bs_openmp_rewrite 10000000 5 128 4096
```

参数依次为：

1. `digits`
2. `task_depth`
3. `leaf_size`
4. `task_min_size`

## 自动调参

先在 1e7 或 5e6 上扫参数，再把最优参数拿去跑 1e8：

```bash
./autotune.sh 10000000
```

## 推荐起点

你现在测出来的好参数附近，建议先从这些开始：

```bash
./pi_bs_openmp_rewrite 10000000 5 128 4096
./pi_bs_openmp_rewrite 10000000 5 256 4096
./pi_bs_openmp_rewrite 10000000 6 128 4096
./pi_bs_openmp_rewrite 100000000 6 256 4096
```

## 关于 y-cruncher

这份程序和 y-cruncher **属于同一算法家族**，但**不是** y-cruncher 本体，也还没有实现 y-cruncher 那些更重的工程优化（例如更完整的精度控制、skewed splitting/backwards summing 的成熟版本、专用大整数内核、磁盘外存模式等）。

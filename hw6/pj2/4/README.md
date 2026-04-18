# 优化版说明

这版在你当前的 Chudnovsky + binary splitting 基础上，主要做了两个更实际的优化：

1. **叶子块化**
   - 原版在整个区间上不断递归，直到最小区间。
   - 这版在小区间内改为顺序合并，减少函数调用和 `mpz_init/mpz_clear` 开销。

2. **叶子大小可调**
   - 单线程版：`./pi_bs_single_opt [digits] [leaf_size]`
   - OpenMP 版：`./pi_bs_openmp_opt [digits] [task_depth] [leaf_size]`
   - 建议测试 `leaf_size = 16, 32, 64, 128`，选机器上最快的值。

## 默认参数建议

- 单线程：`leaf_size=32`
- OpenMP：`task_depth=4`, `leaf_size=64`

## 编译

```bash
gcc -O3 -march=native -std=c11 bs_single/pi_bs_single_opt.c -o pi_bs_single_opt -lmpfr -lgmp
gcc -O3 -march=native -fopenmp -std=c11 bs_openmp/pi_bs_openmp_opt.c -o pi_bs_openmp_opt -lmpfr -lgmp
```

## 运行示例

```bash
./pi_bs_single_opt 1000000 32
export OMP_NUM_THREADS=8
./pi_bs_openmp_opt 1000000 4 64
```

## 你最该调的参数

- `task_depth`：并行任务深度
- `leaf_size`：叶子块大小

通常先固定线程数，再扫一遍：

- `task_depth`: 3, 4, 5, 6
- `leaf_size`: 16, 32, 64, 128

选总时间最短的一组。

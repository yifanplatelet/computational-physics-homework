# GCD-reduced Binary Splitting

这是在你原始 `Chudnovsky + exact binary splitting` 基线上的**保守版 GCD 约分**实验。

## 做了什么

- 保持原来的中点二分和 exact `P/Q/T` 表示。
- 在每个叶子和每次合并后，执行一次
  - `g = gcd(Q, T)`
  - `Q /= g`
  - `T /= g`
- 这样不会改变数学结果，但有机会减小部分中间数。

这不是 fully factored binary splitting；它只是一个容易编译和 benchmark 的安全版本。

## Ubuntu 编译

```bash
sudo apt-get update
sudo apt-get install -y build-essential libgmp-dev libmpfr-dev
gcc -O3 -march=native -std=c11 pi_bs_gcd_reduced.c -lmpfr -lgmp -o pi_bs_gcd_reduced
```

## 运行

```bash
./pi_bs_gcd_reduced 1000000
```

参数：
- 第一个参数：目标小数位数 `digits`

## 输出说明

- `gcd_reductions`：实际发生约分的次数
- `binary_splitting_time`：exact binary splitting 时间
- `final_eval_time`：最后 `sqrt(10005)` 和一次高精度除法时间

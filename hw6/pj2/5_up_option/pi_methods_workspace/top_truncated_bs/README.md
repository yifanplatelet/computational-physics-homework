# Top-truncated / Hybrid Binary Splitting

这是一个**顶层截断风格**的实验版：

- 下层子树仍然用 exact GMP `P/Q/T` binary splitting
- 上层若干层不再保留完整大整数，而是转成 MPFR 的
  - `r = P / Q`
  - `s = T / Q`
- 合并公式变成：
  - `s = s1 + r1 * s2`
  - `r = r1 * r2`
- 最终 `pi = (426880 * sqrt(10005)) / s`

这更接近“顶层截断 / PiFast 风格”的方向：
上层不再让中间整数无限长大，而是直接在目标精度上做近似合并。

## Ubuntu 编译

```bash
sudo apt-get update
sudo apt-get install -y build-essential libgmp-dev libmpfr-dev
gcc -O3 -march=native -std=c11 pi_bs_top_truncated.c -lmpfr -lgmp -o pi_bs_top_truncated
```

## 运行

```bash
./pi_bs_top_truncated 1000000 3 64
```

参数：
- 第一个参数：`digits`
- 第二个参数：`approx_depth`，顶层用 MPFR 近似合并的深度
- 第三个参数：`exact_leaf`，下层 exact 子树切换阈值

## 调参建议

可试：

```bash
./pi_bs_top_truncated 1000000 2 64
./pi_bs_top_truncated 1000000 3 64
./pi_bs_top_truncated 1000000 4 128
```

这个版本是实验型混合实现，重点看：
- 峰值内存是否下降
- 总时间是否在你的机器上受益

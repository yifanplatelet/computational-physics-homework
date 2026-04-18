# 基于你原始代码的 4 个可直接 benchmark 的优化方向

我把最适合直接落到你现有 `Chudnovsky + Binary Splitting` 基线上的 4 个版本分别做成了独立子包：

1. `gcd_reduced_bs`
   - 在 exact binary splitting 的叶子和每次合并后，对 `(Q, T)` 做一次安全的 `gcd` 约分。
   - 这是保守版、易于验证的 GCD reduction。

2. `top_truncated_bs`
   - 下层仍是 exact GMP 子树。
   - 顶层若干层改成 MPFR 上的 `r = P/Q`、`s = T/Q` 近似合并。
   - 用来观察“顶层截断 / 混合精度合并”是否有利于你的机器。

3. `skewed_backward_bs`
   - 改成偏斜切分。
   - 小区间使用 backwards accumulation，从大索引往小索引做块式合并。

4. `newton_endgame_bs`
   - 主干仍是 exact binary splitting。
   - 最终 `Q/T` 用 Newton 倒数迭代替代一次直接大除法。

## 为什么没有把 modular reconstruction 和高阶 Chudnovsky 也一起做成可运行包

这两类都属于“大改”：

- modular / CRT / rational reconstruction
  - 需要一整套模运算、CRT 和重构管线
  - 已经不是在你当前代码上打补丁

- 高阶 Chudnovsky / 153 digits per term 那类扩展
  - 需要换公式、换常量生成、换整套叶子构造逻辑
  - 更像重写一个新项目

如果你愿意，我下一步可以继续补这两类中的其中一个，但我更建议先把这 4 个版本在 Ubuntu 上跑出你的真实 benchmark 数据。

# Newton endgame on top of exact binary splitting

这个版本不改主干的 exact binary splitting，改的是**最后的除法阶段**：

- 先照常算出 exact `Q` 和 `T`
- 再用 Newton 迭代计算 `1 / T`
- 最终 `pi = 426880 * sqrt(10005) * Q * (1 / T)`

它的目标是测试：
- 在你的机器上，最终大除法是否值得换成“倒数 + 乘法”的路线

## Ubuntu 编译

```bash
sudo apt-get update
sudo apt-get install -y build-essential libgmp-dev libmpfr-dev
gcc -O3 -march=native -std=c11 pi_bs_newton_endgame.c -lmpfr -lgmp -o pi_bs_newton_endgame
```

## 运行

```bash
./pi_bs_newton_endgame 1000000
```

参数：
- 第一个参数：`digits`

## 说明

这个版本改的是 `final_eval_time`，不是 `binary_splitting_time`。
如果 profile 显示你的瓶颈几乎全在 `mpz_mul` 和递归树上，那这个版本通常不会有很大优势。

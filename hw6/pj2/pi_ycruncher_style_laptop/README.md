# pi_ycruncher_style_laptop

这不是 y-cruncher 本体，也不是 y-cruncher 的源码复刻。

它是一个 **y-cruncher 风格** 的实现：
- Chudnovsky 公式
- exact binary splitting
- OpenMP task 并行
- 叶子块反向累加
- 偏向笔记本/32G 内存环境的保守参数

## Ubuntu 编译

```bash
sudo apt-get update
sudo apt-get install -y build-essential libgmp-dev libmpfr-dev
make
```

## 运行

推荐先这样：

```bash
./run_32g.sh 10000000 5 128 4096 1
```

或者直接：

```bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=close
export OMP_PLACES=cores
./pi_ycruncher_style 10000000 5 128 4096 1
```

参数：
1. digits
2. task_depth
3. leaf_size
4. task_min_size
5. print_pi (1 打印前 50 位，0 不打印)

## 建议起点

- 1e7 位： `5 128 4096`
- 1e8 位： `6 256 4096`

## 说明

- 如果你的当前最优参数已经更快，这份代码不一定会超越它。
- 它的目的，是给你一份 **可在 32G 笔记本上直接跑的、接近 y-cruncher 主线思路** 的干净版本。

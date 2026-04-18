# pi 高精度优化方案

这个压缩包里给了两套代码方案：

1. `bs_single/pi_bs_single.c`
   - 单线程
   - Chudnovsky + Binary Splitting
   - 适合先验证正确性和速度基线

2. `bs_openmp/pi_bs_openmp.c`
   - OpenMP 并行版本
   - 在二分拆分递归树的前几层并行
   - 适合多核 CPU

## 为什么这两版比“逐项 MPFR 累加”快

原先慢的根本原因是：
- 每一项都做高精度浮点运算
- 每一项都做一次全精度加法
- 每一项都重新参与高精度除法/归一化

这两版改成了：
- 用 **Binary Splitting** 组织 Chudnovsky 求和
- 中间过程尽量放在 **GMP 大整数** 上
- MPFR 只在最后做 `sqrt(10005)` 和一次高精度除法

## 公式

Chudnovsky 公式：

1/pi = 12 / (640320^(3/2)) * sum_k [ (-1)^k (6k)! (13591409 + 545140134k) / ((3k)! (k!)^3 640320^(3k)) ]

Binary Splitting 基本区间：

- P(a,b)
- Q(a,b)
- T(a,b)

叶子：
- P = (6a-5)(2a-1)(6a-1)
- Q = a^3 * 10939058860032000
- T = P * (13591409 + 545140134a)，奇数项取负

合并：
- P = P1 * P2
- Q = Q1 * Q2
- T = T1 * Q2 + P1 * T2

最后：
- pi = (426880 * sqrt(10005) * Q) / T

## 编译示例

### Ubuntu / Debian

先安装依赖：

```bash
sudo apt-get install build-essential libgmp-dev libmpfr-dev
```

### 单线程版

```bash
gcc -O3 -march=native -std=c11 bs_single/pi_bs_single.c -lmpfr -lgmp -o pi_bs_single
```

运行：

```bash
./pi_bs_single 1000000
```

### OpenMP 并行版

```bash
gcc -O3 -march=native -fopenmp -std=c11 bs_openmp/pi_bs_openmp.c -lmpfr -lgmp -o pi_bs_openmp
```

运行：

```bash
./pi_bs_openmp 1000000 4
```

第二个参数 `4` 是并行递归深度，通常取 `3~6` 之间试机器。

## 调参建议

### 位数

第一个命令行参数是目标小数位数，例如：

```bash
./pi_bs_single 100000
./pi_bs_single 1000000
./pi_bs_single 10000000
```

### OpenMP 线程数

```bash
export OMP_NUM_THREADS=8
./pi_bs_openmp 1000000 4
```

### OpenMP 深度

- 深度太小：并行不足
- 深度太大：任务调度开销会明显上升

一般建议：
- 4 核：`3~4`
- 8 核：`4~5`
- 16 核：`5~6`

## 关于“1 分钟内上亿位”

这两份代码已经比“逐项 MPFR 累加”快很多，但如果你的真实目标是：

- `1e8` 位
- `<= 60 s`

那通常还需要：

- 更强的硬件
- 更激进的并行
- 更好的内存管理
- GMP/MPFR 底层对 FFT 乘法的良好支持
- 甚至直接使用专门软件（例如 y-cruncher）

也就是说，这两份代码是**从错误方向切换到正确方向**，但不保证在普通机器上就能直接达到“1 分钟上亿位”。

## 进一步优化方向

如果你还要继续往上冲：

1. 只在高层并行，低层串行，减少小任务
2. 给递归做对象池，减少 `mpz_init/mpz_clear`
3. 用更细的阈值控制叶子区间大小，而不是固定 `32`
4. 把最终输出改成只输出哈希/差值，避免 I/O 干扰
5. 做实际 profiling，再决定瓶颈是乘法、内存还是任务调度

# Skewed splitting + backwards leaf accumulation

这是一个把两个工程优化一起放进去的版本：

1. **skewed splitting**
   - 不再固定 `(a + b) / 2` 中点切分
   - 默认改为 `2/3` 偏斜切分

2. **backwards leaf accumulation**
   - 小区间不再继续递归到单项
   - 直接从大索引往小索引做块式合并

## Ubuntu 编译

```bash
sudo apt-get update
sudo apt-get install -y build-essential libgmp-dev libmpfr-dev
gcc -O3 -march=native -std=c11 pi_bs_skewed_backward.c -lmpfr -lgmp -o pi_bs_skewed_backward
```

## 运行

```bash
./pi_bs_skewed_backward 1000000 32 2 3
```

参数：
- `digits`
- `leaf_size`
- `skew_num`
- `skew_den`

默认 `2/3` 表示偏斜切分点约在 `a + (b-a)*2/3`。

## 建议测试

```bash
./pi_bs_skewed_backward 1000000 16 1 2
./pi_bs_skewed_backward 1000000 32 2 3
./pi_bs_skewed_backward 1000000 64 3 4
```

重点看：
- `leaf_size` 改变时的时间
- 偏斜切分是否比中点切分更适合你的 CPU / 内存层级

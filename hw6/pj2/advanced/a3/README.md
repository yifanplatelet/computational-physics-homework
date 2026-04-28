# advanced/a3

这一版只做三件事：

- 根节点和右脊只算 `Q/T`，不再保留无用的根部 `P`
- 用 `mpz_realloc2()` 给叶子和 merge 结果提前扩容，减少 GMP 反复扩容
- 把 merge 阶段的 task 并行限制在顶层，避免 OpenMP 任务爆炸

算法仍然是：

- Chudnovsky
- exact binary splitting
- GMP 大整数
- MPFR 末端求值

## 编译

```bash
make
```

## 运行

建议配合线程绑定：

```bash
export OMP_PROC_BIND=close
export OMP_PLACES=cores
./pi_bs_qt_tuned 100000000 0 256 16384 1 262144 1
```

参数：

1. `digits`
2. `task_depth`，`0` 表示自动
3. `leaf_size`
4. `task_min_size`
5. `merge_task_depth`
6. `merge_task_min_size`
7. `print_pi`

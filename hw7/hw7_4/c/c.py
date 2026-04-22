import matplotlib
matplotlib.use("Agg")   # Ubuntu 无图形界面时可直接保存

import numpy as np
import matplotlib.pyplot as plt

# 问题设置
x = 0.5
true_val = np.cos(x)
eps_tilde = 7e-17

# 200/201个点
i = np.arange(0, 201)
h = 10.0 ** (-0.1 * i)

# 数值导数（前向差分）
approx = (np.sin(x + h) - np.sin(x)) / h
err = np.abs(approx - true_val)

# 理论误差模型（取绝对值以便画 log10）
theory = np.abs(0.5 * h * np.sin(x) + eps_tilde / h)

# 数值最优 h
best_index = np.argmin(err)
best_h = h[best_index]
best_err = err[best_index]

# 作图
plt.figure(figsize=(8, 5.5))

# 数据点：蓝色
plt.plot(
    np.log10(h),
    np.log10(err),
    'o',
    ms=3,
    color='blue',
    label='numerical error'
)

# 理论线：橙色
plt.plot(
    np.log10(h),
    np.log10(theory),
    '-',
    linewidth=2,
    color='orange',
    label=r'theory: $| \frac{h}{2}\sin(0.5)+\tilde{\varepsilon}/h |$'
)

# 最优位置：红色
plt.axvline(
    np.log10(best_h),
    linestyle='--',
    linewidth=1.5,
    color='red'
)

plt.scatter(
    np.log10(best_h),
    np.log10(best_err),
    color='red',
    s=40,
    zorder=5
)

plt.text(
    np.log10(best_h) + 0.3,
    np.log10(best_err) + 0.3,
    f'best h ≈ {best_h:.1e}',
    color='red',
    fontsize=11
)

# 坐标与标题
plt.xlabel(r'$\log_{10}(h)$')
plt.ylabel(r'$\log_{10}(|\mathrm{error}|)$')
plt.title(r'Forward-difference derivative error for $f(x)=\sin x$ at $x=0.5$')

# 不加网格
# plt.grid(False)

plt.legend()
plt.tight_layout()

# 保存图像
plt.savefig("derivative_error_plot_with_theory.png", dpi=300, bbox_inches="tight")
print(f"best h = {best_h:.6e}")
print("图像已保存为 derivative_error_plot_with_theory.png")
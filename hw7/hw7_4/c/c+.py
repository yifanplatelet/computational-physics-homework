import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt

# =========================================================
# 基本设置
# =========================================================
x0 = 0.5
true_val = np.cos(x0)
eps = np.finfo(float).eps

i = np.arange(0, 201)
h = 10.0 ** (-0.1 * i)
logh = np.log10(h)

def f(x):
    return np.sin(x)

# =========================================================
# 三种差分格式
# =========================================================
# O(h)
D1 = (f(x0 + h) - f(x0)) / h

# O(h^2)
D2 = (f(x0 + h) - f(x0 - h)) / (2.0 * h)

# O(h^4)
D4 = (-f(x0 + 2.0*h) + 8.0*f(x0 + h) - 8.0*f(x0 - h) + f(x0 - 2.0*h)) / (12.0 * h)

err1 = np.abs(D1 - true_val)
err2 = np.abs(D2 - true_val)
err4 = np.abs(D4 - true_val)

# =========================================================
# 理论误差模型
# =========================================================
A1 = abs(np.sin(x0)) / 2.0
A2 = abs(np.cos(x0)) / 6.0
A4 = abs(np.cos(x0)) / 30.0

B1 = 2.0 * abs(np.sin(x0))
B2 = 2.0 * abs(np.sin(x0))
B4 = (18.0 / 12.0) * abs(np.sin(x0))

theory1 = A1 * h + B1 * eps / h
theory2 = A2 * h**2 + B2 * eps / h
theory4 = A4 * h**4 + B4 * eps / h

# =========================================================
# 最优点
# =========================================================
idx_num1 = np.argmin(err1)
idx_num2 = np.argmin(err2)
idx_num4 = np.argmin(err4)

best_h_num1, best_e_num1 = h[idx_num1], err1[idx_num1]
best_h_num2, best_e_num2 = h[idx_num2], err2[idx_num2]
best_h_num4, best_e_num4 = h[idx_num4], err4[idx_num4]

def theoretical_best_h(A, B, p, eps):
    return (B * eps / (p * A)) ** (1.0 / (p + 1.0))

def theoretical_error(A, B, p, eps, hbest):
    return A * hbest**p + B * eps / hbest

best_h_th1 = theoretical_best_h(A1, B1, 1, eps)
best_h_th2 = theoretical_best_h(A2, B2, 2, eps)
best_h_th4 = theoretical_best_h(A4, B4, 4, eps)

best_e_th1 = theoretical_error(A1, B1, 1, eps, best_h_th1)
best_e_th2 = theoretical_error(A2, B2, 2, eps, best_h_th2)
best_e_th4 = theoretical_error(A4, B4, 4, eps, best_h_th4)

# =========================================================
# 配色
# 数值点：蓝色系
# 拟合线：橙色系
# 理论最优：红色系
# 实际最优：紫色系
# =========================================================
point_colors  = ['#94adff', '#6c88e1', '#084594']   # 蓝
line_colors   = ['#ffbe49', '#cd9318', '#6f4300']   # 橙
red_colors    = ['#ffbe49', '#cd9318', '#6f4300']   # 红
purple_colors = ['#94adff', '#6c88e1', '#084594']   # 紫

labels = [r'$O(h)$', r'$O(h^2)$', r'$O(h^4)$']
errs = [err1, err2, err4]
theories = [theory1, theory2, theory4]

best_h_num = [best_h_num1, best_h_num2, best_h_num4]
best_e_num = [best_e_num1, best_e_num2, best_e_num4]

best_h_th = [best_h_th1, best_h_th2, best_h_th4]
best_e_th = [best_e_th1, best_e_th2, best_e_th4]

# =========================================================
# 作图
# =========================================================
fig, ax = plt.subplots(figsize=(12, 8))

for k in range(3):
    # 数值点：蓝色系
    ax.plot(
        logh,
        np.log10(errs[k]),
        'o',
        ms=3.2,
        color=point_colors[k],
        label=rf'numerical error {labels[k]}'
    )

    # 拟合线/理论线：橙色系
    ax.plot(
        logh,
        np.log10(theories[k]),
        '-',
        linewidth=2.4,
        color=line_colors[k],
        label=rf'theory fit {labels[k]}'
    )

    # 理论最优：红色虚线 + 红点
    ax.axvline(
        np.log10(best_h_th[k]),
        linestyle='--',
        linewidth=1.8,
        color=red_colors[k],
        alpha=0.95
    )
    ax.scatter(
        np.log10(best_h_th[k]),
        np.log10(best_e_th[k]),
        s=90,
        color=red_colors[k],
        zorder=6
    )

    # 实际最优：紫色虚线 + 紫点
    ax.axvline(
        np.log10(best_h_num[k]),
        linestyle='--',
        linewidth=1.8,
        color=purple_colors[k],
        alpha=0.95
    )
    ax.scatter(
        np.log10(best_h_num[k]),
        np.log10(best_e_num[k]),
        s=90,
        color=purple_colors[k],
        zorder=6
    )

text_offsets_red = [(0.15, 0.45), (0.15, 0.15), (0.15, -0.15)]
text_offsets_pur = [(-3.65, -0.60), (-3.15, -0.90), (-2.65, -1.20)]

for k in range(3):
    ax.text(
        np.log10(best_h_th[k]) + text_offsets_red[k][0],
        np.log10(best_e_th[k]) + text_offsets_red[k][1],
        rf'theory {labels[k]}:' + '\n' +
        rf'$h_{{best}}\approx {best_h_th[k]:.2e}$' + '\n' +
        rf'$E_{{min}}\approx {best_e_th[k]:.2e}$',
        color=red_colors[k],
        fontsize=10
    )

    ax.text(
        np.log10(best_h_num[k]) + text_offsets_pur[k][0],
        np.log10(best_e_num[k]) + text_offsets_pur[k][1],
        rf'numerical {labels[k]}:' + '\n' +
        rf'$h_{{best}}\approx {best_h_num[k]:.2e}$' + '\n' +
        rf'$E_{{min}}\approx {best_e_num[k]:.2e}$',
        color=purple_colors[k],
        fontsize=10
    )

ax.set_xlabel(r'$\log_{10}(h)$', fontsize=18)
ax.set_ylabel(r'$\log_{10}(|\mathrm{error}|)$', fontsize=18)
ax.set_title(r'Numerical differentiation error for $f(x)=\sin x$ at $x=0.5$', fontsize=22)

ax.grid(False)

for spine in ax.spines.values():
    spine.set_linewidth(1.2)

ax.tick_params(axis='both', labelsize=14, width=1.2, length=7)
ax.legend(fontsize=12, loc='upper right', frameon=True)

plt.tight_layout()
plt.savefig("multi_order_derivative_error_comparison.png", dpi=300, bbox_inches="tight")
print("图像已保存为 multi_order_derivative_error_comparison.png")
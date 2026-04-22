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
# 四种差分格式
# =========================================================
# O(h): forward difference
D1 = (f(x0 + h) - f(x0)) / h

# O(h^2): central difference
D2 = (f(x0 + h) - f(x0 - h)) / (2.0 * h)

# O(h^4): 5-point central difference
D4 = (-f(x0 + 2.0*h) + 8.0*f(x0 + h) - 8.0*f(x0 - h) + f(x0 - 2.0*h)) / (12.0 * h)

# O(h^6): 7-point central difference
D6 = (
    -f(x0 - 3.0*h)
    + 9.0*f(x0 - 2.0*h)
    - 45.0*f(x0 - h)
    + 45.0*f(x0 + h)
    - 9.0*f(x0 + 2.0*h)
    + f(x0 + 3.0*h)
) / (60.0 * h)

# 数值误差
err1 = np.abs(D1 - true_val)
err2 = np.abs(D2 - true_val)
err4 = np.abs(D4 - true_val)
err6 = np.abs(D6 - true_val)

# =========================================================
# 理论误差模型
# E(h) ≈ A*h^p + B*eps/h
# =========================================================
# 截断误差主项系数
A1 = abs(np.sin(x0)) / 2.0          # O(h)
A2 = abs(np.cos(x0)) / 6.0          # O(h^2)
A4 = abs(np.cos(x0)) / 30.0         # O(h^4)
A6 = abs(np.cos(x0)) / 140.0        # O(h^6)

# roundoff 误差量级系数
# 用“格式中系数绝对值和 / 分母”乘函数量级做估计
B1 = 2.0 * abs(np.sin(x0))                  # (1+1)/1
B2 = 2.0 * abs(np.sin(x0))                  # (1+1)/2 * 2 -> 量级上近似仍取 2|sin|
B4 = (1.0 + 8.0 + 8.0 + 1.0) / 12.0 * abs(np.sin(x0))
B6 = (1.0 + 9.0 + 45.0 + 45.0 + 9.0 + 1.0) / 60.0 * abs(np.sin(x0))

theory1 = A1 * h + B1 * eps / h
theory2 = A2 * h**2 + B2 * eps / h
theory4 = A4 * h**4 + B4 * eps / h
theory6 = A6 * h**6 + B6 * eps / h

# =========================================================
# 实际最优点
# =========================================================
idx_num1 = np.argmin(err1)
idx_num2 = np.argmin(err2)
idx_num4 = np.argmin(err4)
idx_num6 = np.argmin(err6)

best_h_num1, best_e_num1 = h[idx_num1], err1[idx_num1]
best_h_num2, best_e_num2 = h[idx_num2], err2[idx_num2]
best_h_num4, best_e_num4 = h[idx_num4], err4[idx_num4]
best_h_num6, best_e_num6 = h[idx_num6], err6[idx_num6]

# =========================================================
# 理论最优点
# 如果 E(h)=A*h^p + B*eps/h
# h_best = (B*eps/(p*A))^(1/(p+1))
# =========================================================
def theoretical_best_h(A, B, p, eps):
    return (B * eps / (p * A)) ** (1.0 / (p + 1.0))

def theoretical_error(A, B, p, eps, hbest):
    return A * hbest**p + B * eps / hbest

best_h_th1 = theoretical_best_h(A1, B1, 1, eps)
best_h_th2 = theoretical_best_h(A2, B2, 2, eps)
best_h_th4 = theoretical_best_h(A4, B4, 4, eps)
best_h_th6 = theoretical_best_h(A6, B6, 6, eps)

best_e_th1 = theoretical_error(A1, B1, 1, eps, best_h_th1)
best_e_th2 = theoretical_error(A2, B2, 2, eps, best_h_th2)
best_e_th4 = theoretical_error(A4, B4, 4, eps, best_h_th4)
best_e_th6 = theoretical_error(A6, B6, 6, eps, best_h_th6)

# =========================================================
# 配色
# =========================================================
# numerical points: 蓝色系
point_colors  = ['#c6dbef', '#6baed6', '#2171b5', '#08306b']

# theory fit lines: 橙色系
line_colors   = ['#fcbba1', '#fc9272', '#ef3b2c', '#99000d']

# theory optimum: 红色系
red_colors    = ['#fcbba1', '#fc9272', '#ef3b2c', '#99000d']

# numerical optimum: 紫色系
purple_colors = ['#c6dbef', '#6baed6', '#2171b5', '#08306b']

labels = [r'$O(h)$', r'$O(h^2)$', r'$O(h^4)$', r'$O(h^6)$']
errs = [err1, err2, err4, err6]
theories = [theory1, theory2, theory4, theory6]

best_h_num = [best_h_num1, best_h_num2, best_h_num4, best_h_num6]
best_e_num = [best_e_num1, best_e_num2, best_e_num4, best_e_num6]

best_h_th = [best_h_th1, best_h_th2, best_h_th4, best_h_th6]
best_e_th = [best_e_th1, best_e_th2, best_e_th4, best_e_th6]

# =========================================================
# 作图
# =========================================================
fig, ax = plt.subplots(figsize=(14, 9))

for k in range(4):
    # numerical points
    ax.plot(
        logh,
        np.log10(errs[k]),
        'o',
        ms=3.0,
        color=point_colors[k],
        label=rf'numerical error {labels[k]}'
    )

    # theory fit lines
    ax.plot(
        logh,
        np.log10(theories[k]),
        '-',
        linewidth=2.3,
        color=line_colors[k],
        label=rf'theory fit {labels[k]}'
    )

    # theory optimum: red dashed line + red point
    ax.axvline(
        np.log10(best_h_th[k]),
        linestyle='--',
        linewidth=1.6,
        color=red_colors[k],
        alpha=0.95
    )
    ax.scatter(
        np.log10(best_h_th[k]),
        np.log10(best_e_th[k]),
        s=80,
        color=red_colors[k],
        zorder=6
    )

    # numerical optimum: purple dashed line + purple point
    ax.axvline(
        np.log10(best_h_num[k]),
        linestyle='--',
        linewidth=1.6,
        color=purple_colors[k],
        alpha=0.95
    )
    ax.scatter(
        np.log10(best_h_num[k]),
        np.log10(best_e_num[k]),
        s=80,
        color=purple_colors[k],
        zorder=6
    )

# =========================================================
# 标注文字
# theory 放右边
# numerical 放左边
# =========================================================
text_offsets_red = [
    (-0.5,  0.70),
    (-0.5,  0.30),
    (1.18,  0.30),
    (0.28, -0.50),
]

text_offsets_pur = [
    (-0.48, -0.60),
    (-0.48, -0.60),
    (-0.48, -0.60),
    (2.28, -0.60),
]

for k in range(4):
    # theory text
    ax.text(
        np.log10(best_h_th[k]) + text_offsets_red[k][0],
        np.log10(best_e_th[k]) + text_offsets_red[k][1],
        rf'theory {labels[k]}:' + '\n' +
        rf'$h_{{best}}\approx {best_h_th[k]:.2e}$' + '\n' +
        rf'$E_{{min}}\approx {best_e_th[k]:.2e}$',
        color=red_colors[k],
        fontsize=10,
        ha='left'
    )

    # numerical text: 左边
    ax.text(
        np.log10(best_h_num[k]) + text_offsets_pur[k][0],
        np.log10(best_e_num[k]) + text_offsets_pur[k][1],
        rf'numerical {labels[k]}:' + '\n' +
        rf'$h_{{best}}\approx {best_h_num[k]:.2e}$' + '\n' +
        rf'$E_{{min}}\approx {best_e_num[k]:.2e}$',
        color=purple_colors[k],
        fontsize=10,
        ha='right'
    )

# =========================================================
# 坐标轴与样式
# =========================================================
ax.set_xlabel(r'$\log_{10}(h)$', fontsize=18)
ax.set_ylabel(r'$\log_{10}(|\mathrm{error}|)$', fontsize=18)
ax.set_title(r'Numerical differentiation error for $f(x)=\sin x$ at $x=0.5$', fontsize=22)

ax.grid(False)

for spine in ax.spines.values():
    spine.set_linewidth(1.2)

ax.tick_params(axis='both', labelsize=14, width=1.2, length=7)

ax.legend(fontsize=11, loc='upper right', frameon=True, ncol=2)

plt.tight_layout()
plt.savefig("multi_order_derivative_error_comparison_with_Oh6.png", dpi=300, bbox_inches="tight")
print("图像已保存为 multi_order_derivative_error_comparison_with_Oh6.png")

# =========================================================
# 输出最优结果
# =========================================================
print("\n===== Theoretical optima =====")
print(f"O(h):   h_best = {best_h_th1:.6e}, E_min = {best_e_th1:.6e}")
print(f"O(h^2): h_best = {best_h_th2:.6e}, E_min = {best_e_th2:.6e}")
print(f"O(h^4): h_best = {best_h_th4:.6e}, E_min = {best_e_th4:.6e}")
print(f"O(h^6): h_best = {best_h_th6:.6e}, E_min = {best_e_th6:.6e}")

print("\n===== Numerical optima =====")
print(f"O(h):   h_best = {best_h_num1:.6e}, E_min = {best_e_num1:.6e}")
print(f"O(h^2): h_best = {best_h_num2:.6e}, E_min = {best_e_num2:.6e}")
print(f"O(h^4): h_best = {best_h_num4:.6e}, E_min = {best_e_num4:.6e}")
print(f"O(h^6): h_best = {best_h_num6:.6e}, E_min = {best_e_num6:.6e}")
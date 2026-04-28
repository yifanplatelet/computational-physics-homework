import numpy as np
import matplotlib.pyplot as plt
import math

# 原函数
def f(x):
    return x * np.exp(x)

# 精确一阶导数
def fp_exact(x):
    return (x + 1) * np.exp(x)

x0 = 2.0
exact = fp_exact(x0)

print("Exact f'(2) =", exact)

# =========================
# Part 1: use table data h=0.1
# =========================

h = 0.1

table = {
    1.8: 10.889365,
    1.9: 12.703199,
    2.0: 14.778112,
    2.1: 17.148957,
    2.2: 19.855030
}

# 一阶差分公式
forward = (table[2.1] - table[2.0]) / h
backward = (table[2.0] - table[1.9]) / h
central = (table[2.1] - table[1.9]) / (2*h)

# 三点前向公式
forward_3 = (
    -3*table[2.0]
    + 4*table[2.1]
    - table[2.2]
) / (2*h)

# 三点后向公式
backward_3 = (
    3*table[2.0]
    - 4*table[1.9]
    + table[1.8]
) / (2*h)

# 五点中心差分公式
central_5 = (
    table[1.8]
    - 8*table[1.9]
    + 8*table[2.1]
    - table[2.2]
) / (12*h)

methods_table = {
    "Forward difference": forward,
    "Backward difference": backward,
    "Central difference": central,
    "Forward 3-point": forward_3,
    "Backward 3-point": backward_3,
    "Central 5-point": central_5
}

print("\nUsing table data, h = 0.1")
print("{:<25s} {:>18s} {:>18s}".format(
    "Method", "Approximation", "Abs Error"
))

for name, value in methods_table.items():
    print("{:<25s} {:>18.8f} {:>18.8f}".format(
        name, value, abs(value - exact)
    ))

# =========================
# Part 2: different h values
# =========================

print("\nErrors for different h")
print("{:<10s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s}".format(
    "h",
    "Forward",
    "Backward",
    "Central",
    "Forward3",
    "Backward3",
    "Central5"
))

hs = [10**(-k) for k in range(1, 13)]

for h in hs:
    forward_h = (f(x0+h) - f(x0)) / h
    backward_h = (f(x0) - f(x0-h)) / h
    central_h = (f(x0+h) - f(x0-h)) / (2*h)

    forward_3_h = (
        -3*f(x0)
        + 4*f(x0+h)
        - f(x0+2*h)
    ) / (2*h)

    backward_3_h = (
        3*f(x0)
        - 4*f(x0-h)
        + f(x0-2*h)
    ) / (2*h)

    central_5_h = (
        f(x0-2*h)
        - 8*f(x0-h)
        + 8*f(x0+h)
        - f(x0+2*h)
    ) / (12*h)

    print("{:<10.1e} {:>15.6e} {:>15.6e} {:>15.6e} {:>15.6e} {:>15.6e} {:>15.6e}".format(
        h,
        abs(forward_h - exact),
        abs(backward_h - exact),
        abs(central_h - exact),
        abs(forward_3_h - exact),
        abs(backward_3_h - exact),
        abs(central_5_h - exact)
    ))

# =========================
# Part 3: plot comparison
# =========================

x = np.linspace(1.8, 2.2, 400)

f0 = table[2.0]

def tangent_line(slope):
    return f0 + slope * (x - x0)

plt.figure(figsize=(10, 6))

# 原函数
plt.plot(x, f(x), label="Original function")

# 各差分公式得到的切线近似
for name, slope in methods_table.items():
    plt.plot(x, tangent_line(slope), linestyle="--", label=name)

# 表格给出的数据点
plt.scatter(
    list(table.keys()),
    list(table.values()),
    label="Given data"
)

plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Comparison of original function and tangent approximations")
plt.legend()
plt.grid(False)
plt.savefig("comparison_plot.png")
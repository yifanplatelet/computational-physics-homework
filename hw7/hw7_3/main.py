import numpy as np
import matplotlib.pyplot as plt
import math

def f(x):
    return x * np.exp(x**2)

def fp_exact(x):
    return (1 + 2*x**2) * np.exp(x**2)

def fpp_exact(x):
    return (6*x + 4*x**3) * np.exp(x**2)

x0 = 2.0
exact = fpp_exact(x0)

print("Exact f''(2) =", exact)

# =========================
# Part 1: use table data h=0.1
# =========================

h = 0.1

table = {
    1.8: 45.960699,
    1.9: 70.235500,
    2.0: 109.196300,
    2.1: 172.765873,
    2.2: 278.232574
}

forward_3 = (table[2.0] - 2*table[2.1] + table[2.2]) / h**2
backward_3 = (table[2.0] - 2*table[1.9] + table[1.8]) / h**2
central_3 = (table[2.1] - 2*table[2.0] + table[1.9]) / h**2

central_5 = (
    -table[1.8]
    + 16*table[1.9]
    - 30*table[2.0]
    + 16*table[2.1]
    - table[2.2]
) / (12*h**2)

methods_table = {
    "Forward 3-point": forward_3,
    "Backward 3-point": backward_3,
    "Central 3-point": central_3,
    "Central 5-point": central_5
}

print("\nUsing table data, h = 0.1")
print("{:<20s} {:>18s} {:>18s}".format("Method", "Approximation", "Abs Error"))

for name, value in methods_table.items():
    print("{:<20s} {:>18.6f} {:>18.6f}".format(
        name, value, abs(value - exact)
    ))

# =========================
# Part 2: different h values
# =========================

print("\nErrors for different h")
print("{:<10s} {:>15s} {:>15s} {:>15s} {:>15s}".format(
    "h", "Forward3", "Backward3", "Central3", "Central5"
))

hs = [10**(-k) for k in range(1, 13)]

for h in hs:
    forward_3_h = (f(x0) - 2*f(x0+h) + f(x0+2*h)) / h**2
    backward_3_h = (f(x0) - 2*f(x0-h) + f(x0-2*h)) / h**2
    central_3_h = (f(x0+h) - 2*f(x0) + f(x0-h)) / h**2

    central_5_h = (
        -f(x0-2*h)
        + 16*f(x0-h)
        - 30*f(x0)
        + 16*f(x0+h)
        - f(x0+2*h)
    ) / (12*h**2)

    print("{:<10.1e} {:>15.6e} {:>15.6e} {:>15.6e} {:>15.6e}".format(
        h,
        abs(forward_3_h - exact),
        abs(backward_3_h - exact),
        abs(central_3_h - exact),
        abs(central_5_h - exact)
    ))

# =========================
# Part 3: plot comparison
# =========================

x = np.linspace(1.8, 2.2, 400)

f0 = f(x0)
fp0 = fp_exact(x0)

def taylor_curve(fpp_approx):
    return f0 + fp0*(x - x0) + 0.5*fpp_approx*(x - x0)**2

plt.figure(figsize=(10, 6))

plt.plot(x, f(x), label="Original function")

for name, fpp_approx in methods_table.items():
    plt.plot(x, taylor_curve(fpp_approx), linestyle="--", label=name)

plt.scatter(list(table.keys()), list(table.values()), label="Given data")

plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Comparison of original function and second-order approximations")
plt.legend()
plt.grid(False)
plt.savefig("comparison_plot.png")
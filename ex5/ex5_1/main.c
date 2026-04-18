#include <stdio.h>
#include <math.h>

#define EPS 1e-10       // 精度要求
#define MAX_ITER 1000   // 最大迭代次数

// 定义函数 f(x) = 2 - x - e^(-x)
double f(double x) {
    return 2.0 - x - exp(-x);
}

// 二分法函数
double bisection(double a, double b, double eps, int max_iter) {
    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0) {
        printf("区间 [%.10f, %.10f] 内没有满足二分法条件（两端异号）的根。\n", a, b);
        return NAN;
    }

    double mid, fmid;
    int iter = 0;

    while ((b - a) / 2.0 > eps && iter < max_iter) {
        mid = (a + b) / 2.0;
        fmid = f(mid);

        // 如果刚好找到根
        if (fabs(fmid) < eps) {
            return mid;
        }

        // 根据符号变化缩小区间
        if (fa * fmid < 0) {
            b = mid;
            fb = fmid;
        } else {
            a = mid;
            fa = fmid;
        }

        iter++;
    }

    return (a + b) / 2.0;
}

int main() {
    double root1, root2;

    // 左侧根：在 [-2, -1]
    root1 = bisection(-2.0, -1.0, EPS, MAX_ITER);

    // 右侧根：在 [1, 2]
    root2 = bisection(1.0, 2.0, EPS, MAX_ITER);

    printf("方程 f(x) = 2 - x - e^(-x) 的两个实根为：\n");
    printf("root1 = %.12f\n", root1);
    printf("root2 = %.12f\n", root2);

    printf("\n验证：\n");
    printf("f(root1) = %.12e\n", f(root1));
    printf("f(root2) = %.12e\n", f(root2));

    return 0;
}
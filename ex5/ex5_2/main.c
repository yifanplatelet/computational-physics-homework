#include <stdio.h>
#include <math.h>

#define EPS 1e-10
#define MAX_ITER 100

// 定义函数 f(x) = 2 - x - e^(-x)
double f(double x) {
    return 2.0 - x - exp(-x);
}

// 割线法
double secant(double x0, double x1) {
    double x2;
    int iter = 0;

    while (fabs(x1 - x0) > EPS && iter < MAX_ITER) {
        double f0 = f(x0);
        double f1 = f(x1);

        if (fabs(f1 - f0) < 1e-15) {
            printf("分母过小，迭代失败！\n");
            return x1;
        }

        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        printf("迭代 %2d: x = %.12f, f(x) = %.12e\n", iter + 1, x2, f(x2));

        x0 = x1;
        x1 = x2;
        iter++;
    }

    return x1;
}

int main() {
    double root1, root2;

    printf("求第一个根（初值取 x0=-2, x1=-1）:\n");
    root1 = secant(-2.0, -1.0);
    printf("第一个根约为: %.12f\n\n", root1);

    printf("求第二个根（初值取 x0=1, x1=2）:\n");
    root2 = secant(1.0, 2.0);
    printf("第二个根约为: %.12f\n", root2);

    return 0;
}
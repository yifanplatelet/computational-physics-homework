#include <stdio.h>
#include <math.h>

double f(double x) {
    return x * (x - 1.0);
}

double derivative(double (*func)(double), double x, double delta) {
    return (func(x + delta) - func(x)) / delta;
}

int main() {
    double x = 1.0;
    double deltas[] = {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16, 1e-18};
    int n = sizeof(deltas) / sizeof(deltas[0]);
    double true_value = 2.0 * x - 1.0;

    for (int i = 0; i < n; i++) {
        double d = deltas[i];
        double num = derivative(f, x, d);
        double err = fabs(num - true_value);
        printf("delta = %e, numerical = %.16f, error = %e\n", d, num, err);
    }

    return 0;
}
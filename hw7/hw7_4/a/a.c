#include <stdio.h>

double f(double x) {
    return x * (x - 1.0);
}

double derivative(double (*func)(double), double x, double delta) {
    return (func(x + delta) - func(x)) / delta;
}

int main() {
    double x = 1.0;
    double delta = 1e-2;

    double numerical = derivative(f, x, delta);
    double analytical = 2.0 * x - 1.0;

    printf("numerical derivative = %.10f\n", numerical);
    printf("analytical derivative = %.10f\n", analytical);
    printf("error = %.10f\n", numerical - analytical);

    return 0;
}
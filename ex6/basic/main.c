#include <stdio.h>

double lagrange(double x[], double y[], int n, double t) {
    double result = 0.0;

    for (int i = 0; i < n; i++) {
        double term = y[i];
        for (int j = 0; j < n; j++) {
            if (j != i) {
                term *= (t - x[j]) / (x[i] - x[j]);
            }
        }
        result += term;
    }

    return result;
}

int main() {
    int n = 8;

    /* 用 0~7 表示 1920~1990 */
    double x[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    double y[8] = {106.46, 123.08, 132.12, 152.27, 180.67, 205.05, 227.23, 249.46};

    double p2000 = lagrange(x, y, n, 8);   // 2000
    double p2010 = lagrange(x, y, n, 9);   // 2010
    double p2020 = lagrange(x, y, n, 10);  // 2020

    printf("Prediction for 2000: %.2f\n", p2000);
    printf("Prediction for 2010: %.2f\n", p2010);
    printf("Prediction for 2020: %.2f\n", p2020);

    if (p2020 < 0) {
        printf("Conclusion: The 2020 prediction is NOT reasonable.\n");
    } else {
        printf("Conclusion: The 2020 prediction may be reasonable.\n");
    }

    return 0;
}
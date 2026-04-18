#include <stdio.h>
#include <math.h>

#define TOTAL_YEARS 16

typedef struct {
    int year;
    double pop;
} DataPoint;

/* 牛顿插值：根据 n 个点构造插值多项式，并在 target_x 处求值 */
double newton_interpolation(double x[], double y[], int n, double target_x) {
    double dd[32][32];   // 足够放本题最多 15 个点
    int i, j;

    for (i = 0; i < n; i++) {
        dd[i][0] = y[i];
    }

    for (j = 1; j < n; j++) {
        for (i = 0; i < n - j; i++) {
            dd[i][j] = (dd[i + 1][j - 1] - dd[i][j - 1]) / (x[i + j] - x[i]);
        }
    }

    double result = dd[0][0];
    double term = 1.0;

    for (j = 1; j < n; j++) {
        term *= (target_x - x[j - 1]);
        result += dd[0][j] * term;
    }

    return result;
}

int main() {
    DataPoint data[TOTAL_YEARS] = {
        {2010, 309321666},
        {2011, 311556874},
        {2012, 313830990},
        {2013, 315993715},
        {2014, 318301008},
        {2015, 320635163},
        {2016, 322941311},
        {2017, 324985539},
        {2018, 326687501},
        {2019, 328239523},
        {2020, 331578104},
        {2021, 332100166},
        {2022, 333996304},
        {2023, 336755052},
        {2024, 340003797},
        {2025, 341784857}
    };

    const double actual_2025 = 341784857.0;

    printf("Actual population in 2025 = %.0f\n\n", actual_2025);
    printf("%-20s %-8s %-18s %-18s %-12s\n",
           "Range", "N", "Predicted_2025", "Abs_Error", "Rel_Error");
    printf("-------------------------------------------------------------------------------\n");

    double best_error = 1e100;
    int best_start = -1;
    double best_pred = 0.0;

    /* 分别做：
       2010-2024 -> 预测2025
       2011-2024 -> 预测2025
       ...
       2021-2024 -> 预测2025
    */
    for (int start = 2010; start <= 2021; start++) {
        int n = 2024 - start + 1;   // 用 start..2024，共 n 个点
        double x[32], y[32];

        for (int i = 0; i < n; i++) {
            x[i] = (double)(start + i);
            y[i] = data[(start - 2010) + i].pop;
        }

        double pred = newton_interpolation(x, y, n, 2025.0);
        double abs_error = fabs(pred - actual_2025);
        double rel_error = abs_error / actual_2025 * 100.0;

        printf("%d-%d -> 2025     %-8d %-18.0f %-18.0f %-11.4f%%\n",
               start, 2024, n, pred, abs_error, rel_error);

        if (abs_error < best_error) {
            best_error = abs_error;
            best_start = start;
            best_pred = pred;
        }
    }

    printf("\nBest result: use %d-2024 data to predict 2025\n", best_start);
    printf("Predicted 2025 = %.0f\n", best_pred);
    printf("Absolute error = %.0f\n", best_error);

    return 0;
}
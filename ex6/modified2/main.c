#include <stdio.h>
#include <math.h>

#define TOTAL_YEARS 16

typedef struct {
    int year;
    double pop;
} DataPoint;

/* 牛顿插值 */
double newton_interpolation(double x[], double y[], int n, double target_x) {
    double dd[20][20];
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

    printf("Prediction of USA population in 2026 using N-point interpolation\n");
    printf("---------------------------------------------------------------\n");
    printf("%-18s %-8s %-18s\n", "Range", "N", "Predicted_2026");
    printf("---------------------------------------------------------------\n");

    for (int start = 2010; start <= 2021; start++) {
        int n = 2025 - start + 1;   // start..2025
        double x[20], y[20];

        for (int i = 0; i < n; i++) {
            x[i] = (double)(start + i);
            y[i] = data[start - 2010 + i].pop;
        }

        double pred = newton_interpolation(x, y, n, 2026.0);

        printf("%d-%d            %-8d %-18.0f\n", start, 2025, n, pred);
    }

    printf("\nAnalysis:\n");
    printf("1. Long-interval high-order interpolation gives severe oscillation.\n");
    printf("2. Results from 2010-2025 to 2016-2025 are clearly unreasonable.\n");
    printf("3. Results from 2019-2025, 2020-2025, 2021-2025 are much more reasonable.\n");
    printf("4. The most reasonable estimate is likely from 2020-2025 or 2021-2025.\n");
    printf("5. If one must be chosen, 2020-2025 is a good compromise between stability and data amount.\n");

    return 0;
}
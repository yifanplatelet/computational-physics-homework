#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOTAL_YEARS 16
#define MAXN 20
#define CURVE_POINTS 321
#define CASES 15   /* 2010-2025 到 2024-2025 */

typedef struct {
    int year;
    double value;
} DataPoint;

typedef struct {
    int start_year;
    int end_year;
    int n;
    double pred2026;
    double score;
} Result;

/* ---------- 牛顿插值部分 ---------- */

/* 构造牛顿插值差商系数 */
void build_newton_coeff(double x[], double y[], int n, double coeff[]) {
    double dd[MAXN][MAXN];
    int i, j;

    for (i = 0; i < n; i++) {
        dd[i][0] = y[i];
    }

    for (j = 1; j < n; j++) {
        for (i = 0; i < n - j; i++) {
            dd[i][j] = (dd[i + 1][j - 1] - dd[i][j - 1]) / (x[i + j] - x[i]);
        }
    }

    for (i = 0; i < n; i++) {
        coeff[i] = dd[0][i];
    }
}

/* 牛顿插值求值 */
double newton_eval(double x[], double coeff[], int n, double t) {
    double result = coeff[0];
    double term = 1.0;
    int i;

    for (i = 1; i < n; i++) {
        term *= (t - x[i - 1]);
        result += coeff[i] * term;
    }

    return result;
}

/* 简单评分：只是为了从很多坏结果里找“相对没那么夸张”的 */
double calc_score(double pred, double y2023, double y2024, double y2025) {
    double recent_trend = y2025 + (y2025 - y2024);
    double smooth_ref = (y2023 + y2024 + y2025) / 3.0;

    double s1 = fabs(pred - recent_trend);
    double s2 = fabs(pred - smooth_ref);

    if (pred < 0) {
        s1 += 1000000.0;
        s2 += 1000000.0;
    }

    return s1 + 0.5 * s2;
}

void find_best_k(Result results[], int m, int top_idx[], int k) {
    int used[CASES] = {0};
    int i, j, best;

    for (i = 0; i < k; i++) {
        best = -1;
        for (j = 0; j < m; j++) {
            if (!used[j]) {
                if (best == -1 || results[j].score < results[best].score) {
                    best = j;
                }
            }
        }
        top_idx[i] = best;
        if (best != -1) used[best] = 1;
    }
}

/* ---------- 二次最小二乘拟合部分 ---------- */
/* 拟合形式：y = a*x^2 + b*x + c，其中 x 用相对年份 0,1,...,15 */

void solve_3x3(double A[3][4], double ans[3]) {
    int i, j, k;
    for (i = 0; i < 3; i++) {
        int pivot = i;
        for (j = i + 1; j < 3; j++) {
            if (fabs(A[j][i]) > fabs(A[pivot][i])) {
                pivot = j;
            }
        }

        if (pivot != i) {
            for (k = 0; k < 4; k++) {
                double tmp = A[i][k];
                A[i][k] = A[pivot][k];
                A[pivot][k] = tmp;
            }
        }

        if (fabs(A[i][i]) < 1e-12) {
            printf("Quadratic fitting failed: singular matrix.\n");
            exit(1);
        }

        for (j = i + 1; j < 3; j++) {
            double factor = A[j][i] / A[i][i];
            for (k = i; k < 4; k++) {
                A[j][k] -= factor * A[i][k];
            }
        }
    }

    for (i = 2; i >= 0; i--) {
        ans[i] = A[i][3];
        for (j = i + 1; j < 3; j++) {
            ans[i] -= A[i][j] * ans[j];
        }
        ans[i] /= A[i][i];
    }
}

void quadratic_fit(double x[], double y[], int n, double coef[3]) {
    double s0 = n;
    double s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0;
    double t0 = 0.0, t1 = 0.0, t2 = 0.0;
    int i;

    for (i = 0; i < n; i++) {
        double xi = x[i];
        double yi = y[i];
        s1 += xi;
        s2 += xi * xi;
        s3 += xi * xi * xi;
        s4 += xi * xi * xi * xi;

        t0 += yi;
        t1 += xi * yi;
        t2 += xi * xi * yi;
    }

    /* 正规方程 */
    /* [ s4 s3 s2 ] [a] = [t2] */
    /* [ s3 s2 s1 ] [b] = [t1] */
    /* [ s2 s1 s0 ] [c] = [t0] */
    double A[3][4] = {
        {s4, s3, s2, t2},
        {s3, s2, s1, t1},
        {s2, s1, s0, t0}
    };

    solve_3x3(A, coef);
}

double quadratic_eval(double coef[3], double x) {
    return coef[0] * x * x + coef[1] * x + coef[2];
}

/* ---------- 主程序 ---------- */

int main() {
    DataPoint data[TOTAL_YEARS] = {
        {2010, 1588},
        {2011, 1600},
        {2012, 1635},
        {2013, 1640},
        {2014, 1687},
        {2015, 1655},
        {2016, 1786},
        {2017, 1723},
        {2018, 1523},
        {2019, 1465},
        {2020, 1200},
        {2021, 1062},
        {2022, 956},
        {2023, 902},
        {2024, 954},
        {2025, 792}
    };

    FILE *fp_actual, *fp_all_pred, *fp_best_pred, *fp_gp;
    FILE *fp_fit_curve, *fp_fit_pred;
    Result results[CASES];
    int result_count = 0;
    int start, i, k;

    /* 输出真实数据点 */
    fp_actual = fopen("actual_data.dat", "w");
    if (fp_actual == NULL) {
        printf("Cannot open actual_data.dat\n");
        return 1;
    }
    for (i = 0; i < TOTAL_YEARS; i++) {
        fprintf(fp_actual, "%d %.6f\n", data[i].year, data[i].value);
    }
    fclose(fp_actual);

    /* 所有插值预测点 */
    fp_all_pred = fopen("all_predictions.dat", "w");
    if (fp_all_pred == NULL) {
        printf("Cannot open all_predictions.dat\n");
        return 1;
    }

    printf("Interpolation results:\n");
    printf("Range                N      Predicted_2026      Score\n");
    printf("----------------------------------------------------------\n");

    for (start = 2010; start <= 2024; start++) {
        int n = 2025 - start + 1;
        double x[MAXN], y[MAXN], coeff[MAXN];
        double pred2026;
        double score;
        char filename[64];
        FILE *fp_curve;

        for (i = 0; i < n; i++) {
            x[i] = (double)i;
            y[i] = data[(start - 2010) + i].value;
        }

        build_newton_coeff(x, y, n, coeff);
        pred2026 = newton_eval(x, coeff, n, (double)n);

        score = calc_score(pred2026, data[13].value, data[14].value, data[15].value);

        results[result_count].start_year = start;
        results[result_count].end_year = 2025;
        results[result_count].n = n;
        results[result_count].pred2026 = pred2026;
        results[result_count].score = score;
        result_count++;

        fprintf(fp_all_pred, "2026 %.6f \"%d-%d\"\n", pred2026, start, 2025);

        printf("%d-%d          %-6d %-18.6f %.6f\n",
               start, 2025, n, pred2026, score);

        sprintf(filename, "curve_%d_2025.dat", start);
        fp_curve = fopen(filename, "w");
        if (fp_curve == NULL) {
            printf("Cannot open %s\n", filename);
            fclose(fp_all_pred);
            return 1;
        }

        for (k = 0; k < CURVE_POINTS; k++) {
            double year = 2010.0 + 16.0 * k / (CURVE_POINTS - 1); /* 2010~2026 */
            double t = year - start;
            double val = newton_eval(x, coeff, n, t);
            fprintf(fp_curve, "%.6f %.6f\n", year, val);
        }
        fclose(fp_curve);
    }
    fclose(fp_all_pred);

    /* 相对不那么离谱的插值点 */
    {
        int top_idx[3];
        fp_best_pred = fopen("best_predictions.dat", "w");
        if (fp_best_pred == NULL) {
            printf("Cannot open best_predictions.dat\n");
            return 1;
        }

        find_best_k(results, result_count, top_idx, 3);

        printf("\nRelatively less extreme interpolation predictions:\n");
        for (i = 0; i < 3; i++) {
            int idx = top_idx[i];
            if (idx != -1) {
                fprintf(fp_best_pred, "2026 %.6f \"%d-%d\"\n",
                        results[idx].pred2026,
                        results[idx].start_year,
                        results[idx].end_year);

                printf("Top %d: %d-%d -> %.6f\n",
                       i + 1,
                       results[idx].start_year,
                       results[idx].end_year,
                       results[idx].pred2026);
            }
        }
        fclose(fp_best_pred);
    }

    /* ---------- 二次拟合 ---------- */
    {
        double x[TOTAL_YEARS], y[TOTAL_YEARS];
        double coef[3];
        double pred2026;

        for (i = 0; i < TOTAL_YEARS; i++) {
            x[i] = (double)i;   /* 2010->0, ..., 2025->15 */
            y[i] = data[i].value;
        }

        quadratic_fit(x, y, TOTAL_YEARS, coef);
        pred2026 = quadratic_eval(coef, 16.0); /* 2026 -> 16 */

        printf("\nQuadratic least-squares fitting result:\n");
        printf("y = %.6f x^2 + %.6f x + %.6f\n", coef[0], coef[1], coef[2]);
        printf("Predicted 2026 by quadratic fit = %.6f\n", pred2026);

        fp_fit_curve = fopen("quadratic_fit_curve.dat", "w");
        if (fp_fit_curve == NULL) {
            printf("Cannot open quadratic_fit_curve.dat\n");
            return 1;
        }

        for (k = 0; k < CURVE_POINTS; k++) {
            double year = 2010.0 + 16.0 * k / (CURVE_POINTS - 1);
            double t = year - 2010.0;
            double val = quadratic_eval(coef, t);
            fprintf(fp_fit_curve, "%.6f %.6f\n", year, val);
        }
        fclose(fp_fit_curve);

        fp_fit_pred = fopen("fit_prediction.dat", "w");
        if (fp_fit_pred == NULL) {
            printf("Cannot open fit_prediction.dat\n");
            return 1;
        }
        fprintf(fp_fit_pred, "2026 %.6f \"QuadraticFit\"\n", pred2026);
        fclose(fp_fit_pred);
    }

    /* ---------- 生成 gnuplot 脚本 ---------- */
    fp_gp = fopen("plot_population.gnuplot", "w");
    if (fp_gp == NULL) {
        printf("Cannot open plot_population.gnuplot\n");
        return 1;
    }

    fprintf(fp_gp, "set terminal pngcairo size 1500,900 enhanced font 'Arial,14'\n");
    fprintf(fp_gp, "set output 'interpolation_extrapolation.png'\n");
    fprintf(fp_gp, "set title 'Interpolation vs Fitting for 2026 Prediction'\n");
    fprintf(fp_gp, "set xlabel 'Year'\n");
    fprintf(fp_gp, "set ylabel 'Value'\n");
    fprintf(fp_gp, "set grid\n");
    fprintf(fp_gp, "set key outside\n");
    fprintf(fp_gp, "set xrange [2010:2026.8]\n");
    fprintf(fp_gp, "set yrange [0:2000]\n");
    fprintf(fp_gp, "set style data lines\n");

    fprintf(fp_gp, "plot \\\n");

    for (start = 2010; start <= 2024; start++) {
        fprintf(fp_gp,
                "'curve_%d_2025.dat' using 1:2 with lines lw 1 title '%d-2025 interp', \\\n",
                start, start);
    }

    fprintf(fp_gp,
            "'quadratic_fit_curve.dat' using 1:2 with lines lw 3 title 'Quadratic fit curve', \\\n");
    fprintf(fp_gp,
            "'actual_data.dat' using 1:2 with points pt 7 ps 1.8 title 'Actual data points', \\\n");
    fprintf(fp_gp,
            "'all_predictions.dat' using 1:2 with points pt 5 ps 1.4 title 'All interpolation predictions', \\\n");
    fprintf(fp_gp,
            "'all_predictions.dat' using 1:2:3 with labels offset 0.8,0.3 notitle, \\\n");
    fprintf(fp_gp,
            "'best_predictions.dat' using 1:2 with points pt 9 ps 2.0 title 'Less extreme interpolation predictions', \\\n");
    fprintf(fp_gp,
            "'best_predictions.dat' using 1:2:3 with labels offset 0.8,0.8 notitle, \\\n");
    fprintf(fp_gp,
            "'fit_prediction.dat' using 1:2 with points pt 11 ps 2.5 title 'Quadratic-fit prediction', \\\n");
    fprintf(fp_gp,
            "'fit_prediction.dat' using 1:2:3 with labels offset 0.8,1.2 notitle\n");

    fclose(fp_gp);

    printf("\nFiles generated:\n");
    printf("1. actual_data.dat\n");
    printf("2. curve_2010_2025.dat ... curve_2024_2025.dat\n");
    printf("3. all_predictions.dat\n");
    printf("4. best_predictions.dat\n");
    printf("5. quadratic_fit_curve.dat\n");
    printf("6. fit_prediction.dat\n");
    printf("7. plot_population.gnuplot\n");
    printf("8. interpolation_extrapolation.png (after running gnuplot)\n");

    system("gnuplot plot_population.gnuplot");

    printf("\nDone.\n");
    return 0;
}
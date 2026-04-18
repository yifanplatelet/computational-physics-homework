#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOTAL_YEARS 16
#define MAXN 20
#define CURVE_POINTS 321
#define CASES 15   /* 2010-2025 到 2024-2025，共 15 组 */

typedef struct {
    int year;
    double value;
} DataPoint;

typedef struct {
    int start_year;
    int end_year;
    int n;
    double pred2026;
    double score;   /* 越小表示相对越“不离谱” */
} Result;

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

/* 简单评分：越接近“最近趋势延续”，分数越低
   这里只是为了“相对比较”，不是判定真好坏 */
double calc_score(double pred, double y2023, double y2024, double y2025) {
    double recent_trend = y2025 + (y2025 - y2024);   /* 用最后一年斜率做粗略参考 */
    double smooth_ref = (y2023 + y2024 + y2025) / 3.0; /* 最近三年均值 */

    double s1 = fabs(pred - recent_trend);
    double s2 = fabs(pred - smooth_ref);

    if (pred < 0) {
        s1 += 1000000.0;
        s2 += 1000000.0;
    }

    return s1 + 0.5 * s2;
}

/* 选出分数最小的前k个 */
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
    Result results[CASES];
    int result_count = 0;
    int start, i, k;

    fp_actual = fopen("actual_data.dat", "w");
    if (fp_actual == NULL) {
        printf("Cannot open actual_data.dat\n");
        return 1;
    }
    for (i = 0; i < TOTAL_YEARS; i++) {
        fprintf(fp_actual, "%d %.6f\n", data[i].year, data[i].value);
    }
    fclose(fp_actual);

    fp_all_pred = fopen("all_predictions.dat", "w");
    if (fp_all_pred == NULL) {
        printf("Cannot open all_predictions.dat\n");
        return 1;
    }

    printf("Range                N      Predicted_2026      Score\n");
    printf("----------------------------------------------------------\n");

    /* 这里改成 2024，就会包含：
       2022-2025, 2023-2025, 2024-2025 */
    for (start = 2010; start <= 2024; start++) {
        int n = 2025 - start + 1;
        double x[MAXN], y[MAXN], coeff[MAXN];
        double pred2026;
        double score;
        char filename[64];
        FILE *fp_curve;

        for (i = 0; i < n; i++) {
            x[i] = (double)i;   /* 相对年份 */
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
            double t = year - start;  /* 相对该窗口起点 */
            double val = newton_eval(x, coeff, n, t);
            fprintf(fp_curve, "%.6f %.6f\n", year, val);
        }
        fclose(fp_curve);
    }
    fclose(fp_all_pred);

    {
        int top_idx[3];
        fp_best_pred = fopen("best_predictions.dat", "w");
        if (fp_best_pred == NULL) {
            printf("Cannot open best_predictions.dat\n");
            return 1;
        }

        find_best_k(results, result_count, top_idx, 3);

        printf("\nRelatively less extreme predictions:\n");
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

    fp_gp = fopen("plot_population.gnuplot", "w");
    if (fp_gp == NULL) {
        printf("Cannot open plot_population.gnuplot\n");
        return 1;
    }

    fprintf(fp_gp, "set terminal pngcairo size 1500,900 enhanced font 'Arial,14'\n");
    fprintf(fp_gp, "set output 'interpolation_extrapolation.png'\n");
    fprintf(fp_gp, "set title 'Interpolation Curves and 2026 Extrapolated Points'\n");
    fprintf(fp_gp, "set xlabel 'Year'\n");
    fprintf(fp_gp, "set ylabel 'Value'\n");
    fprintf(fp_gp, "set grid\n");
    fprintf(fp_gp, "set key outside\n");
    fprintf(fp_gp, "set xrange [2010:2026.8]\n");
    fprintf(fp_gp, "set style data lines\n");

    fprintf(fp_gp, "plot \\\n");
    for (start = 2010; start <= 2024; start++) {
        fprintf(fp_gp,
                "'curve_%d_2025.dat' using 1:2 with lines lw 1 title '%d-2025 fit', \\\n",
                start, start);
    }

    fprintf(fp_gp,
            "'actual_data.dat' using 1:2 with points pt 7 ps 1.6 title 'Actual data', \\\n");
    fprintf(fp_gp,
            "'all_predictions.dat' using 1:2 with points pt 5 ps 1.6 title 'All 2026 predictions', \\\n");
    fprintf(fp_gp,
            "'all_predictions.dat' using 1:2:3 with labels offset 1,0.5 notitle, \\\n");
    fprintf(fp_gp,
            "'best_predictions.dat' using 1:2 with points pt 9 ps 2.2 title 'Less extreme predictions', \\\n");
    fprintf(fp_gp,
            "'best_predictions.dat' using 1:2:3 with labels offset 1,1.2 notitle\n");

    fclose(fp_gp);

    printf("\nFiles generated:\n");
    printf("1. actual_data.dat\n");
    printf("2. curve_2010_2025.dat ... curve_2024_2025.dat\n");
    printf("3. all_predictions.dat\n");
    printf("4. best_predictions.dat\n");
    printf("5. plot_population.gnuplot\n");
    printf("6. interpolation_extrapolation.png (after running gnuplot)\n");

    system("gnuplot plot_population.gnuplot");

    printf("\nDone.\n");
    return 0;
}
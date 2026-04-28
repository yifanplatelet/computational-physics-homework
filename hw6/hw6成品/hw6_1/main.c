#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 7
#define M 1200   /* 绘图采样点数 */

/* 原始数据 */
double X[N] = {-1.00, 0.00, 1.27, 2.55, 3.82, 4.92, 5.02};
double Y[N] = {-14.58, 0.00, 0.00, 0.00, 0.00, 0.88, 11.17};

/* =========================================================
   1. 6次整体插值多项式：拉格朗日形式
   ========================================================= */
double poly6_eval(double x)
{
    int i, j;
    double sum = 0.0;

    for (i = 0; i < N; i++) {
        double term = Y[i];
        for (j = 0; j < N; j++) {
            if (j != i) {
                term *= (x - X[j]) / (X[i] - X[j]);
            }
        }
        sum += term;
    }
    return sum;
}

/* =========================================================
   2. 自然三次样条
   每段:
   S_i(x)=a[i]+b[i](x-X[i])+c[i](x-X[i])^2+d[i](x-X[i])^3
   ========================================================= */
void natural_spline_coeff(double a[], double b[], double c[], double d[])
{
    double h[N - 1], alpha[N], l[N], mu[N], z[N];
    int i, j;

    for (i = 0; i < N - 1; i++) {
        h[i] = X[i + 1] - X[i];
        a[i] = Y[i];
    }

    alpha[0] = 0.0;
    alpha[N - 1] = 0.0;

    for (i = 1; i < N - 1; i++) {
        alpha[i] = 3.0 / h[i] * (Y[i + 1] - Y[i])
                 - 3.0 / h[i - 1] * (Y[i] - Y[i - 1]);
    }

    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;

    for (i = 1; i < N - 1; i++) {
        l[i] = 2.0 * (X[i + 1] - X[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[N - 1] = 1.0;
    z[N - 1] = 0.0;
    c[N - 1] = 0.0;

    for (j = N - 2; j >= 0; j--) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (Y[j + 1] - Y[j]) / h[j]
             - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
    }
}

double natural_spline_eval(double x, double a[], double b[], double c[], double d[])
{
    int i;
    if (x <= X[0]) i = 0;
    else if (x >= X[N - 1]) i = N - 2;
    else {
        for (i = 0; i < N - 1; i++) {
            if (x >= X[i] && x <= X[i + 1]) break;
        }
    }

    double dx = x - X[i];
    return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
}

/* =========================================================
   3. 分段线性插值
   ========================================================= */
double linear_eval(double x)
{
    int i;
    if (x <= X[0]) return Y[0];
    if (x >= X[N - 1]) return Y[N - 1];

    for (i = 0; i < N - 1; i++) {
        if (x >= X[i] && x <= X[i + 1]) {
            double t = (x - X[i]) / (X[i + 1] - X[i]);
            return Y[i] + t * (Y[i + 1] - Y[i]);
        }
    }
    return 0.0;
}

/* =========================================================
   4. PCHIP 保形三次 Hermite 插值
   ========================================================= */
void pchip_slopes(double m[])
{
    double h[N - 1], delta[N - 1];
    int k;

    for (k = 0; k < N - 1; k++) {
        h[k] = X[k + 1] - X[k];
        delta[k] = (Y[k + 1] - Y[k]) / h[k];
    }

    /* 内点 */
    for (k = 1; k < N - 1; k++) {
        if (delta[k - 1] == 0.0 || delta[k] == 0.0 ||
            (delta[k - 1] > 0 && delta[k] < 0) ||
            (delta[k - 1] < 0 && delta[k] > 0)) {
            m[k] = 0.0;
        } else {
            double w1 = 2.0 * h[k] + h[k - 1];
            double w2 = h[k] + 2.0 * h[k - 1];
            m[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k]);
        }
    }

    /* 左端点 */
    {
        double m0 = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1]);

        if ((m0 > 0 && delta[0] < 0) || (m0 < 0 && delta[0] > 0))
            m0 = 0.0;
        else if (((delta[0] > 0 && delta[1] < 0) || (delta[0] < 0 && delta[1] > 0))
                  && fabs(m0) > fabs(3.0 * delta[0]))
            m0 = 3.0 * delta[0];

        m[0] = m0;
    }

    /* 右端点 */
    {
        double mn = ((2.0 * h[N - 2] + h[N - 3]) * delta[N - 2]
                    - h[N - 2] * delta[N - 3]) / (h[N - 2] + h[N - 3]);

        if ((mn > 0 && delta[N - 2] < 0) || (mn < 0 && delta[N - 2] > 0))
            mn = 0.0;
        else if (((delta[N - 2] > 0 && delta[N - 3] < 0) ||
                  (delta[N - 2] < 0 && delta[N - 3] > 0))
                 && fabs(mn) > fabs(3.0 * delta[N - 2]))
            mn = 3.0 * delta[N - 2];

        m[N - 1] = mn;
    }
}

double pchip_eval(double x)
{
    double m[N];
    double h, t;
    int i;

    pchip_slopes(m);

    if (x <= X[0]) i = 0;
    else if (x >= X[N - 1]) i = N - 2;
    else {
        for (i = 0; i < N - 1; i++) {
            if (x >= X[i] && x <= X[i + 1]) break;
        }
    }

    h = X[i + 1] - X[i];
    t = (x - X[i]) / h;

    /* Hermite 基函数 */
    double h00 =  2 * t * t * t - 3 * t * t + 1;
    double h10 =      t * t * t - 2 * t * t + t;
    double h01 = -2 * t * t * t + 3 * t * t;
    double h11 =      t * t * t -     t * t;

    return h00 * Y[i]
         + h10 * h * m[i]
         + h01 * Y[i + 1]
         + h11 * h * m[i + 1];
}

/* =========================================================
   5. 分段物理启发模型
   反向区：线性近似
   截止区：I=0
   正向导通区：指数模型
   ========================================================= */
double physics_eval(double x)
{
    const double x0 = 3.82;
    const double A = 0.005436635421263133;
    const double B = 7.26920018899648;

    if (x <= 0.0) {
        return 14.58 * x;
    } else if (x > 0.0 && x < x0) {
        return 0.0;
    } else {
        return A * (exp(B * (x - x0)) - 1.0);
    }
}

/* =========================================================
   输出绘图数据
   compare_data.txt:
   x poly6 spline linear pchip physics
   ========================================================= */
void write_data_and_plot_script()
{
    FILE *fp, *gp;
    double a[N - 1], b[N - 1], c[N], d[N - 1];
    int i;

    natural_spline_coeff(a, b, c, d);

    fp = fopen("compare_data.txt", "w");
    if (!fp) {
        printf("无法创建 compare_data.txt\n");
        exit(1);
    }

    {
        double xmin = X[0];
        double xmax = X[N - 1];
        double step = (xmax - xmin) / (M - 1);

        for (i = 0; i < M; i++) {
            double x = xmin + i * step;
            double y_poly   = poly6_eval(x);
            double y_spline = natural_spline_eval(x, a, b, c, d);
            double y_linear = linear_eval(x);
            double y_pchip  = pchip_eval(x);
            double y_phys   = physics_eval(x);

            fprintf(fp, "%.10f %.10f %.10f %.10f %.10f %.10f\n",
                    x, y_poly, y_spline, y_linear, y_pchip, y_phys);
        }
    }
    fclose(fp);

    fp = fopen("raw_points.txt", "w");
    if (!fp) {
        printf("无法创建 raw_points.txt\n");
        exit(1);
    }
    for (i = 0; i < N; i++) {
        fprintf(fp, "%.10f %.10f\n", X[i], Y[i]);
    }
    fclose(fp);

    gp = fopen("plot.gp", "w");
    if (!gp) {
        printf("无法创建 plot.gp\n");
        exit(1);
    }

    fprintf(gp, "set terminal pngcairo size 1200,800\n");
    fprintf(gp, "set output 'result.png'\n");
    fprintf(gp, "set title 'Zener diode V-I: fitting/interpolation comparison'\n");
    fprintf(gp, "set xlabel 'Voltage'\n");
    fprintf(gp, "set ylabel 'Current'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "plot \\\n");
    fprintf(gp, "'compare_data.txt' using 1:2 with lines lw 2 title 'Polynomial degree 6', \\\n");
    fprintf(gp, "'compare_data.txt' using 1:3 with lines lw 2 title 'Natural cubic spline', \\\n");
    fprintf(gp, "'compare_data.txt' using 1:4 with lines lw 2 title 'Piecewise linear', \\\n");
    fprintf(gp, "'compare_data.txt' using 1:5 with lines lw 2 title 'PCHIP', \\\n");
    fprintf(gp, "'compare_data.txt' using 1:6 with lines lw 2 title 'Physics-inspired', \\\n");
    fprintf(gp, "'raw_points.txt' using 1:2 with points pt 7 ps 1.5 title 'Data points'\n");

    fclose(gp);
}

/* =========================================================
   打印部分测试值
   ========================================================= */
void print_sample_values()
{
    double a[N - 1], b[N - 1], c[N], d[N - 1];
    double test_x[] = {-1.0, -0.5, 0.0, 1.27, 2.0, 3.82, 4.5, 4.92, 5.0, 5.02};
    int k, cnt = sizeof(test_x) / sizeof(test_x[0]);

    natural_spline_coeff(a, b, c, d);

    printf("x\t\tpoly6\t\tspline\t\tlinear\t\tpchip\t\tphysics\n");
    for (k = 0; k < cnt; k++) {
        double x = test_x[k];
        printf("%8.3f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\n",
               x,
               poly6_eval(x),
               natural_spline_eval(x, a, b, c, d),
               linear_eval(x),
               pchip_eval(x),
               physics_eval(x));
    }
}

int main()
{
    int i;
    double a[N - 1], b[N - 1], c[N], d[N - 1];

    natural_spline_coeff(a, b, c, d);

    printf("====================================================\n");
    printf("稳压二极管 V-I 数据多种拟合/插值比较\n");
    printf("====================================================\n\n");

    printf("原始数据点：\n");
    for (i = 0; i < N; i++) {
        printf("(%.2f, %.2f)\n", X[i], Y[i]);
    }

    printf("\n自然三次样条各分段函数：\n\n");
    for (i = 0; i < N - 1; i++) {
        printf("S%d(x), x in [%.2f, %.2f]\n", i, X[i], X[i + 1]);
        printf("= %.8f + %.8f*(x-%.2f) + %.8f*(x-%.2f)^2 + %.8f*(x-%.2f)^3\n\n",
               a[i], b[i], X[i], c[i], X[i], d[i], X[i]);
    }

    print_sample_values();
    write_data_and_plot_script();

    printf("\n已生成文件：\n");
    printf("1. compare_data.txt  -> 各方法绘图数据\n");
    printf("2. raw_points.txt    -> 原始数据点\n");
    printf("3. plot.gp           -> gnuplot 绘图脚本\n");
    printf("\n绘图命令：\n");
    printf("gnuplot plot.gp\n");
    printf("生成图像文件：result.png\n");

    return 0;
}
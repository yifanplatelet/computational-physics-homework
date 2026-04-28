#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUM_CASES 6
#define SAMPLE_POINTS 2000

/* 题目要求的节点数 */
int node_counts[NUM_CASES] = {5, 7, 9, 17, 19, 21};

/* Runge 函数 */
double runge_function(double x)
{
    return 1.0 / (1.0 + 16.0 * x * x);
}

/* 生成 [-1,1] 上 n 个等距节点 */
void generate_nodes(int n, double x[], double y[])
{
    int i;
    double h = 2.0 / (n - 1);

    for (i = 0; i < n; i++) {
        x[i] = -1.0 + i * h;
        y[i] = runge_function(x[i]);
    }
}

/* 拉格朗日插值 */
double lagrange_interpolation(double t, int n, double x[], double y[])
{
    int i, j;
    double p = 0.0;

    for (i = 0; i < n; i++) {
        double L = 1.0;
        for (j = 0; j < n; j++) {
            if (j != i) {
                L *= (t - x[j]) / (x[i] - x[j]);
            }
        }
        p += y[i] * L;
    }

    return p;
}

/* 为某个节点数 n 生成数据文件 */
void write_case_data(int n)
{
    int i;
    double *x = (double *)malloc(n * sizeof(double));
    double *y = (double *)malloc(n * sizeof(double));
    FILE *fp_curve, *fp_nodes;
    char curve_file[64], node_file[64];

    if (x == NULL || y == NULL) {
        printf("内存分配失败\n");
        exit(1);
    }

    generate_nodes(n, x, y);

    sprintf(curve_file, "runge_%d_curve.txt", n);
    sprintf(node_file,  "runge_%d_nodes.txt", n);

    fp_curve = fopen(curve_file, "w");
    if (!fp_curve) {
        printf("无法创建文件 %s\n", curve_file);
        free(x);
        free(y);
        exit(1);
    }

    fp_nodes = fopen(node_file, "w");
    if (!fp_nodes) {
        printf("无法创建文件 %s\n", node_file);
        fclose(fp_curve);
        free(x);
        free(y);
        exit(1);
    }

    /* 保存节点 */
    for (i = 0; i < n; i++) {
        fprintf(fp_nodes, "%.12f %.12f\n", x[i], y[i]);
    }

    /* 保存原函数与插值函数的采样结果 */
    for (i = 0; i < SAMPLE_POINTS; i++) {
        double t = -1.0 + 2.0 * i / (SAMPLE_POINTS - 1);
        double f = runge_function(t);
        double p = lagrange_interpolation(t, n, x, y);
        fprintf(fp_curve, "%.12f %.12f %.12f\n", t, f, p);
    }

    fclose(fp_curve);
    fclose(fp_nodes);
    free(x);
    free(y);
}

/* 生成 gnuplot 绘图脚本 */
void write_gnuplot_script(void)
{
    FILE *gp = fopen("plot_runge.gp", "w");
    if (!gp) {
        printf("无法创建 gnuplot 脚本\n");
        exit(1);
    }

    fprintf(gp, "set terminal pngcairo size 1400,1800\n");
    fprintf(gp, "set output 'runge_comparison.png'\n");
    fprintf(gp, "set multiplot layout 3,2 title 'Runge Function: Exact vs Polynomial Interpolation'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "set key outside\n");

    for (int k = 0; k < NUM_CASES; k++) {
        int n = node_counts[k];
        fprintf(gp, "set title 'n = %d nodes, degree = %d'\n", n, n - 1);
        fprintf(gp,
                "plot 'runge_%d_curve.txt' using 1:2 with lines lw 2 title 'Exact f(x)', "
                "'runge_%d_curve.txt' using 1:3 with lines lw 2 title 'Interpolation', "
                "'runge_%d_nodes.txt' using 1:2 with points pt 7 ps 1.2 title 'Nodes'\n",
                n, n, n);
    }

    fprintf(gp, "unset multiplot\n");
    fclose(gp);
}

/* 打印部分节点信息 */
void print_nodes_info(int n)
{
    int i;
    double *x = (double *)malloc(n * sizeof(double));
    double *y = (double *)malloc(n * sizeof(double));

    if (x == NULL || y == NULL) {
        printf("内存分配失败\n");
        exit(1);
    }

    generate_nodes(n, x, y);

    printf("n = %d 个等距节点：\n", n);
    for (i = 0; i < n; i++) {
        printf("x[%2d] = % .6f,  f(x) = %.6f\n", i, x[i], y[i]);
    }
    printf("\n");

    free(x);
    free(y);
}

int main(void)
{
    int k;

    printf("=============================================\n");
    printf("Runge 函数高次插值实验\n");
    printf("f(x) = 1 / (1 + 16x^2)\n");
    printf("区间: [-1, 1]\n");
    printf("=============================================\n\n");

    for (k = 0; k < NUM_CASES; k++) {
        int n = node_counts[k];
        print_nodes_info(n);
        write_case_data(n);
    }

    write_gnuplot_script();

    printf("已生成以下文件：\n");
    for (k = 0; k < NUM_CASES; k++) {
        int n = node_counts[k];
        printf("  runge_%d_curve.txt\n", n);
        printf("  runge_%d_nodes.txt\n", n);
    }
    printf("  plot_runge.gp\n\n");

    printf("编译运行后，执行下面命令绘图：\n");
    printf("  gnuplot plot_runge.gp\n\n");
    printf("生成图像文件：\n");
    printf("  runge_comparison.png\n");

    return 0;
}
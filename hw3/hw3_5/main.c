#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288L
#endif

/* --------------------------------------------------
 * 原始递推：故意保留，用来观察误差恶化
 * -------------------------------------------------- */
void compute_naive(int N, double *z) {
    for (int i = 0; i <= N; ++i) z[i] = 0.0;
    z[2] = 2.0;

    for (int n = 2; n < N; ++n) {
        double a = ldexp(z[n] * z[n], 2 - 2 * n);   /* a = 4^(1-n) z_n^2 */
        double inner = 1.0 - sqrt(1.0 - a);         /* 不稳定 */
        z[n + 1] = pow(2.0, n - 0.5) * sqrt(inner);
    }
}

/* --------------------------------------------------
 * 稳定递推：作为参考值
 * z_{n+1} = sqrt(2) * z_n / sqrt(1 + sqrt(1-a_n))
 * -------------------------------------------------- */
void compute_stable(int N, long double *z) {
    for (int i = 0; i <= N; ++i) z[i] = 0.0L;
    z[2] = 2.0L;

    for (int n = 2; n < N; ++n) {
        long double a = ldexpl(z[n] * z[n], 2 - 2 * n);
        long double s = sqrtl(1.0L - a);
        z[n + 1] = sqrtl(2.0L) * z[n] / sqrtl(1.0L + s);
    }
}

/* --------------------------------------------------
 * 稳定计算 t = 1 - sqrt(1-a)
 * 用有理化避免消去
 * -------------------------------------------------- */
long double stable_t(long double a) {
    long double s = sqrtl(1.0L - a);
    return a / (1.0L + s);
}

int main(void) {
    int N;
    printf("请输入 N (建议 40~60): ");
    if (scanf("%d", &N) != 1 || N < 2) {
        printf("输入无效。\n");
        return 1;
    }

    double *z_naive = (double *)malloc((N + 1) * sizeof(double));
    long double *z_stable = (long double *)malloc((N + 1) * sizeof(long double));

    if (!z_naive || !z_stable) {
        printf("内存分配失败。\n");
        free(z_naive);
        free(z_stable);
        return 1;
    }

    compute_naive(N, z_naive);
    compute_stable(N, z_stable);

    const long double eps = (long double)DBL_EPSILON;

    printf("\n==================== Part 1 ====================\n");
    printf("原始递推 vs 稳定递推\n\n");
    printf("%3s %20s %20s %18s %18s\n",
           "n", "z_naive", "z_stable", "|naive-pi|", "|stable-pi|");

    for (int n = 2; n <= N; ++n) {
        long double err_naive_pi = fabsl((long double)z_naive[n] - M_PI);
        long double err_stable_pi = fabsl(z_stable[n] - M_PI);

        printf("%3d %20.12Lf %20.12Lf %18.10Le %18.10Le\n",
               n,
               (long double)z_naive[n],
               z_stable[n],
               err_naive_pi,
               err_stable_pi);
    }

    printf("\n==================== Part 2 ====================\n");
    printf("真实误差 + 稳定误差估计\n\n");
    printf("%3s %18s %18s %18s %18s %18s\n",
           "n", "a_n", "t_n(stable)", "real_abs_err",
           "rel_local_est", "log10(rel_est)");

    for (int n = 2; n < N; ++n) {
        long double a = ldexpl(z_stable[n] * z_stable[n], 2 - 2 * n);
        long double t = stable_t(a);

        long double real_abs_err = fabsl((long double)z_naive[n + 1] - z_stable[n + 1]);

        /* 局部相对误差估计：
           rel(z_{n+1}) ~ 0.5 * eps / t ~ eps / a
        */
        long double rel_local_est = 0.5L * eps / t;
        long double log_rel = log10l(rel_local_est);

        printf("%3d %18.10Le %18.10Le %18.10Le %18.10Le %18.10Lf\n",
               n, a, t, real_abs_err, rel_local_est, log_rel);
    }

    printf("\n==================== Part 3 ====================\n");
    printf("与渐近公式比较\n\n");
    printf("%3s %18s %18s %18s\n",
           "n", "rel_local_est", "eps*4^(n-1)/pi^2", "ratio");

    for (int n = 2; n < N; ++n) {
        long double a = ldexpl(z_stable[n] * z_stable[n], 2 - 2 * n);
        long double t = stable_t(a);
        long double rel_local_est = 0.5L * eps / t;

        long double asym = eps * powl(4.0L, (long double)(n - 1)) / (M_PI * M_PI);
        long double ratio = rel_local_est / asym;

        printf("%3d %18.10Le %18.10Le %18.10Lf\n",
               n, rel_local_est, asym, ratio);
    }

    printf("\n==================== 结论 ====================\n");
    printf("1. 原始递推在 n 较大时因 1 - sqrt(1-a) 发生灾难性消去而失真。\n");
    printf("2. 稳定递推用有理化形式，避免了这个问题，因此可稳定逼近 pi。\n");
    printf("3. 误差分析中若也直接算 1 - sqrt(1-a)，模型本身会先坏掉；\n");
    printf("   必须同样使用稳定公式 t = a / (1 + sqrt(1-a))。\n");
    printf("4. 局部相对误差量级约为 eps / a_n，随 n 近似按 4^(n-1) 增长。\n");

    free(z_naive);
    free(z_stable);
    return 0;
}
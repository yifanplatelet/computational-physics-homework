#include <stdio.h>
#include <math.h>

/* 相对误差 */
static long double rel_err(long double approx, long double exact) {
    return fabsl((approx - exact) / exact);
}

/* 参考值：用 long double 计算，并用稳定方式得到小根 */
static long double x2_ref(long double b) {
    long double r = sqrtl(b * b - 4.0L);
    long double x1 = (b + r) / 2.0L;   // 大根，稳定
    return 1.0L / x1;                  // 因为 x1*x2 = 1
}

/* ---------- double 版本 ---------- */

/* 标准公式计算小根 */
static double x2_std_double(double b) {
    double r = sqrt(b * b - 4.0);
    return (b - r) / 2.0;
}

/* 有理化公式计算小根 */
static double x2_rat_double(double b) {
    double r = sqrt(b * b - 4.0);
    return 2.0 / (b + r);
}

/* ---------- float 版本 ---------- */

/* 标准公式计算小根 */
static float x2_std_float(float b) {
    float r = sqrtf(b * b - 4.0f);
    return (b - r) / 2.0f;
}

/* 有理化公式计算小根 */
static float x2_rat_float(float b) {
    float r = sqrtf(b * b - 4.0f);
    return 2.0f / (b + r);
}

int main(void) {
    long double bs[] = {
        100.0L, 1000.0L, 10000.0L, 100000.0L,
        1000000.0L, 10000000.0L, 100000000.0L
    };
    int n = (int)(sizeof(bs) / sizeof(bs[0]));

    /* ---------- double 实验 ---------- */
    printf("===== double 实验 =====\n");
    printf("%12s %20s %20s %14s %14s\n",
           "b", "x2_standard", "x2_rational", "relerr_std", "relerr_rat");

    for (int i = 0; i < n; ++i) {
        long double b = bs[i];
        long double ref = x2_ref(b);

        double xs = x2_std_double((double)b);
        double xr = x2_rat_double((double)b);

        printf("%12.0Lf %20.12e %20.12e %14.6Le %14.6Le\n",
               b, xs, xr,
               rel_err((long double)xs, ref),
               rel_err((long double)xr, ref));
    }

    printf("\n");

    /* ---------- float 实验 ---------- */
    printf("===== float 实验（误差更早暴露） =====\n");
    printf("%12s %20s %20s %14s %14s\n",
           "b", "x2_standard", "x2_rational", "relerr_std", "relerr_rat");

    for (int i = 0; i < n; ++i) {
        long double b = bs[i];
        long double ref = x2_ref(b);

        float xs = x2_std_float((float)b);
        float xr = x2_rat_float((float)b);

        printf("%12.0Lf %20.12e %20.12e %14.6Le %14.6Le\n",
               b, xs, xr,
               rel_err((long double)xs, ref),
               rel_err((long double)xr, ref));
    }

    return 0;
}
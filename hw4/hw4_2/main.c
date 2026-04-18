#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <quadmath.h>

typedef __float128 quad_t;

/* (x-1)^10 = x^10 -10x^9 +45x^8 -120x^7 +210x^6 -252x^5
              +210x^4 -120x^3 +45x^2 -10x +1 */
static const int C[11] = {1, -10, 45, -120, 210, -252, 210, -120, 45, -10, 1};

/* ---------- 直接计算 (x-1)^10 ---------- */
static float direct_f(float x) {
    float t = x - 1.0f;
    float y = 1.0f;
    for (int i = 0; i < 10; ++i) y *= t;
    return y;
}

static double direct_d(double x) {
    double t = x - 1.0;
    double y = 1.0;
    for (int i = 0; i < 10; ++i) y *= t;
    return y;
}

static quad_t direct_q(quad_t x) {
    quad_t t = x - (quad_t)1.0;
    quad_t y = (quad_t)1.0;
    for (int i = 0; i < 10; ++i) y *= t;
    return y;
}

/* ---------- 直接按展开式计算 ---------- */
static float expand_f(float x) {
    float p[11];
    p[0] = 1.0f;
    for (int i = 1; i <= 10; ++i) p[i] = p[i - 1] * x; /* p[i] = x^i */

    float y = 0.0f;
    for (int i = 0; i <= 10; ++i) {
        y += (float)C[i] * p[10 - i];
    }
    return y;
}

static double expand_d(double x) {
    double p[11];
    p[0] = 1.0;
    for (int i = 1; i <= 10; ++i) p[i] = p[i - 1] * x;

    double y = 0.0;
    for (int i = 0; i <= 10; ++i) {
        y += (double)C[i] * p[10 - i];
    }
    return y;
}

static quad_t expand_q(quad_t x) {
    quad_t p[11];
    p[0] = (quad_t)1.0;
    for (int i = 1; i <= 10; ++i) p[i] = p[i - 1] * x;

    quad_t y = (quad_t)0.0;
    for (int i = 0; i <= 10; ++i) {
        y += (quad_t)C[i] * p[10 - i];
    }
    return y;
}

/* ---------- Horner 方法 ---------- */
static float horner_f(float x) {
    float y = (float)C[0];
    for (int i = 1; i <= 10; ++i) {
        y = y * x + (float)C[i];
    }
    return y;
}

static double horner_d(double x) {
    double y = (double)C[0];
    for (int i = 1; i <= 10; ++i) {
        y = y * x + (double)C[i];
    }
    return y;
}

static quad_t horner_q(quad_t x) {
    quad_t y = (quad_t)C[0];
    for (int i = 1; i <= 10; ++i) {
        y = y * x + (quad_t)C[i];
    }
    return y;
}

/* ---------- 工具函数 ---------- */
static void qtoa(char *buf, size_t n, quad_t x) {
    quadmath_snprintf(buf, n, "%.36Qg", x);
}

static void fprint_q(FILE *fp, quad_t x) {
    char buf[128];
    qtoa(buf, sizeof(buf), x);
    fputs(buf, fp);
}

static quad_t relerr_q(quad_t approx, quad_t ref) {
    if (ref == 0) return NAN;  /* x=1 时真值为 0，相对误差无定义 */
    return fabsq((approx - ref) / ref);
}

int main(void) {
    const int N = 6000;              /* 6001 个点 */
    const double xmin = 0.7;
    const double xmax = 1.3;

    FILE *fv = fopen("values.csv", "w");
    FILE *fe = fopen("relerr.csv", "w");
    if (!fv || !fe) {
        fprintf(stderr, "Cannot open output files.\n");
        return 1;
    }

    fprintf(fv,
        "x,"
        "direct_f,expand_f,horner_f,"
        "direct_d,expand_d,horner_d,"
        "direct_q,expand_q,horner_q\n");

    fprintf(fe,
        "x,"
        "rel_direct_f,rel_expand_f,rel_horner_f,"
        "rel_direct_d,rel_expand_d,rel_horner_d,"
        "rel_direct_q,rel_expand_q,rel_horner_q\n");

    quad_t max_rel[9] = {0}; /* f:3, d:3, q:3 */

    for (int i = 0; i <= N; ++i) {
        double x = xmin + (xmax - xmin) * (double)i / (double)N;

        float  xf = (float)x;
        double xd = x;
        quad_t xq = (quad_t)x;

        /* 三种算法 */
        float  df = direct_f(xf);
        float  ef = expand_f(xf);
        float  hf = horner_f(xf);

        double dd = direct_d(xd);
        double ed = expand_d(xd);
        double hd = horner_d(xd);

        quad_t dq = direct_q(xq);
        quad_t eq = expand_q(xq);
        quad_t hq = horner_q(xq);

        /* 用 quad 的直接计算作为参考 */
        quad_t ref = dq;

        quad_t r_df = relerr_q((quad_t)df, ref);
        quad_t r_ef = relerr_q((quad_t)ef, ref);
        quad_t r_hf = relerr_q((quad_t)hf, ref);

        quad_t r_dd = relerr_q((quad_t)dd, ref);
        quad_t r_ed = relerr_q((quad_t)ed, ref);
        quad_t r_hd = relerr_q((quad_t)hd, ref);

        quad_t r_dq = relerr_q(dq, ref);
        quad_t r_eq = relerr_q(eq, ref);
        quad_t r_hq = relerr_q(hq, ref);

        /* 写 values.csv */
        fprintf(fv, "%.17g,", x);
        fprintf(fv, "%.9g,%.9g,%.9g,", df, ef, hf);
        fprintf(fv, "%.17g,%.17g,%.17g,", dd, ed, hd);
        fprint_q(fv, dq); fputc(',', fv);
        fprint_q(fv, eq); fputc(',', fv);
        fprint_q(fv, hq); fputc('\n', fv);

        /* 写 relerr.csv */
        fprintf(fe, "%.17g,", x);
        fprint_q(fe, r_df); fputc(',', fe);
        fprint_q(fe, r_ef); fputc(',', fe);
        fprint_q(fe, r_hf); fputc(',', fe);
        fprint_q(fe, r_dd); fputc(',', fe);
        fprint_q(fe, r_ed); fputc(',', fe);
        fprint_q(fe, r_hd); fputc(',', fe);
        fprint_q(fe, r_dq); fputc(',', fe);
        fprint_q(fe, r_eq); fputc(',', fe);
        fprint_q(fe, r_hq); fputc('\n', fe);

        /* 更新最大相对误差（跳过 NaN） */
        quad_t arr[9] = {r_df, r_ef, r_hf, r_dd, r_ed, r_hd, r_dq, r_eq, r_hq};
        for (int k = 0; k < 9; ++k) {
            if (!isnan((double)arr[k]) && arr[k] > max_rel[k]) {
                max_rel[k] = arr[k];
            }
        }
    }

    fclose(fv);
    fclose(fe);

    const char *name[9] = {
        "direct_f", "expand_f", "horner_f",
        "direct_d", "expand_d", "horner_d",
        "direct_q", "expand_q", "horner_q"
    };

    printf("Done.\n");
    printf("Generated files:\n");
    printf("  values.csv   -> 用来画函数曲线\n");
    printf("  relerr.csv   -> 用来画相对误差分布\n\n");

    printf("Maximum relative errors (vs direct_q):\n");
    for (int i = 0; i < 9; ++i) {
        char buf[128];
        qtoa(buf, sizeof(buf), max_rel[i]);
        printf("  %-10s : %s\n", name[i], buf);
    }

    printf("\nNote:\n");
    printf("  x = 1 时真值为 0，所以相对误差无定义，CSV 中会写成 nan。\n");
    printf("  quad 使用 __float128，需要 GCC/Clang + libquadmath。\n");

    return 0;
}
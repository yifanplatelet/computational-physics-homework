#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <quadmath.h>

static void print_q(const char *label, __float128 x) {
    char buf[128];
    quadmath_snprintf(buf, sizeof(buf), "%.36Qg", x);
    printf("%-24s = %s\n", label, buf);
}

static void sort2f(float *x, float *y) {
    if (*x > *y) {
        float t = *x; *x = *y; *y = t;
    }
}

static void sort2d(double *x, double *y) {
    if (*x > *y) {
        double t = *x; *x = *y; *y = t;
    }
}

static void sort2q(__float128 *x, __float128 *y) {
    if (*x > *y) {
        __float128 t = *x; *x = *y; *y = t;
    }
}

/* Round to a given number of decimal significant digits. */
static double round_sig(double x, int sig) {
    if (x == 0.0) {
        return 0.0;
    }

    double ax = fabs(x);
    double e = floor(log10(ax));
    double scale = pow(10.0, sig - 1 - e);
    return round(x * scale) / scale;
}

/* Simulate fl(x op y) in 3-significant-digit decimal arithmetic. */
static double fl_add(double x, double y) { return round_sig(x + y, 3); }
static double fl_sub(double x, double y) { return round_sig(x - y, 3); }
static double fl_mul(double x, double y) { return round_sig(x * y, 3); }

static void discriminant_three_digit_demo(void) {
    double a = 1.22, b = 3.34, c = 2.28;
    double exact_disc = b * b - 4.0 * a * c;

    double b2_exact = b * b;
    double b2_fl = round_sig(b2_exact, 3);

    /* Compute 4ac with rounding after each multiplication. */
    double foura_fl = fl_mul(4.0, a);
    double fourac_fl = fl_mul(foura_fl, c);
    double fourac_exact = 4.0 * a * c;

    double disc_fl = fl_sub(b2_fl, fourac_fl);
    double abs_err = fabs(disc_fl - exact_disc);
    double rel_err = abs_err / fabs(exact_disc);

    printf("\n[3-significant-digit arithmetic demo for b^2 - 4ac]\n");
    printf("a = %.2f, b = %.2f, c = %.2f\n", a, b, c);
    printf("exact b^2                 = %.10g\n", b2_exact);
    printf("fl(b^2)                   = %.10g\n", b2_fl);
    printf("exact 4ac                 = %.10g\n", fourac_exact);
    printf("fl(4a)                    = %.10g\n", foura_fl);
    printf("fl(fl(4a)*c) = fl(4ac)    = %.10g\n", fourac_fl);
    printf("exact b^2 - 4ac           = %.10g\n", exact_disc);
    printf("fl(fl(b^2) - fl(4ac))     = %.10g\n", disc_fl);
    printf("absolute error            = %.10g\n", abs_err);
    printf("relative error            = %.10g  (%.2f%%)\n", rel_err, 100.0 * rel_err);

    printf("\nInterpretation:\n");
    printf("  b^2  = 11.1556  -> 11.2  (3 significant digits)\n");
    printf("  4ac  = 11.1264  -> 11.1  (3 significant digits)\n");
    printf("  11.2 - 11.1 = 0.1, while the exact value is 0.0292\n");
    printf("  This is a classic cancellation error.\n");
}

static float direct7f(float x) {
    float t = x - 1.0f;
    float t2 = t * t;
    float t4 = t2 * t2;
    return t4 * t2 * t;
}

static double direct7d(double x) {
    double t = x - 1.0;
    double t2 = t * t;
    double t4 = t2 * t2;
    return t4 * t2 * t;
}

static __float128 direct7q(__float128 x) {
    __float128 t = x - 1.0Q;
    __float128 t2 = t * t;
    __float128 t4 = t2 * t2;
    return t4 * t2 * t;
}

/* Expanded polynomial written explicitly (not Horner),
 * so the cancellation is easier to observe near x = 1.
 */
static float expanded7f(float x) {
    float x2 = x * x;
    float x3 = x2 * x;
    float x4 = x3 * x;
    float x5 = x4 * x;
    float x6 = x5 * x;
    float x7 = x6 * x;
    return x7 - 7.0f * x6 + 21.0f * x5 - 35.0f * x4 +
           35.0f * x3 - 21.0f * x2 + 7.0f * x - 1.0f;
}

static double expanded7d(double x) {
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double x5 = x4 * x;
    double x6 = x5 * x;
    double x7 = x6 * x;
    return x7 - 7.0 * x6 + 21.0 * x5 - 35.0 * x4 +
           35.0 * x3 - 21.0 * x2 + 7.0 * x - 1.0;
}

static __float128 expanded7q(__float128 x) {
    __float128 x2 = x * x;
    __float128 x3 = x2 * x;
    __float128 x4 = x3 * x;
    __float128 x5 = x4 * x;
    __float128 x6 = x5 * x;
    __float128 x7 = x6 * x;
    return x7 - 7.0Q * x6 + 21.0Q * x5 - 35.0Q * x4 +
           35.0Q * x3 - 21.0Q * x2 + 7.0Q * x - 1.0Q;
}

static void quadratic_float(void) {
    float a = 1.22f, b = 3.34f, c = 2.28f;
    float disc = b * b - 4.0f * a * c;
    float s = sqrtf(disc);

    float r1_std = (-b + s) / (2.0f * a);
    float r2_std = (-b - s) / (2.0f * a);

    float q = -0.5f * (b + (b >= 0.0f ? s : -s));
    float r1_stable = q / a;
    float r2_stable = c / q;

    sort2f(&r1_std, &r2_std);
    sort2f(&r1_stable, &r2_stable);

    printf("\n[float]\n");
    printf("discriminant              = %.9g\n", disc);
    printf("standard root 1           = %.9g\n", r1_std);
    printf("standard root 2           = %.9g\n", r2_std);
    printf("stable   root 1           = %.9g\n", r1_stable);
    printf("stable   root 2           = %.9g\n", r2_stable);
}

static void quadratic_double(void) {
    double a = 1.22, b = 3.34, c = 2.28;
    double disc = b * b - 4.0 * a * c;
    double s = sqrt(disc);

    double r1_std = (-b + s) / (2.0 * a);
    double r2_std = (-b - s) / (2.0 * a);

    double q = -0.5 * (b + (b >= 0.0 ? s : -s));
    double r1_stable = q / a;
    double r2_stable = c / q;

    sort2d(&r1_std, &r2_std);
    sort2d(&r1_stable, &r2_stable);

    printf("\n[double]\n");
    printf("discriminant              = %.17g\n", disc);
    printf("standard root 1           = %.17g\n", r1_std);
    printf("standard root 2           = %.17g\n", r2_std);
    printf("stable   root 1           = %.17g\n", r1_stable);
    printf("stable   root 2           = %.17g\n", r2_stable);
}

static void quadratic_quad(void) {
    __float128 a = 1.22Q, b = 3.34Q, c = 2.28Q;
    __float128 disc = b * b - 4.0Q * a * c;
    __float128 s = sqrtq(disc);

    __float128 r1_std = (-b + s) / (2.0Q * a);
    __float128 r2_std = (-b - s) / (2.0Q * a);

    __float128 q = -0.5Q * (b + (b >= 0.0Q ? s : -s));
    __float128 r1_stable = q / a;
    __float128 r2_stable = c / q;

    sort2q(&r1_std, &r2_std);
    sort2q(&r1_stable, &r2_stable);

    printf("\n[quad: __float128]\n");
    print_q("discriminant", disc);
    print_q("standard root 1", r1_std);
    print_q("standard root 2", r2_std);
    print_q("stable   root 1", r1_stable);
    print_q("stable   root 2", r2_stable);
}

static void sample_point_report(void) {
    float xf = 1.0f + 1.0e-3f;
    double xd = 1.0 + 1.0e-9;
    __float128 xq = 1.0Q + 1.0e-12Q;

    printf("\nRepresentative sample points near x = 1:\n");
    printf("float  : x = %.9g, direct = %.9g, expanded = %.9g\n",
           xf, direct7f(xf), expanded7f(xf));
    printf("double : x = %.17g, direct = %.17g, expanded = %.17g\n",
           xd, direct7d(xd), expanded7d(xd));

    char xb[128], d1[128], d2[128];
    quadmath_snprintf(xb, sizeof(xb), "%.36Qg", xq);
    quadmath_snprintf(d1, sizeof(d1), "%.36Qg", direct7q(xq));
    quadmath_snprintf(d2, sizeof(d2), "%.36Qg", expanded7q(xq));
    printf("quad   : x = %s, direct = %s, expanded = %s\n", xb, d1, d2);
}

static void write_csv(const char *filename, long n, long double left, long double right) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror(filename);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "x,float_direct,float_expanded,double_direct,double_expanded,quad_direct,quad_expanded\n");

    for (long i = 0; i < n; ++i) {
        long double xl = left + (right - left) * (long double)i / (long double)(n - 1);

        float xf = (float)xl;
        double xd = (double)xl;
        __float128 xq = (__float128)xl;

        char qd[128], qe[128];
        quadmath_snprintf(qd, sizeof(qd), "%.36Qg", direct7q(xq));
        quadmath_snprintf(qe, sizeof(qe), "%.36Qg", expanded7q(xq));

        fprintf(fp,
                "%.21Lg,%.9g,%.9g,%.17g,%.17g,%s,%s\n",
                xl,
                direct7f(xf), expanded7f(xf),
                direct7d(xd), expanded7d(xd),
                qd, qe);
    }

    fclose(fp);
}

static void write_gnuplot_script(void) {
    FILE *gp = fopen("plot.gp", "w");
    if (!gp) {
        perror("plot.gp");
        exit(EXIT_FAILURE);
    }

    fprintf(gp,
        "set datafile separator ','\n"
        "set term pngcairo size 1400,900\n"
        "set grid\n"
        "set key outside\n"
        "\n"
        "set output 'zoom1_float.png'\n"
        "set title 'Zoom 1: x in [0.99, 1.01] (float)'\n"
        "plot 'poly_zoom1.csv' using 1:2 with lines title '(x-1)^7', \\\n"
        "     'poly_zoom1.csv' using 1:3 with lines title 'expanded'\n"
        "\n"
        "set output 'zoom1_double.png'\n"
        "set title 'Zoom 1: x in [0.99, 1.01] (double)'\n"
        "plot 'poly_zoom1.csv' using 1:4 with lines title '(x-1)^7', \\\n"
        "     'poly_zoom1.csv' using 1:5 with lines title 'expanded'\n"
        "\n"
        "set output 'zoom1_quad.png'\n"
        "set title 'Zoom 1: x in [0.99, 1.01] (quad)'\n"
        "plot 'poly_zoom1.csv' using 1:6 with lines title '(x-1)^7', \\\n"
        "     'poly_zoom1.csv' using 1:7 with lines title 'expanded'\n"
        "\n"
        "set output 'zoom2_float.png'\n"
        "set title 'Zoom 2: x in [0.9999, 1.0001] (float)'\n"
        "plot 'poly_zoom2.csv' using 1:2 with lines title '(x-1)^7', \\\n"
        "     'poly_zoom2.csv' using 1:3 with lines title 'expanded'\n"
        "\n"
        "set output 'zoom2_double.png'\n"
        "set title 'Zoom 2: x in [0.9999, 1.0001] (double)'\n"
        "plot 'poly_zoom2.csv' using 1:4 with lines title '(x-1)^7', \\\n"
        "     'poly_zoom2.csv' using 1:5 with lines title 'expanded'\n"
        "\n"
        "set output 'zoom2_quad.png'\n"
        "set title 'Zoom 2: x in [0.9999, 1.0001] (quad)'\n"
        "plot 'poly_zoom2.csv' using 1:6 with lines title '(x-1)^7', \\\n"
        "     'poly_zoom2.csv' using 1:7 with lines title 'expanded'\n");

    fclose(gp);
}

int main(void) {
    printf("Quadratic equation: 1.22 x^2 + 3.34 x + 2.28 = 0\n");

    discriminant_three_digit_demo();

    quadratic_float();
    quadratic_double();
    quadratic_quad();

    sample_point_report();

    write_csv("poly_zoom1.csv", 4001, 0.99L, 1.01L);
    write_csv("poly_zoom2.csv", 4001, 0.9999L, 1.0001L);
    write_gnuplot_script();

    printf("\nGenerated files:\n");
    printf("  poly_zoom1.csv   (x in [0.99,   1.01])\n");
    printf("  poly_zoom2.csv   (x in [0.9999, 1.0001])\n");
    printf("  plot.gp          (gnuplot script)\n");
    printf("\nRun: gnuplot plot.gp\n");

    return 0;
}
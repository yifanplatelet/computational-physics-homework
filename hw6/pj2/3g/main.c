#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mpfr.h>

#define DIGITS 1000000
#define FORCE_M     -1   // -1: 自动估计；>=3: 手动指定
#define FORCE_TERMS -1   // -1: 自动估计；>=1: 手动指定
#define EXTRA_TERMS 4

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static int auto_choose_m(int digits) {
    // 文章中给出的近似最优：
    // m ≈ sqrt(N / (2 * log10(2)))
    const double log10_2 = 0.30102999566398119521;
    double m = sqrt((double)digits / (2.0 * log10_2));
    int mi = (int)(m + 0.5);
    if (mi < 3) mi = 3;
    return mi;
}

static int auto_choose_terms(int digits, int m) {
    // 文章里给出的估算：
    // terms ≈ N / (2 m log10(2))
    const double log10_2 = 0.30102999566398119521;
    double t = (double)digits / (2.0 * (double)m * log10_2);
    int ti = (int)ceil(t) + EXTRA_TERMS;
    if (ti < 1) ti = 1;
    return ti;
}

int main(void) {
    mpfr_prec_t prec = (mpfr_prec_t)(DIGITS * 3.3219280948873626) + 256;

    int m = (FORCE_M >= 3) ? FORCE_M : auto_choose_m(DIGITS);
    int terms = (FORCE_TERMS >= 1) ? FORCE_TERMS : auto_choose_terms(DIGITS, m);

    mpfr_t T, tmp, x, x2, term, sum;
    mpfr_t pi_calc, pi_ref, diff, log10diff;

    mpfr_inits2(
        prec,
        T, tmp, x, x2, term, sum,
        pi_calc, pi_ref, diff, log10diff,
        (mpfr_ptr)0
    );

    double t0 = now_seconds();

    // ---------- 1) 构造 T_{m-3} ----------
    // T0 = sqrt(2)
    mpfr_set_ui(T, 2, MPFR_RNDN);
    mpfr_sqrt(T, T, MPFR_RNDN);

    // T_n = sqrt(2 + T_{n-1})
    for (int n = 1; n <= m - 3; ++n) {
        mpfr_add_ui(T, T, 2, MPFR_RNDN);
        mpfr_sqrt(T, T, MPFR_RNDN);
    }

    // ---------- 2) X_m = 1/2 * sqrt(2 - T_{m-3}) ----------
    mpfr_ui_sub(x, 2, T, MPFR_RNDN);   // 2 - T
    mpfr_sqrt(x, x, MPFR_RNDN);        // sqrt(2 - T)
    mpfr_div_2ui(x, x, 1, MPFR_RNDN);  // /2

    // ---------- 3) arcsin(x) 级数展开 ----------
    // arcsin(x) = x + (1/2)x^3/3 + (1*3)/(2*4)x^5/5 + ...
    // 用递推：
    // term_{k+1} = term_k * x^2 * (2k+1)^2 / ((2k+2)(2k+3))
    mpfr_mul(x2, x, x, MPFR_RNDN);

    mpfr_set(term, x, MPFR_RNDN);
    mpfr_set(sum, term, MPFR_RNDN);

    for (int k = 0; k < terms - 1; ++k) {
        unsigned long a = (unsigned long)(2 * k + 1);
        unsigned long b = (unsigned long)(2 * k + 2);
        unsigned long c = (unsigned long)(2 * k + 3);

        mpfr_mul(term, term, x2, MPFR_RNDN);      // * x^2
        mpfr_mul_ui(term, term, a * a, MPFR_RNDN);
        mpfr_div_ui(term, term, b * c, MPFR_RNDN);

        mpfr_add(sum, sum, term, MPFR_RNDN);
    }

    // ---------- 4) pi = 2^m * arcsin(x) ----------
    mpfr_mul_2si(pi_calc, sum, m, MPFR_RNDN);

    double t1 = now_seconds();

    // ---------- 5) 参考值与误差 ----------
    mpfr_const_pi(pi_ref, MPFR_RNDN);

    mpfr_sub(diff, pi_calc, pi_ref, MPFR_RNDN);
    mpfr_abs(diff, diff, MPFR_RNDN);

    printf("算法: nested-sqrt + arcsin-series\n");
    printf("目标小数位数: %d\n", DIGITS);
    printf("m = %d\n", m);
    printf("级数展开项数: %d\n", terms);
    printf("计算耗时: %.6f 秒\n", t1 - t0);

    mpfr_printf("pi_calc (前50位) = %.50RNf\n", pi_calc);
    mpfr_printf("pi_ref  (前50位) = %.50RNf\n", pi_ref);
    mpfr_printf("|pi_calc - pi_ref| = %.20Re\n", diff);

    if (mpfr_zero_p(diff)) {
        printf("在当前精度设置下，误差为 0\n");
    } else {
        mpfr_log10(log10diff, diff, MPFR_RNDN);
        mpfr_printf("log10(误差) = %.10Rf\n", log10diff);
    }

    mpfr_clears(
        T, tmp, x, x2, term, sum,
        pi_calc, pi_ref, diff, log10diff,
        (mpfr_ptr)0
    );
    mpfr_free_cache();

    return 0;
}
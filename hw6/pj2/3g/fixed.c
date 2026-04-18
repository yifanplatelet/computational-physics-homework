#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpfr.h>

#define DIGITS        1000000
#define GUARD_DIGITS  20
#define FORCE_M       -1      // -1: auto; >=3: manual
#define M_SCAN_RADIUS 8

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static int base_m_estimate(int digits) {
    const double log10_2 = 0.30102999566398119521;
    double m = sqrt((double)digits / (2.0 * log10_2));
    int mi = (int)(m + 0.5);
    if (mi < 3) mi = 3;
    return mi;
}

static double estimated_terms_for_m(int digits, int m) {
    const double log10_2 = 0.30102999566398119521;
    return (double)digits / (2.0 * (double)m * log10_2);
}

static int choose_best_m(int digits) {
    int m0 = base_m_estimate(digits);
    int best_m = m0;
    double best_cost = 1e300;

    for (int m = m0 - M_SCAN_RADIUS; m <= m0 + M_SCAN_RADIUS; ++m) {
        if (m < 3) continue;
        double terms = estimated_terms_for_m(digits, m);
        double cost = (double)m + terms;
        if (cost < best_cost) {
            best_cost = cost;
            best_m = m;
        }
    }
    return best_m;
}

int main(void) {
    mpfr_prec_t prec = (mpfr_prec_t)(DIGITS * 3.3219280948873626) + 256;

    int m = (FORCE_M >= 3) ? FORCE_M : choose_best_m(DIGITS);

    mpfr_t T, x, x2, term, sum;
    mpfr_t pi_calc, pi_ref, diff, log10diff;
    mpfr_t threshold, scaled_term;

    mpfr_inits2(
        prec,
        T, x, x2, term, sum,
        pi_calc, pi_ref, diff, log10diff,
        threshold, scaled_term,
        (mpfr_ptr)0
    );

    // threshold = 10^-(DIGITS + GUARD_DIGITS)
    mpfr_set_ui(threshold, 10, MPFR_RNDN);
    mpfr_pow_si(threshold, threshold, -(long)(DIGITS + GUARD_DIGITS), MPFR_RNDN);

    double t0 = now_seconds();

    // 1) 构造 T_{m-3}
    // T0 = sqrt(2)
    mpfr_set_ui(T, 2, MPFR_RNDN);
    mpfr_sqrt(T, T, MPFR_RNDN);

    for (int n = 1; n <= m - 3; ++n) {
        mpfr_add_ui(T, T, 2, MPFR_RNDN);
        mpfr_sqrt(T, T, MPFR_RNDN);
    }

    // 2) x = X_m = 1/2 * sqrt(2 - T_{m-3})
    mpfr_ui_sub(x, 2, T, MPFR_RNDN);
    mpfr_sqrt(x, x, MPFR_RNDN);
    mpfr_div_2ui(x, x, 1, MPFR_RNDN);

    // 3) arcsin(x) 级数
    // arcsin(x)=sum_{k>=0} a_k
    // a_0 = x
    // a_{k+1}=a_k * x^2 * (2k+1)^2 / ((2k+2)(2k+3))
    mpfr_mul(x2, x, x, MPFR_RNDN);

    mpfr_set(term, x, MPFR_RNDN);
    mpfr_set(sum, term, MPFR_RNDN);

    int terms_used = 1;

    for (int k = 0;; ++k) {
        unsigned long a = (unsigned long)(2 * k + 1);
        unsigned long b = (unsigned long)(2 * k + 2);
        unsigned long c = (unsigned long)(2 * k + 3);

        mpfr_mul(term, term, x2, MPFR_RNDN);
        mpfr_mul_ui(term, term, a * a, MPFR_RNDN);
        mpfr_div_ui(term, term, b * c, MPFR_RNDN);

        mpfr_add(sum, sum, term, MPFR_RNDN);
        ++terms_used;

        // 判断最终对 pi 的贡献是否已足够小
        mpfr_mul_2si(scaled_term, term, m, MPFR_RNDN);
        mpfr_abs(scaled_term, scaled_term, MPFR_RNDN);

        if (mpfr_cmp(scaled_term, threshold) < 0) {
            break;
        }

        // 防止异常情况下无限循环
        if (terms_used > 100000000) {
            fprintf(stderr, "Too many terms, aborting.\n");
            break;
        }
    }

    // 4) pi = 2^m * arcsin(x)
    mpfr_mul_2si(pi_calc, sum, m, MPFR_RNDN);

    double t1 = now_seconds();

    // 5) 参考值与误差
    mpfr_const_pi(pi_ref, MPFR_RNDN);

    mpfr_sub(diff, pi_calc, pi_ref, MPFR_RNDN);
    mpfr_abs(diff, diff, MPFR_RNDN);

    printf("算法: nested-sqrt + arcsin-series (optimized)\n");
    printf("目标小数位数: %d\n", DIGITS);
    printf("m = %d\n", m);
    printf("级数展开项数: %d\n", terms_used);
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
        T, x, x2, term, sum,
        pi_calc, pi_ref, diff, log10diff,
        threshold, scaled_term,
        (mpfr_ptr)0
    );
    mpfr_free_cache();

    return 0;
}
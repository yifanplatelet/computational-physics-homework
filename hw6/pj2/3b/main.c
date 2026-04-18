#include <stdio.h>
#include <time.h>
#include <mpfr.h>

#define DIGITS 10000000
#define EXTRA_ITERS 2
#define FORCE_ITERS_NONIC 8   // -1 表示自动估计；>=0 表示手动指定

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static int auto_iters_for_base(int digits, int base, int extra) {
    int iters = 0;
    double covered = 1.0;
    while (covered < (double)digits) {
        covered *= (double)base;
        ++iters;
    }
    return iters + extra;
}

int main(void) {
    mpfr_prec_t prec = (mpfr_prec_t)(DIGITS * 3.3219280948873626) + 256;
    int iters = (FORCE_ITERS_NONIC >= 0)
        ? FORCE_ITERS_NONIC
        : auto_iters_for_base(DIGITS, 9, EXTRA_ITERS);

    mpfr_t a, r, s;
    mpfr_t t, u, v, w;
    mpfr_t s_next, r_next;
    mpfr_t pi_calc, pi_ref, diff, log10diff;

    mpfr_t tmp1, tmp2, tmp3, tmp4, tmp5;
    mpfr_t pow3term;

    mpfr_inits2(
        prec,
        a, r, s, t, u, v, w, s_next, r_next,
        pi_calc, pi_ref, diff, log10diff,
        tmp1, tmp2, tmp3, tmp4, tmp5, pow3term,
        (mpfr_ptr)0
    );

    // a0 = 1/3
    mpfr_set_ui(a, 1, MPFR_RNDN);
    mpfr_div_ui(a, a, 3, MPFR_RNDN);

    // r0 = (sqrt(3) - 1) / 2
    mpfr_set_ui(r, 3, MPFR_RNDN);
    mpfr_sqrt(r, r, MPFR_RNDN);
    mpfr_sub_ui(r, r, 1, MPFR_RNDN);
    mpfr_div_ui(r, r, 2, MPFR_RNDN);

    // s0 = cbrt(1 - r0^3)
    mpfr_mul(s, r, r, MPFR_RNDN);      // r^2
    mpfr_mul(s, s, r, MPFR_RNDN);      // r^3
    mpfr_ui_sub(s, 1, s, MPFR_RNDN);   // 1 - r^3
    mpfr_cbrt(s, s, MPFR_RNDN);        // s0

    double t0_sec = now_seconds();

    for (int n = 0; n < iters; ++n) {
        // t_{n+1} = 1 + 2r_n
        mpfr_mul_ui(t, r, 2, MPFR_RNDN);
        mpfr_add_ui(t, t, 1, MPFR_RNDN);

        // u_{n+1} = cbrt(9 r_n (1 + r_n + r_n^2))
        mpfr_mul(tmp1, r, r, MPFR_RNDN);       // r^2
        mpfr_add(tmp2, r, tmp1, MPFR_RNDN);    // r + r^2
        mpfr_add_ui(tmp2, tmp2, 1, MPFR_RNDN); // 1 + r + r^2
        mpfr_mul(tmp1, r, tmp2, MPFR_RNDN);    // r(1+r+r^2)
        mpfr_mul_ui(tmp1, tmp1, 9, MPFR_RNDN); // 9r(1+r+r^2)
        mpfr_cbrt(u, tmp1, MPFR_RNDN);

        // v_{n+1} = t^2 + t*u + u^2
        mpfr_mul(tmp1, t, t, MPFR_RNDN);       // t^2
        mpfr_mul(tmp2, t, u, MPFR_RNDN);       // t*u
        mpfr_mul(tmp3, u, u, MPFR_RNDN);       // u^2
        mpfr_add(v, tmp1, tmp2, MPFR_RNDN);
        mpfr_add(v, v, tmp3, MPFR_RNDN);

        // w_{n+1} = 27 (1 + s_n + s_n^2) / v_{n+1}
        mpfr_mul(tmp1, s, s, MPFR_RNDN);       // s^2
        mpfr_add(tmp2, s, tmp1, MPFR_RNDN);    // s + s^2
        mpfr_add_ui(tmp2, tmp2, 1, MPFR_RNDN); // 1 + s + s^2
        mpfr_mul_ui(tmp2, tmp2, 27, MPFR_RNDN);
        mpfr_div(w, tmp2, v, MPFR_RNDN);

        // a_{n+1} = w a_n + 3^(2n-1) (1 - w)
        // 这里按公式直接实现；n=0 时指数为 -1，即 1/3
        mpfr_set_ui(pow3term, 1, MPFR_RNDN);
        if (2 * n - 1 >= 0) {
            mpfr_mul_ui(pow3term, pow3term, 3, MPFR_RNDN);
            for (int i = 1; i < 2 * n - 1; ++i) {
                mpfr_mul_ui(pow3term, pow3term, 3, MPFR_RNDN);
            }
        } else {
            mpfr_div_ui(pow3term, pow3term, 3, MPFR_RNDN);
        }

        mpfr_ui_sub(tmp1, 1, w, MPFR_RNDN);    // 1 - w
        mpfr_mul(tmp1, tmp1, pow3term, MPFR_RNDN);
        mpfr_mul(tmp2, w, a, MPFR_RNDN);
        mpfr_add(a, tmp2, tmp1, MPFR_RNDN);

        // s_{n+1} = (1 - r_n)^3 / ((t_{n+1} + 2u_{n+1}) v_{n+1})
        mpfr_ui_sub(tmp1, 1, r, MPFR_RNDN);    // 1 - r
        mpfr_mul(tmp2, tmp1, tmp1, MPFR_RNDN); // (1-r)^2
        mpfr_mul(tmp1, tmp2, tmp1, MPFR_RNDN); // (1-r)^3

        mpfr_mul_ui(tmp2, u, 2, MPFR_RNDN);    // 2u
        mpfr_add(tmp2, t, tmp2, MPFR_RNDN);    // t + 2u
        mpfr_mul(tmp2, tmp2, v, MPFR_RNDN);    // (t+2u)v
        mpfr_div(s_next, tmp1, tmp2, MPFR_RNDN);

        // r_{n+1} = cbrt(1 - s_{n+1}^3)
        mpfr_mul(tmp1, s_next, s_next, MPFR_RNDN);
        mpfr_mul(tmp1, tmp1, s_next, MPFR_RNDN); // s_next^3
        mpfr_ui_sub(tmp1, 1, tmp1, MPFR_RNDN);
        mpfr_cbrt(r_next, tmp1, MPFR_RNDN);

        mpfr_set(s, s_next, MPFR_RNDN);
        mpfr_set(r, r_next, MPFR_RNDN);
    }

    double t1_sec = now_seconds();

    // a_n -> 1/pi，所以 pi ≈ 1 / a_n
    mpfr_ui_div(pi_calc, 1, a, MPFR_RNDN);

    // 参考值
    mpfr_const_pi(pi_ref, MPFR_RNDN);

    // diff = |pi_calc - pi_ref|
    mpfr_sub(diff, pi_calc, pi_ref, MPFR_RNDN);
    mpfr_abs(diff, diff, MPFR_RNDN);

    printf("算法: Borwein nonic (9次收敛)\n");
    printf("目标小数位数: %d\n", DIGITS);
    printf("迭代次数: %d\n", iters);
    printf("计算耗时: %.6f 秒\n", t1_sec - t0_sec);

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
        a, r, s, t, u, v, w, s_next, r_next,
        pi_calc, pi_ref, diff, log10diff,
        tmp1, tmp2, tmp3, tmp4, tmp5, pow3term,
        (mpfr_ptr)0
    );
    mpfr_free_cache();

    return 0;
}
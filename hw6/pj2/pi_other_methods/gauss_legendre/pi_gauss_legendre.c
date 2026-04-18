#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpfr.h>

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static int auto_iters(long digits) {
    long covered = 1;
    int iters = 0;
    while (covered < digits) {
        covered <<= 1;
        ++iters;
    }
    return iters + 3;
}

int main(int argc, char **argv) {
    long digits = 1000000;
    int iters = -1;
    if (argc >= 2) digits = atol(argv[1]);
    if (argc >= 3) iters = atoi(argv[2]);
    if (digits <= 0) {
        fprintf(stderr, "Usage: %s [digits] [iters=-1(auto)]\n", argv[0]);
        return 1;
    }
    if (iters < 0) iters = auto_iters(digits);

    mpfr_prec_t prec = (mpfr_prec_t)(digits * 3.32192809488736234787) + 256;
    mpfr_t a, b, t, p, an, bn, tn, tmp1, tmp2, pi;
    mpfr_inits2(prec, a, b, t, p, an, bn, tn, tmp1, tmp2, pi, (mpfr_ptr)0);

    mpfr_set_ui(a, 1, MPFR_RNDN);
    mpfr_set_ui(b, 2, MPFR_RNDN);
    mpfr_sqrt(b, b, MPFR_RNDN);
    mpfr_ui_div(b, 1, b, MPFR_RNDN);
    mpfr_set_ui(t, 1, MPFR_RNDN);
    mpfr_div_2ui(t, t, 2, MPFR_RNDN);
    mpfr_div_2ui(t, t, 1, MPFR_RNDN);   /* t = 1/4 */
    mpfr_set_ui(p, 1, MPFR_RNDN);

    double t0 = now_seconds();

    for (int i = 0; i < iters; ++i) {
        mpfr_add(an, a, b, MPFR_RNDN);
        mpfr_div_2ui(an, an, 1, MPFR_RNDN);

        mpfr_mul(bn, a, b, MPFR_RNDN);
        mpfr_sqrt(bn, bn, MPFR_RNDN);

        mpfr_sub(tmp1, a, an, MPFR_RNDN);
        mpfr_mul(tmp1, tmp1, tmp1, MPFR_RNDN);
        mpfr_mul(tmp1, tmp1, p, MPFR_RNDN);
        mpfr_sub(tn, t, tmp1, MPFR_RNDN);

        mpfr_mul_2ui(p, p, 1, MPFR_RNDN);
        mpfr_set(a, an, MPFR_RNDN);
        mpfr_set(b, bn, MPFR_RNDN);
        mpfr_set(t, tn, MPFR_RNDN);
    }

    mpfr_add(tmp1, a, b, MPFR_RNDN);
    mpfr_mul(tmp1, tmp1, tmp1, MPFR_RNDN);
    mpfr_mul_ui(tmp2, t, 4, MPFR_RNDN);
    mpfr_div(pi, tmp1, tmp2, MPFR_RNDN);

    double t1 = now_seconds();

    printf("method=gauss_legendre\n");
    printf("digits=%ld\n", digits);
    printf("iterations=%d\n", iters);
    printf("compute_time=%.6f s\n", t1 - t0);
    printf("total_time=%.6f s\n", t1 - t0);
    mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);

    mpfr_clears(a, b, t, p, an, bn, tn, tmp1, tmp2, pi, (mpfr_ptr)0);
    mpfr_free_cache();
    return 0;
}

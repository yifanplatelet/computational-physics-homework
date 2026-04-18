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
    long double covered = 1.0L;
    int iters = 0;
    while (covered < (long double)digits) {
        covered *= 4.0L;
        ++iters;
    }
    return iters + 2;
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
    mpfr_t a, y, yn, t, t2, t3, t4, t5, t6, pi;
    mpfr_inits2(prec, a, y, yn, t, t2, t3, t4, t5, t6, pi, (mpfr_ptr)0);

    mpfr_set_ui(t, 2, MPFR_RNDN);
    mpfr_sqrt(t, t, MPFR_RNDN);          /* sqrt(2) */
    mpfr_sub_ui(y, t, 1, MPFR_RNDN);     /* y0 = sqrt(2)-1 */
    mpfr_mul_ui(a, t, 4, MPFR_RNDN);
    mpfr_ui_sub(a, 6, a, MPFR_RNDN);     /* a0 = 6-4sqrt(2) */

    double t0s = now_seconds();

    for (int n = 0; n < iters; ++n) {
        mpfr_mul(t, y, y, MPFR_RNDN);
        mpfr_mul(t, t, t, MPFR_RNDN);    /* y^4 */
        mpfr_ui_sub(t, 1, t, MPFR_RNDN);
        mpfr_sqrt(t, t, MPFR_RNDN);
        mpfr_sqrt(t, t, MPFR_RNDN);      /* (1-y^4)^(1/4) */

        mpfr_ui_sub(t2, 1, t, MPFR_RNDN);
        mpfr_add_ui(t3, t, 1, MPFR_RNDN);
        mpfr_div(yn, t2, t3, MPFR_RNDN);

        mpfr_add_ui(t2, yn, 1, MPFR_RNDN);
        mpfr_mul(t3, t2, t2, MPFR_RNDN);
        mpfr_mul(t3, t3, t3, MPFR_RNDN); /* (1+y)^4 */

        mpfr_mul(t4, a, t3, MPFR_RNDN);

        mpfr_mul(t5, yn, yn, MPFR_RNDN);
        mpfr_add(t6, yn, t5, MPFR_RNDN);
        mpfr_add_ui(t6, t6, 1, MPFR_RNDN);
        mpfr_mul(t5, yn, t6, MPFR_RNDN);
        mpfr_mul_2si(t5, t5, 2 * n + 3, MPFR_RNDN);

        mpfr_sub(a, t4, t5, MPFR_RNDN);
        mpfr_set(y, yn, MPFR_RNDN);
    }

    mpfr_ui_div(pi, 1, a, MPFR_RNDN);
    double t1s = now_seconds();

    printf("method=borwein_quartic\n");
    printf("digits=%ld\n", digits);
    printf("iterations=%d\n", iters);
    printf("compute_time=%.6f s\n", t1s - t0s);
    printf("total_time=%.6f s\n", t1s - t0s);
    mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);

    mpfr_clears(a, y, yn, t, t2, t3, t4, t5, t6, pi, (mpfr_ptr)0);
    mpfr_free_cache();
    return 0;
}

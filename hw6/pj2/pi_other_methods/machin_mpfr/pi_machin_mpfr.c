#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpfr.h>

#define GUARD_DIGITS 20

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static long atan_series(mpfr_t out, unsigned long q, long digits, mpfr_prec_t prec) {
    mpfr_t x, x2, term, sum, threshold, tmp;
    mpfr_inits2(prec, x, x2, term, sum, threshold, tmp, (mpfr_ptr)0);

    mpfr_set_ui(x, 1, MPFR_RNDN);
    mpfr_div_ui(x, x, q, MPFR_RNDN);
    mpfr_mul(x2, x, x, MPFR_RNDN);
    mpfr_set(term, x, MPFR_RNDN);
    mpfr_set(sum, term, MPFR_RNDN);

    mpfr_set_ui(threshold, 10, MPFR_RNDN);
    mpfr_pow_si(threshold, threshold, -(long)(digits + GUARD_DIGITS), MPFR_RNDN);

    long used = 1;
    for (long k = 0;; ++k) {
        mpfr_mul(term, term, x2, MPFR_RNDN);
        mpfr_neg(term, term, MPFR_RNDN);
        mpfr_mul_si(term, term, 2 * k + 1, MPFR_RNDN);
        mpfr_div_si(term, term, 2 * k + 3, MPFR_RNDN);
        mpfr_add(sum, sum, term, MPFR_RNDN);
        ++used;

        mpfr_abs(tmp, term, MPFR_RNDN);
        if (mpfr_cmp(tmp, threshold) < 0) break;
        if (used > 100000000L) break;
    }

    mpfr_set(out, sum, MPFR_RNDN);
    mpfr_clears(x, x2, term, sum, threshold, tmp, (mpfr_ptr)0);
    return used;
}

int main(int argc, char **argv) {
    long digits = 1000000;
    if (argc >= 2) digits = atol(argv[1]);
    if (digits <= 0) {
        fprintf(stderr, "Usage: %s [digits]\n", argv[0]);
        return 1;
    }

    mpfr_prec_t prec = (mpfr_prec_t)(digits * 3.32192809488736234787) + 256;
    mpfr_t a, b, pi;
    mpfr_inits2(prec, a, b, pi, (mpfr_ptr)0);

    double t0 = now_seconds();
    long ta = atan_series(a, 5, digits, prec);
    long tb = atan_series(b, 239, digits, prec);
    mpfr_mul_ui(a, a, 16, MPFR_RNDN);
    mpfr_mul_ui(b, b, 4, MPFR_RNDN);
    mpfr_sub(pi, a, b, MPFR_RNDN);
    double t1 = now_seconds();

    printf("method=machin_mpfr\n");
    printf("digits=%ld\n", digits);
    printf("atan_terms_q5=%ld\n", ta);
    printf("atan_terms_q239=%ld\n", tb);
    printf("compute_time=%.6f s\n", t1 - t0);
    printf("total_time=%.6f s\n", t1 - t0);
    mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);

    mpfr_clears(a, b, pi, (mpfr_ptr)0);
    mpfr_free_cache();
    return 0;
}

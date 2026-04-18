#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>

#define A 13591409L
#define B 545140134L
#define C3_OVER_24 10939058860032000ULL

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static void bs_exact(mpz_t P, mpz_t Q, mpz_t T, long a, long b) {
    if (b - a == 1) {
        if (a == 0) {
            mpz_set_ui(P, 1);
            mpz_set_ui(Q, 1);
            mpz_set_ui(T, A);
        } else {
            mpz_set_si(P, 6 * a - 5);
            mpz_mul_si(P, P, 2 * a - 1);
            mpz_mul_si(P, P, 6 * a - 1);

            mpz_set_si(Q, a);
            mpz_mul(Q, Q, Q);
            mpz_mul_ui(Q, Q, C3_OVER_24);

            mpz_set_si(T, B);
            mpz_mul_si(T, T, a);
            mpz_add_ui(T, T, A);
            mpz_mul(T, T, P);
            if (a & 1L) {
                mpz_neg(T, T);
            }
        }
        return;
    }

    long m = (a + b) / 2;
    mpz_t P1, Q1, T1, P2, Q2, T2, tmp;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, tmp, NULL);

    bs_exact(P1, Q1, T1, a, m);
    bs_exact(P2, Q2, T2, m, b);

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(tmp, T1, Q2);
    mpz_mul(T, P1, T2);
    mpz_add(T, T, tmp);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, tmp, NULL);
}

static void reciprocal_newton(mpfr_t out, const mpfr_t x, mpfr_prec_t target_prec) {
    mpfr_prec_t cur = target_prec < 256 ? target_prec : 256;
    if (cur < 64) cur = 64;

    mpfr_t y, t;
    mpfr_inits2(cur, y, t, (mpfr_ptr)0);

    /* Low-precision seed */
    mpfr_ui_div(y, 1, x, MPFR_RNDN);

    while (cur < target_prec) {
        mpfr_prec_t next = cur * 2;
        if (next > target_prec) next = target_prec;
        mpfr_prec_round(y, next, MPFR_RNDN);
        mpfr_prec_round(t, next, MPFR_RNDN);

        /* y = y * (2 - x*y) */
        mpfr_mul(t, x, y, MPFR_RNDN);
        mpfr_ui_sub(t, 2, t, MPFR_RNDN);
        mpfr_mul(y, y, t, MPFR_RNDN);
        cur = next;
    }

    mpfr_set(out, y, MPFR_RNDN);
    mpfr_clears(y, t, (mpfr_ptr)0);
}

int main(int argc, char **argv) {
    long digits = 1000000;
    if (argc >= 2) digits = atol(argv[1]);
    if (digits <= 0) {
        fprintf(stderr, "Usage: %s [digits]\n", argv[0]);
        return 1;
    }

    long terms = digits / 14 + 10;
    mpfr_prec_t prec = (mpfr_prec_t)(digits * 3.32192809488736234787) + 256;

    mpz_t P, Q, T;
    mpz_inits(P, Q, T, NULL);

    double t0 = now_seconds();
    bs_exact(P, Q, T, 0, terms);
    double t1 = now_seconds();

    mpfr_t qf, tf, inv_tf, sqrt10005, factor, pi;
    mpfr_inits2(prec, qf, tf, inv_tf, sqrt10005, factor, pi, (mpfr_ptr)0);
    mpfr_set_z(qf, Q, MPFR_RNDN);
    mpfr_set_z(tf, T, MPFR_RNDN);

    reciprocal_newton(inv_tf, tf, prec);

    mpfr_set_ui(sqrt10005, 10005, MPFR_RNDN);
    mpfr_sqrt(sqrt10005, sqrt10005, MPFR_RNDN);
    mpfr_mul_ui(factor, sqrt10005, 426880, MPFR_RNDN);

    mpfr_mul(pi, factor, qf, MPFR_RNDN);
    mpfr_mul(pi, pi, inv_tf, MPFR_RNDN);

    double t2 = now_seconds();

    printf("method=newton_endgame_bs\n");
    printf("digits=%ld\n", digits);
    printf("terms=%ld\n", terms);
    printf("binary_splitting_time=%.6f s\n", t1 - t0);
    printf("final_eval_time=%.6f s\n", t2 - t1);
    printf("total_time=%.6f s\n", t2 - t0);
    mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);

    mpz_clears(P, Q, T, NULL);
    mpfr_clears(qf, tf, inv_tf, sqrt10005, factor, pi, (mpfr_ptr)0);
    mpfr_free_cache();
    return 0;
}

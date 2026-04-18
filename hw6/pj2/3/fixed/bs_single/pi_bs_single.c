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
    return (double) ts.tv_sec + (double) ts.tv_nsec / 1e9;
}

static void bs(mpz_t P, mpz_t Q, mpz_t T, long a, long b) {
    if (b - a == 1) {
        if (a == 0) {
            mpz_set_ui(P, 1);
            mpz_set_ui(Q, 1);
            mpz_set_ui(T, A);
        } else {
            mpz_t t1, t2;
            mpz_inits(t1, t2, NULL);

            mpz_set_si(P, 6 * a - 5);
            mpz_mul_si(P, P, 2 * a - 1);
            mpz_mul_si(P, P, 6 * a - 1);

            mpz_set_si(Q, a);
            mpz_mul(Q, Q, Q);
            mpz_mul_ui(Q, Q, C3_OVER_24);

            mpz_set_si(t1, B);
            mpz_mul_si(t1, t1, a);
            mpz_add_ui(t1, t1, A);

            mpz_mul(T, P, t1);
            if (a & 1L) {
                mpz_neg(T, T);
            }

            mpz_clears(t1, t2, NULL);
        }
        return;
    }

    long m = (a + b) / 2;

    mpz_t P1, Q1, T1, P2, Q2, T2, tmp;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, tmp, NULL);

    bs(P1, Q1, T1, a, m);
    bs(P2, Q2, T2, m, b);

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);

    mpz_mul(tmp, T1, Q2);
    mpz_mul(T, P1, T2);
    mpz_add(T, T, tmp);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, tmp, NULL);
}

int main(int argc, char **argv) {
    long digits = 1000000;
    if (argc >= 2) {
        digits = atol(argv[1]);
        if (digits <= 0) {
            fprintf(stderr, "Usage: %s [digits]\n", argv[0]);
            return 1;
        }
    }

    long terms = digits / 14 + 10;
    mpfr_prec_t prec = (mpfr_prec_t) (digits * 3.32192809488736234787) + 256;

    mpz_t P, Q, T;
    mpz_inits(P, Q, T, NULL);

    double t0 = now_seconds();
    bs(P, Q, T, 0, terms);
    double t1 = now_seconds();

    mpfr_t pi, sqrt10005, factor, qf, tf;
    mpfr_inits2(prec, pi, sqrt10005, factor, qf, tf, (mpfr_ptr) 0);

    mpfr_set_ui(sqrt10005, 10005, MPFR_RNDN);
    mpfr_sqrt(sqrt10005, sqrt10005, MPFR_RNDN);
    mpfr_mul_ui(factor, sqrt10005, 426880, MPFR_RNDN);

    mpfr_set_z(qf, Q, MPFR_RNDN);
    mpfr_set_z(tf, T, MPFR_RNDN);
    mpfr_mul(pi, factor, qf, MPFR_RNDN);
    mpfr_div(pi, pi, tf, MPFR_RNDN);

    double t2 = now_seconds();

    printf("digits=%ld\n", digits);
    printf("terms=%ld\n", terms);
    printf("binary_splitting_time=%.6f s\n", t1 - t0);
    printf("final_eval_time=%.6f s\n", t2 - t1);
    printf("total_time=%.6f s\n", t2 - t0);
    mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);

    mpz_clears(P, Q, T, NULL);
    mpfr_clears(pi, sqrt10005, factor, qf, tf, (mpfr_ptr) 0);
    mpfr_free_cache();
    return 0;
}

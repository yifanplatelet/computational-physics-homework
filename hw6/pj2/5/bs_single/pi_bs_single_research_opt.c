#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>

#define A 13591409UL
#define B 545140134UL
#define C3_OVER_24 10939058860032000UL

static long g_leaf_size = 64;

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static inline void term_pqt(mpz_t p, mpz_t q, mpz_t t, unsigned long k) {
    mpz_set_ui(p, 6 * k - 5);
    mpz_mul_ui(p, p, 2 * k - 1);
    mpz_mul_ui(p, p, 6 * k - 1);

    mpz_set_ui(q, k);
    mpz_mul(q, q, q);
    mpz_mul_ui(q, q, k);
    mpz_mul_ui(q, q, C3_OVER_24);

    mpz_set_ui(t, B);
    mpz_mul_ui(t, t, k);
    mpz_add_ui(t, t, A);
    mpz_mul(t, t, p);
    if (k & 1UL) {
        mpz_neg(t, t);
    }
}

// Backward block summation over [a, b). This reduces recursion overhead and
// matches the "sum low terms first" intuition often used in optimized BS code.
static void bs_leaf(mpz_t P, mpz_t Q, mpz_t T, long a, long b) {
    mpz_t p, q, t, tmp;
    mpz_inits(p, q, t, tmp, NULL);

    mpz_set_ui(P, 1);
    mpz_set_ui(Q, 1);
    mpz_set_ui(T, 0);

    for (long k = b - 1; k >= a; --k) {
        if (k == 0) {
            // Prepend the k = 0 term: P and Q stay unchanged, T <- A*Q + T.
            mpz_addmul_ui(T, Q, A);
            continue;
        }

        term_pqt(p, q, t, (unsigned long)k);

        // Prepend term k to the already accumulated suffix [k+1, b):
        // P <- p * P
        // Q <- q * Q
        // T <- t * Q + p * T
        mpz_mul(tmp, p, T);
        mpz_mul(T, t, Q);
        mpz_add(T, T, tmp);

        mpz_mul(P, P, p);
        mpz_mul(Q, Q, q);
    }

    mpz_clears(p, q, t, tmp, NULL);
}

static void bs(mpz_t P, mpz_t Q, mpz_t T, long a, long b) {
    long n = b - a;
    if (n <= g_leaf_size) {
        bs_leaf(P, Q, T, a, b);
        return;
    }

    long m = a + n / 2;

    mpz_t P1, Q1, T1, P2, Q2, T2;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, NULL);

    bs(P1, Q1, T1, a, m);
    bs(P2, Q2, T2, m, b);

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);

    // T = P1*T2 + T1*Q2
    mpz_mul(T, P1, T2);
    mpz_addmul(T, T1, Q2);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, NULL);
}

int main(int argc, char **argv) {
    long digits = 1000000;
    if (argc >= 2) {
        digits = atol(argv[1]);
    }
    if (argc >= 3) {
        g_leaf_size = atol(argv[2]);
    }
    if (digits <= 0 || g_leaf_size <= 0) {
        fprintf(stderr, "Usage: %s [digits] [leaf_size]\n", argv[0]);
        return 1;
    }

    long terms = digits / 14 + 10;
    mpfr_prec_t prec = (mpfr_prec_t)(digits * 3.32192809488736234787) + 256;

    mpz_t P, Q, T;
    mpz_inits(P, Q, T, NULL);

    double t0 = now_seconds();
    bs(P, Q, T, 0, terms);
    double t1 = now_seconds();

    mpfr_t pi, sqrt10005, factor, qf, tf;
    mpfr_inits2(prec, pi, sqrt10005, factor, qf, tf, (mpfr_ptr)0);

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
    printf("leaf_size=%ld\n", g_leaf_size);
    printf("binary_splitting_time=%.6f s\n", t1 - t0);
    printf("final_eval_time=%.6f s\n", t2 - t1);
    printf("total_time=%.6f s\n", t2 - t0);
    mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);

    mpz_clears(P, Q, T, NULL);
    mpfr_clears(pi, sqrt10005, factor, qf, tf, (mpfr_ptr)0);
    mpfr_free_cache();
    return 0;
}

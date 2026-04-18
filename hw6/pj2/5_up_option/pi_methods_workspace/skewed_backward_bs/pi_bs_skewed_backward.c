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

static void leaf_term(mpz_t P, mpz_t Q, mpz_t T, long k) {
    if (k == 0) {
        mpz_set_ui(P, 1);
        mpz_set_ui(Q, 1);
        mpz_set_ui(T, A);
        return;
    }

    mpz_set_si(P, 6 * k - 5);
    mpz_mul_si(P, P, 2 * k - 1);
    mpz_mul_si(P, P, 6 * k - 1);

    mpz_set_si(Q, k);
    mpz_mul(Q, Q, Q);
    mpz_mul_ui(Q, Q, C3_OVER_24);

    mpz_set_si(T, B);
    mpz_mul_si(T, T, k);
    mpz_add_ui(T, T, A);
    mpz_mul(T, T, P);
    if (k & 1L) {
        mpz_neg(T, T);
    }
}

static void block_backward(mpz_t P, mpz_t Q, mpz_t T, long a, long b) {
    mpz_set_ui(P, 1);
    mpz_set_ui(Q, 1);
    mpz_set_ui(T, 0);

    mpz_t pk, qk, tk, tmp;
    mpz_inits(pk, qk, tk, tmp, NULL);

    for (long k = b - 1; k >= a; --k) {
        leaf_term(pk, qk, tk, k);

        mpz_mul(tmp, tk, Q);      /* tk * Q_acc */
        mpz_mul(T, pk, T);        /* pk * T_acc */
        mpz_add(T, T, tmp);       /* tk*Q + pk*T */

        mpz_mul(P, pk, P);
        mpz_mul(Q, qk, Q);
    }

    mpz_clears(pk, qk, tk, tmp, NULL);
}

static long choose_skew_mid(long a, long b, long skew_num, long skew_den) {
    long len = b - a;
    long m = a + (len * skew_num) / skew_den;
    if (m <= a) m = a + 1;
    if (m >= b) m = b - 1;
    return m;
}

static void bs_skew(mpz_t P, mpz_t Q, mpz_t T, long a, long b,
                    long leaf_size, long skew_num, long skew_den,
                    unsigned long *block_leaves) {
    if (b - a <= leaf_size) {
        block_backward(P, Q, T, a, b);
        if (block_leaves) {
            (*block_leaves)++;
        }
        return;
    }

    long m = choose_skew_mid(a, b, skew_num, skew_den);

    mpz_t P1, Q1, T1, P2, Q2, T2, tmp;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, tmp, NULL);

    bs_skew(P1, Q1, T1, a, m, leaf_size, skew_num, skew_den, block_leaves);
    bs_skew(P2, Q2, T2, m, b, leaf_size, skew_num, skew_den, block_leaves);

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(tmp, T1, Q2);
    mpz_mul(T, P1, T2);
    mpz_add(T, T, tmp);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, tmp, NULL);
}

int main(int argc, char **argv) {
    long digits = 1000000;
    long leaf_size = 32;
    long skew_num = 2;
    long skew_den = 3;

    if (argc >= 2) digits = atol(argv[1]);
    if (argc >= 3) leaf_size = atol(argv[2]);
    if (argc >= 4) skew_num = atol(argv[3]);
    if (argc >= 5) skew_den = atol(argv[4]);

    if (digits <= 0 || leaf_size <= 0 || skew_num <= 0 || skew_den <= 0) {
        fprintf(stderr, "Usage: %s [digits] [leaf_size] [skew_num] [skew_den]\n", argv[0]);
        return 1;
    }

    long terms = digits / 14 + 10;
    mpfr_prec_t prec = (mpfr_prec_t)(digits * 3.32192809488736234787) + 256;

    mpz_t P, Q, T;
    mpz_inits(P, Q, T, NULL);

    unsigned long block_leaves = 0;
    double t0 = now_seconds();
    bs_skew(P, Q, T, 0, terms, leaf_size, skew_num, skew_den, &block_leaves);
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

    printf("method=skewed_backward_bs\n");
    printf("digits=%ld\n", digits);
    printf("terms=%ld\n", terms);
    printf("leaf_size=%ld\n", leaf_size);
    printf("skew_ratio=%ld/%ld\n", skew_num, skew_den);
    printf("block_leaves=%lu\n", block_leaves);
    printf("binary_splitting_time=%.6f s\n", t1 - t0);
    printf("final_eval_time=%.6f s\n", t2 - t1);
    printf("total_time=%.6f s\n", t2 - t0);
    mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);

    mpz_clears(P, Q, T, NULL);
    mpfr_clears(pi, sqrt10005, factor, qf, tf, (mpfr_ptr)0);
    mpfr_free_cache();
    return 0;
}

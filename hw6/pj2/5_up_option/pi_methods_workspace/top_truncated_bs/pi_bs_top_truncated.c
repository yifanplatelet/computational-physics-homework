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

static void exact_to_ratio_sum(mpfr_t r, mpfr_t s, long a, long b, mpfr_prec_t prec, unsigned long *exact_subtrees) {
    mpz_t P, Q, T;
    mpz_inits(P, Q, T, NULL);
    bs_exact(P, Q, T, a, b);
    mpfr_t pf, qf, tf;
    mpfr_inits2(prec, pf, qf, tf, (mpfr_ptr)0);
    mpfr_set_z(pf, P, MPFR_RNDN);
    mpfr_set_z(qf, Q, MPFR_RNDN);
    mpfr_set_z(tf, T, MPFR_RNDN);
    mpfr_div(r, pf, qf, MPFR_RNDN);  /* r = P / Q */
    mpfr_div(s, tf, qf, MPFR_RNDN);  /* s = T / Q */
    if (exact_subtrees) {
        (*exact_subtrees)++;
    }
    mpfr_clears(pf, qf, tf, (mpfr_ptr)0);
    mpz_clears(P, Q, T, NULL);
}

static void bs_hybrid(mpfr_t r, mpfr_t s, long a, long b, int depth, int approx_depth,
                      long exact_leaf, mpfr_prec_t prec,
                      unsigned long *exact_subtrees, unsigned long *approx_merges) {
    if (b - a <= exact_leaf || depth >= approx_depth) {
        exact_to_ratio_sum(r, s, a, b, prec, exact_subtrees);
        return;
    }

    long m = (a + b) / 2;
    mpfr_t r1, s1, r2, s2, tmp;
    mpfr_inits2(prec, r1, s1, r2, s2, tmp, (mpfr_ptr)0);

    bs_hybrid(r1, s1, a, m, depth + 1, approx_depth, exact_leaf, prec, exact_subtrees, approx_merges);
    bs_hybrid(r2, s2, m, b, depth + 1, approx_depth, exact_leaf, prec, exact_subtrees, approx_merges);

    /* s = s1 + r1 * s2 */
    mpfr_mul(tmp, r1, s2, MPFR_RNDN);
    mpfr_add(s, s1, tmp, MPFR_RNDN);

    /* r = r1 * r2 */
    mpfr_mul(r, r1, r2, MPFR_RNDN);

    if (approx_merges) {
        (*approx_merges)++;
    }

    mpfr_clears(r1, s1, r2, s2, tmp, (mpfr_ptr)0);
}

int main(int argc, char **argv) {
    long digits = 1000000;
    int approx_depth = 3;
    long exact_leaf = 64;

    if (argc >= 2) digits = atol(argv[1]);
    if (argc >= 3) approx_depth = atoi(argv[2]);
    if (argc >= 4) exact_leaf = atol(argv[3]);

    if (digits <= 0 || approx_depth < 0 || exact_leaf <= 0) {
        fprintf(stderr, "Usage: %s [digits] [approx_depth] [exact_leaf]\n", argv[0]);
        return 1;
    }

    long terms = digits / 14 + 10;
    mpfr_prec_t prec = (mpfr_prec_t)(digits * 3.32192809488736234787) + 256;

    mpfr_t r, s;
    mpfr_inits2(prec, r, s, (mpfr_ptr)0);

    unsigned long exact_subtrees = 0;
    unsigned long approx_merges = 0;

    double t0 = now_seconds();
    bs_hybrid(r, s, 0, terms, 0, approx_depth, exact_leaf, prec, &exact_subtrees, &approx_merges);
    double t1 = now_seconds();

    mpfr_t pi, sqrt10005, factor;
    mpfr_inits2(prec, pi, sqrt10005, factor, (mpfr_ptr)0);
    mpfr_set_ui(sqrt10005, 10005, MPFR_RNDN);
    mpfr_sqrt(sqrt10005, sqrt10005, MPFR_RNDN);
    mpfr_mul_ui(factor, sqrt10005, 426880, MPFR_RNDN);
    mpfr_div(pi, factor, s, MPFR_RNDN);
    double t2 = now_seconds();

    printf("method=top_truncated_hybrid_bs\n");
    printf("digits=%ld\n", digits);
    printf("terms=%ld\n", terms);
    printf("approx_depth=%d\n", approx_depth);
    printf("exact_leaf=%ld\n", exact_leaf);
    printf("exact_subtrees=%lu\n", exact_subtrees);
    printf("approx_merges=%lu\n", approx_merges);
    printf("hybrid_build_time=%.6f s\n", t1 - t0);
    printf("final_eval_time=%.6f s\n", t2 - t1);
    printf("total_time=%.6f s\n", t2 - t0);
    mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);

    mpfr_clears(r, s, pi, sqrt10005, factor, (mpfr_ptr)0);
    mpfr_free_cache();
    return 0;
}

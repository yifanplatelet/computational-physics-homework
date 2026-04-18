#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define A 13591409L
#define B 545140134L
#define C3_OVER_24 10939058860032000ULL
#define DEFAULT_LEAF_SIZE 64L

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double) ts.tv_sec + (double) ts.tv_nsec / 1e9;
}

static void bs_leaf(mpz_t P, mpz_t Q, mpz_t T, long a, long b) {
    mpz_set_ui(P, 1);
    mpz_set_ui(Q, 1);
    mpz_set_ui(T, 0);

    mpz_t p, q, t, tmp, tmp2;
    mpz_inits(p, q, t, tmp, tmp2, NULL);

    for (long k = a; k < b; ++k) {
        if (k == 0) {
            mpz_set_ui(p, 1);
            mpz_set_ui(q, 1);
            mpz_set_ui(t, A);
        } else {
            mpz_set_si(p, 6 * k - 5);
            mpz_mul_si(p, p, 2 * k - 1);
            mpz_mul_si(p, p, 6 * k - 1);

            mpz_set_si(q, k);
            mpz_mul(q, q, q);
            mpz_mul_ui(q, q, C3_OVER_24);

            mpz_set_si(t, B);
            mpz_mul_si(t, t, k);
            mpz_add_ui(t, t, A);
            mpz_mul(t, t, p);
            if (k & 1L) {
                mpz_neg(t, t);
            }
        }

        mpz_mul(tmp, T, q);
        mpz_mul(tmp2, P, t);
        mpz_add(T, tmp, tmp2);
        mpz_mul(P, P, p);
        mpz_mul(Q, Q, q);
    }

    mpz_clears(p, q, t, tmp, tmp2, NULL);
}

static void bs_serial(mpz_t P, mpz_t Q, mpz_t T, long a, long b, long leaf_size) {
    long len = b - a;
    if (len <= leaf_size) {
        bs_leaf(P, Q, T, a, b);
        return;
    }

    long m = a + len / 2;
    mpz_t P1, Q1, T1, P2, Q2, T2, tmp;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, tmp, NULL);

    bs_serial(P1, Q1, T1, a, m, leaf_size);
    bs_serial(P2, Q2, T2, m, b, leaf_size);

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(tmp, T1, Q2);
    mpz_mul(T, P1, T2);
    mpz_add(T, T, tmp);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, tmp, NULL);
}

static void bs_parallel(mpz_t P, mpz_t Q, mpz_t T, long a, long b, int depth, int max_depth, long leaf_size) {
    long len = b - a;
    if (len <= leaf_size) {
        bs_leaf(P, Q, T, a, b);
        return;
    }
    if (depth >= max_depth) {
        bs_serial(P, Q, T, a, b, leaf_size);
        return;
    }

    long m = a + len / 2;
    mpz_t P1, Q1, T1, P2, Q2, T2, tmp;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, tmp, NULL);

    #pragma omp task shared(P1, Q1, T1)
    bs_parallel(P1, Q1, T1, a, m, depth + 1, max_depth, leaf_size);

    #pragma omp task shared(P2, Q2, T2)
    bs_parallel(P2, Q2, T2, m, b, depth + 1, max_depth, leaf_size);

    #pragma omp taskwait

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(tmp, T1, Q2);
    mpz_mul(T, P1, T2);
    mpz_add(T, T, tmp);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, tmp, NULL);
}

int main(int argc, char **argv) {
    long digits = 1000000;
    int max_depth = 4;
    long leaf_size = DEFAULT_LEAF_SIZE;

    if (argc >= 2) digits = atol(argv[1]);
    if (argc >= 3) max_depth = atoi(argv[2]);
    if (argc >= 4) leaf_size = atol(argv[3]);
    if (digits <= 0 || max_depth < 0 || leaf_size <= 0) {
        fprintf(stderr, "Usage: %s [digits] [task_depth] [leaf_size]\n", argv[0]);
        return 1;
    }

    long terms = digits / 14 + 10;
    mpfr_prec_t prec = (mpfr_prec_t) (digits * 3.32192809488736234787) + 256;

    mpz_t P, Q, T;
    mpz_inits(P, Q, T, NULL);

    double t0 = now_seconds();
    #pragma omp parallel
    {
        #pragma omp single nowait
        bs_parallel(P, Q, T, 0, terms, 0, max_depth, leaf_size);
    }
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
    printf("task_depth=%d\n", max_depth);
    printf("leaf_size=%ld\n", leaf_size);
#ifdef _OPENMP
    printf("omp_threads=%d\n", omp_get_max_threads());
#else
    printf("omp_threads=1 (OpenMP disabled at compile time)\n");
#endif
    printf("binary_splitting_time=%.6f s\n", t1 - t0);
    printf("final_eval_time=%.6f s\n", t2 - t1);
    printf("total_time=%.6f s\n", t2 - t0);
    mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);

    mpz_clears(P, Q, T, NULL);
    mpfr_clears(pi, sqrt10005, factor, qf, tf, (mpfr_ptr) 0);
    mpfr_free_cache();
    return 0;
}

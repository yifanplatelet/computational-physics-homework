#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define A0 1103UL
#define B0 26390UL
#define CONST_3964 24591257856UL

static long g_leaf_size = 64;
static long g_task_min_size = 2048;

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static inline void term_pqt(mpz_t p, mpz_t q, mpz_t t, unsigned long k) {
    mpz_set_ui(p, 4 * k - 3);
    mpz_mul_ui(p, p, 4 * k - 2);
    mpz_mul_ui(p, p, 4 * k - 1);
    mpz_mul_ui(p, p, 4 * k);

    mpz_set_ui(q, k);
    mpz_mul(q, q, q);
    mpz_mul(q, q, q);
    mpz_mul_ui(q, q, CONST_3964);

    mpz_set_ui(t, B0);
    mpz_mul_ui(t, t, k);
    mpz_add_ui(t, t, A0);
    mpz_mul(t, t, p);
}

static void bs_leaf(mpz_t P, mpz_t Q, mpz_t T, long a, long b) {
    mpz_t p, q, t, tmp;
    mpz_inits(p, q, t, tmp, NULL);

    mpz_set_ui(P, 1);
    mpz_set_ui(Q, 1);
    mpz_set_ui(T, 0);

    for (long k = b - 1; k >= a; --k) {
        if (k == 0) {
            mpz_addmul_ui(T, Q, A0);
            continue;
        }
        term_pqt(p, q, t, (unsigned long)k);

        mpz_mul(tmp, p, T);
        mpz_mul(T, t, Q);
        mpz_add(T, T, tmp);

        mpz_mul(P, P, p);
        mpz_mul(Q, Q, q);
    }

    mpz_clears(p, q, t, tmp, NULL);
}

static void bs_serial(mpz_t P, mpz_t Q, mpz_t T, long a, long b) {
    long n = b - a;
    if (n <= g_leaf_size) {
        bs_leaf(P, Q, T, a, b);
        return;
    }
    long m = a + n / 2;
    mpz_t P1, Q1, T1, P2, Q2, T2;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, NULL);

    bs_serial(P1, Q1, T1, a, m);
    bs_serial(P2, Q2, T2, m, b);

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(T, P1, T2);
    mpz_addmul(T, T1, Q2);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, NULL);
}

static void bs_parallel(mpz_t P, mpz_t Q, mpz_t T, long a, long b, int depth, int max_depth) {
    long n = b - a;
    if (n <= g_leaf_size) {
        bs_leaf(P, Q, T, a, b);
        return;
    }
    if (depth >= max_depth || n < g_task_min_size) {
        bs_serial(P, Q, T, a, b);
        return;
    }

    long m = a + n / 2;
    mpz_t P1, Q1, T1, P2, Q2, T2;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, NULL);

    #pragma omp task shared(P1, Q1, T1)
    bs_parallel(P1, Q1, T1, a, m, depth + 1, max_depth);

    #pragma omp task shared(P2, Q2, T2)
    bs_parallel(P2, Q2, T2, m, b, depth + 1, max_depth);

    #pragma omp taskwait

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(T, P1, T2);
    mpz_addmul(T, T1, Q2);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, NULL);
}

int main(int argc, char **argv) {
    long digits = 1000000;
    int max_depth = 4;
    if (argc >= 2) digits = atol(argv[1]);
    if (argc >= 3) max_depth = atoi(argv[2]);
    if (argc >= 4) g_leaf_size = atol(argv[3]);
    if (argc >= 5) g_task_min_size = atol(argv[4]);
    if (digits <= 0 || max_depth < 0 || g_leaf_size <= 0 || g_task_min_size <= 0) {
        fprintf(stderr, "Usage: %s [digits] [task_depth] [leaf_size] [task_min_size]\n", argv[0]);
        return 1;
    }

    /* Ramanujan gives about 8 digits per term. */
    long terms = digits / 8 + 10;
    mpfr_prec_t prec = (mpfr_prec_t)(digits * 3.32192809488736234787) + 256;

    mpz_t P, Q, T;
    mpz_inits(P, Q, T, NULL);

    double t0 = now_seconds();

    #pragma omp parallel
    {
        #pragma omp single nowait
        bs_parallel(P, Q, T, 0, terms, 0, max_depth);
    }

    double t1 = now_seconds();

    mpfr_t pi, factor, qf, tf, sqrt2;
    mpfr_inits2(prec, pi, factor, qf, tf, sqrt2, (mpfr_ptr)0);
    mpfr_set_ui(sqrt2, 2, MPFR_RNDN);
    mpfr_sqrt(sqrt2, sqrt2, MPFR_RNDN);
    mpfr_mul_ui(sqrt2, sqrt2, 2, MPFR_RNDN);      /* 2*sqrt(2) */
    mpfr_ui_div(factor, 9801, sqrt2, MPFR_RNDN);  /* 9801 / (2*sqrt(2)) */
    mpfr_set_z(qf, Q, MPFR_RNDN);
    mpfr_set_z(tf, T, MPFR_RNDN);
    mpfr_mul(pi, factor, qf, MPFR_RNDN);
    mpfr_div(pi, pi, tf, MPFR_RNDN);

    double t2 = now_seconds();

    printf("method=ramanujan_openmp_exact_bs\n");
    printf("digits=%ld\n", digits);
    printf("terms=%ld\n", terms);
    printf("task_depth=%d\n", max_depth);
    printf("leaf_size=%ld\n", g_leaf_size);
    printf("task_min_size=%ld\n", g_task_min_size);
#ifdef _OPENMP
    printf("omp_threads=%d\n", omp_get_max_threads());
#else
    printf("omp_threads=1\n");
#endif
    printf("binary_splitting_time=%.6f s\n", t1 - t0);
    printf("final_eval_time=%.6f s\n", t2 - t1);
    printf("total_time=%.6f s\n", t2 - t0);
    mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);

    mpz_clears(P, Q, T, NULL);
    mpfr_clears(pi, factor, qf, tf, sqrt2, (mpfr_ptr)0);
    mpfr_free_cache();
    return 0;
}

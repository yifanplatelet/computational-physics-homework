#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define A_CONST 13591409UL
#define B_CONST 545140134UL
#define C3_OVER_24_CONST 10939058860032000UL

static long g_leaf_size = 256;
static long g_task_min_size = 4096;
static long g_merge_task_min_size = 65536;
static int g_merge_task_depth = 2;

static double now_seconds(void){
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

/*
 * Leaf recurrence, fully exact.
 *
 * Descending formulation:
 *   T <- p(k) * (T +/- (A + Bk) * Q)
 *   Q <- Q * q(k)
 *   P <- P * p(k)                  [PQT only]
 *
 * This avoids the explicit lhs = p*T and rhs = p*(A+Bk)Q temporaries from the
 * older version and removes several full mpz temporaries and ui-muls per term.
 */
static void bs_leaf_pqt(mpz_t P, mpz_t Q, mpz_t T, long a, long b){
    mpz_set_ui(P, 1);
    mpz_set_ui(Q, 1);
    mpz_set_ui(T, 0);

    for (long k = b - 1; k >= a; --k){
        if (k == 0){
            mpz_addmul_ui(T, Q, A_CONST);
            continue;
        }

        unsigned long uk = (unsigned long)k;
        unsigned long f1 = 6UL * uk - 5UL;
        unsigned long f2 = 2UL * uk - 1UL;
        unsigned long f3 = 6UL * uk - 1UL;
        unsigned long linear = A_CONST + B_CONST * uk;

        if (uk & 1UL){
            mpz_submul_ui(T, Q, linear);
        }else{
            mpz_addmul_ui(T, Q, linear);
        }

        mpz_mul_ui(T, T, f1);
        mpz_mul_ui(T, T, f2);
        mpz_mul_ui(T, T, f3);

        mpz_mul_ui(P, P, f1);
        mpz_mul_ui(P, P, f2);
        mpz_mul_ui(P, P, f3);

        mpz_mul_ui(Q, Q, uk);
        mpz_mul_ui(Q, Q, uk);
        mpz_mul_ui(Q, Q, uk);
        mpz_mul_ui(Q, Q, C3_OVER_24_CONST);
    }
}

/* Same leaf as above, but skip P entirely. This is used on the right spine where
 * the aggregate P of the subtree is never consumed by the parent. */
static void bs_leaf_qt(mpz_t Q, mpz_t T, long a, long b){
    mpz_set_ui(Q, 1);
    mpz_set_ui(T, 0);

    for (long k = b - 1; k >= a; --k){
        if (k == 0){
            mpz_addmul_ui(T, Q, A_CONST);
            continue;
        }

        unsigned long uk = (unsigned long)k;
        unsigned long f1 = 6UL * uk - 5UL;
        unsigned long f2 = 2UL * uk - 1UL;
        unsigned long f3 = 6UL * uk - 1UL;
        unsigned long linear = A_CONST + B_CONST * uk;

        if (uk & 1UL){
            mpz_submul_ui(T, Q, linear);
        }else{
            mpz_addmul_ui(T, Q, linear);
        }

        mpz_mul_ui(T, T, f1);
        mpz_mul_ui(T, T, f2);
        mpz_mul_ui(T, T, f3);

        mpz_mul_ui(Q, Q, uk);
        mpz_mul_ui(Q, Q, uk);
        mpz_mul_ui(Q, Q, uk);
        mpz_mul_ui(Q, Q, C3_OVER_24_CONST);
    }
}

static void bs_serial_pqt(mpz_t P, mpz_t Q, mpz_t T, long a, long b);
static void bs_serial_qt(mpz_t Q, mpz_t T, long a, long b);

static void bs_serial_pqt(mpz_t P, mpz_t Q, mpz_t T, long a, long b){
    long n = b - a;
    if (n <= g_leaf_size){
        bs_leaf_pqt(P, Q, T, a, b);
        return;
    }

    long m = a + n / 2;
    mpz_t P1, Q1, T1, P2, Q2, T2;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, NULL);

    bs_serial_pqt(P1, Q1, T1, a, m);
    bs_serial_pqt(P2, Q2, T2, m, b);

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(T, T1, Q2);
    mpz_addmul(T, P1, T2);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, NULL);
}

static void bs_serial_qt(mpz_t Q, mpz_t T, long a, long b){
    long n = b - a;
    if (n <= g_leaf_size){
        bs_leaf_qt(Q, T, a, b);
        return;
    }

    long m = a + n / 2;
    mpz_t P1, Q1, T1, Q2, T2;
    mpz_inits(P1, Q1, T1, Q2, T2, NULL);

    bs_serial_pqt(P1, Q1, T1, a, m);
    bs_serial_qt(Q2, T2, m, b);

    mpz_mul(Q, Q1, Q2);
    mpz_mul(T, T1, Q2);
    mpz_addmul(T, P1, T2);

    mpz_clears(P1, Q1, T1, Q2, T2, NULL);
}

static void bs_parallel_pqt(mpz_t P, mpz_t Q, mpz_t T, long a, long b, int depth, int max_depth);
static void bs_parallel_qt(mpz_t Q, mpz_t T, long a, long b, int depth, int max_depth);

static void merge_pqt(mpz_t P, mpz_t Q, mpz_t T,
                      mpz_t P1, mpz_t Q1, mpz_t T1,
                      mpz_t P2, mpz_t Q2, mpz_t T2,
                      long n, int depth){
    int do_merge_tasks = (depth < g_merge_task_depth) && (n >= g_merge_task_min_size);
    if (!do_merge_tasks){
        mpz_mul(P, P1, P2);
        mpz_mul(Q, Q1, Q2);
        mpz_mul(T, T1, Q2);
        mpz_addmul(T, P1, T2);
        return;
    }

    mpz_t tmp;
    mpz_init(tmp);

    #pragma omp task shared(P, P1, P2)
    mpz_mul(P, P1, P2);

    #pragma omp task shared(Q, Q1, Q2)
    mpz_mul(Q, Q1, Q2);

    #pragma omp task shared(T, T1, Q2)
    mpz_mul(T, T1, Q2);

    #pragma omp task shared(tmp, P1, T2)
    mpz_mul(tmp, P1, T2);

    #pragma omp taskwait
    mpz_add(T, T, tmp);

    mpz_clear(tmp);
}

static void merge_qt(mpz_t Q, mpz_t T,
                     mpz_t P1, mpz_t Q1, mpz_t T1,
                     mpz_t Q2, mpz_t T2,
                     long n, int depth){
    int do_merge_tasks = (depth < g_merge_task_depth) && (n >= g_merge_task_min_size);
    if (!do_merge_tasks){
        mpz_mul(Q, Q1, Q2);
        mpz_mul(T, T1, Q2);
        mpz_addmul(T, P1, T2);
        return;
    }

    mpz_t tmp;
    mpz_init(tmp);

    #pragma omp task shared(Q, Q1, Q2)
    mpz_mul(Q, Q1, Q2);

    #pragma omp task shared(T, T1, Q2)
    mpz_mul(T, T1, Q2);

    #pragma omp task shared(tmp, P1, T2)
    mpz_mul(tmp, P1, T2);

    #pragma omp taskwait
    mpz_add(T, T, tmp);

    mpz_clear(tmp);
}

static void bs_parallel_pqt(mpz_t P, mpz_t Q, mpz_t T, long a, long b, int depth, int max_depth){
    long n = b - a;
    if (n <= g_leaf_size){
        bs_leaf_pqt(P, Q, T, a, b);
        return;
    }
    if (depth >= max_depth || n < g_task_min_size){
        bs_serial_pqt(P, Q, T, a, b);
        return;
    }

    long m = a + n / 2;
    mpz_t P1, Q1, T1, P2, Q2, T2;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, NULL);

    #pragma omp task shared(P1, Q1, T1) firstprivate(a, m, depth, max_depth)
    bs_parallel_pqt(P1, Q1, T1, a, m, depth + 1, max_depth);

    bs_parallel_pqt(P2, Q2, T2, m, b, depth + 1, max_depth);

    #pragma omp taskwait

    merge_pqt(P, Q, T, P1, Q1, T1, P2, Q2, T2, n, depth);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, NULL);
}

static void bs_parallel_qt(mpz_t Q, mpz_t T, long a, long b, int depth, int max_depth){
    long n = b - a;
    if (n <= g_leaf_size){
        bs_leaf_qt(Q, T, a, b);
        return;
    }
    if (depth >= max_depth || n < g_task_min_size){
        bs_serial_qt(Q, T, a, b);
        return;
    }

    long m = a + n / 2;
    mpz_t P1, Q1, T1, Q2, T2;
    mpz_inits(P1, Q1, T1, Q2, T2, NULL);

    #pragma omp task shared(P1, Q1, T1) firstprivate(a, m, depth, max_depth)
    bs_parallel_pqt(P1, Q1, T1, a, m, depth + 1, max_depth);

    bs_parallel_qt(Q2, T2, m, b, depth + 1, max_depth);

    #pragma omp taskwait

    merge_qt(Q, T, P1, Q1, T1, Q2, T2, n, depth);

    mpz_clears(P1, Q1, T1, Q2, T2, NULL);
}

int main(int argc, char **argv){
    long digits = 100000000;
    int max_depth = 12;

    if (argc >= 2) digits = atol(argv[1]);
    if (argc >= 3) max_depth = atoi(argv[2]);
    if (argc >= 4) g_leaf_size = atol(argv[3]);
    if (argc >= 5) g_task_min_size = atol(argv[4]);
    if (argc >= 6) g_merge_task_depth = atoi(argv[5]);
    if (argc >= 7) g_merge_task_min_size = atol(argv[6]);

    if (digits <= 0 || max_depth < 0 || g_leaf_size <= 0 || g_task_min_size <= 0 || g_merge_task_depth < 0 || g_merge_task_min_size <= 0){
        fprintf(stderr, "Usage: %s [digits] [task_depth] [leaf_size] [task_min_size] [merge_task_depth] [merge_task_min_size]\n", argv[0]);
        return 1;
    }

#ifdef _OPENMP
    omp_set_dynamic(0);
#endif

    long terms = digits / 14 + 10;
    mpfr_prec_t prec = (mpfr_prec_t)(digits * 3.32192809488736234787) + 256;

    mpz_t Q, T;
    mpz_inits(Q, T, NULL);

    double t0 = now_seconds();

    #pragma omp parallel
    {
        #pragma omp single nowait
        bs_parallel_qt(Q, T, 0, terms, 0, max_depth);
    }

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

    printf("method=chudnovsky_exact_bs_openmp_yc_learned\n");
    printf("digits=%ld\n", digits);
    printf("terms=%ld\n", terms);
    printf("task_depth=%d\n", max_depth);
    printf("leaf_size=%ld\n", g_leaf_size);
    printf("task_min_size=%ld\n", g_task_min_size);
    printf("merge_task_depth=%d\n", g_merge_task_depth);
    printf("merge_task_min_size=%ld\n", g_merge_task_min_size);
#ifdef _OPENMP
    printf("omp_threads=%d\n", omp_get_max_threads());
#else
    printf("omp_threads=1 (OpenMP disabled at compile time)\n");
#endif
    printf("binary_splitting_time=%.6f s\n", t1 - t0);
    printf("final_eval_time=%.6f s\n", t2 - t1);
    printf("total_time=%.6f s\n", t2 - t0);
    mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);

    mpz_clears(Q, T, NULL);
    mpfr_clears(pi, sqrt10005, factor, qf, tf, (mpfr_ptr)0);
    mpfr_free_cache();
    return 0;
}

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

static long g_leaf_size = 128;
static long g_task_min_size = 4096;
static int  g_merge_task_depth_limit = 2;
static long g_merge_task_min_size = 65536;

typedef struct {
    mpz_t lhs;
    mpz_t rhs;
    int initialized;
} thread_scratch_t;

static thread_scratch_t* g_leaf_scratch = NULL;
static int g_scratch_count = 0;

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static void init_leaf_scratch(int threads) {
    if (threads <= 0) {
        threads = 1;
    }
    g_leaf_scratch = (thread_scratch_t*)calloc((size_t)threads, sizeof(thread_scratch_t));
    if (!g_leaf_scratch) {
        fprintf(stderr, "failed to allocate leaf scratch table\n");
        exit(1);
    }
    g_scratch_count = threads;
    for (int i = 0; i < threads; ++i) {
        mpz_init(g_leaf_scratch[i].lhs);
        mpz_init(g_leaf_scratch[i].rhs);
        g_leaf_scratch[i].initialized = 1;
    }
}

static void destroy_leaf_scratch(void) {
    if (!g_leaf_scratch) {
        return;
    }
    for (int i = 0; i < g_scratch_count; ++i) {
        if (g_leaf_scratch[i].initialized) {
            mpz_clear(g_leaf_scratch[i].lhs);
            mpz_clear(g_leaf_scratch[i].rhs);
        }
    }
    free(g_leaf_scratch);
    g_leaf_scratch = NULL;
    g_scratch_count = 0;
}

static inline thread_scratch_t* get_leaf_scratch(void) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    if (tid < 0 || tid >= g_scratch_count) {
        tid = 0;
    }
    return &g_leaf_scratch[tid];
}

static void bs_leaf_fast(mpz_t P, mpz_t Q, mpz_t T, long a, long b) {
    thread_scratch_t* scratch = get_leaf_scratch();
    mpz_t lhs, rhs;
    lhs[0] = scratch->lhs[0];
    rhs[0] = scratch->rhs[0];

    mpz_set_ui(P, 1);
    mpz_set_ui(Q, 1);
    mpz_set_ui(T, 0);

    for (long k = b - 1; k >= a; --k) {
        if (k == 0) {
            mpz_addmul_ui(T, Q, A_CONST);
            continue;
        }

        unsigned long uk = (unsigned long)k;
        unsigned long f1 = 6UL * uk - 5UL;
        unsigned long f2 = 2UL * uk - 1UL;
        unsigned long f3 = 6UL * uk - 1UL;
        unsigned long linear = A_CONST + B_CONST * uk;

        mpz_set(lhs, T);
        mpz_mul_ui(lhs, lhs, f1);
        mpz_mul_ui(lhs, lhs, f2);
        mpz_mul_ui(lhs, lhs, f3);

        mpz_set(rhs, Q);
        mpz_mul_ui(rhs, rhs, linear);
        mpz_mul_ui(rhs, rhs, f1);
        mpz_mul_ui(rhs, rhs, f2);
        mpz_mul_ui(rhs, rhs, f3);
        if (uk & 1UL) {
            mpz_neg(rhs, rhs);
        }

        mpz_add(T, lhs, rhs);

        mpz_mul_ui(P, P, f1);
        mpz_mul_ui(P, P, f2);
        mpz_mul_ui(P, P, f3);

        mpz_mul_ui(Q, Q, uk);
        mpz_mul_ui(Q, Q, uk);
        mpz_mul_ui(Q, Q, uk);
        mpz_mul_ui(Q, Q, C3_OVER_24_CONST);
    }
}

static void bs_serial(mpz_t P, mpz_t Q, mpz_t T, long a, long b) {
    long n = b - a;
    if (n <= g_leaf_size) {
        bs_leaf_fast(P, Q, T, a, b);
        return;
    }

    long m = a + n / 2;
    mpz_t P1, Q1, T1, P2, Q2, T2, tmp;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, tmp, NULL);

    bs_serial(P1, Q1, T1, a, m);
    bs_serial(P2, Q2, T2, m, b);

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(T, T1, Q2);
    mpz_mul(tmp, P1, T2);
    mpz_add(T, T, tmp);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, tmp, NULL);
}

static void bs_parallel(mpz_t P, mpz_t Q, mpz_t T, long a, long b, int depth, int max_depth) {
    long n = b - a;
    if (n <= g_leaf_size) {
        bs_leaf_fast(P, Q, T, a, b);
        return;
    }
    if (depth >= max_depth || n < g_task_min_size) {
        bs_serial(P, Q, T, a, b);
        return;
    }

    long m = a + n / 2;
    mpz_t P1, Q1, T1, P2, Q2, T2, tmp;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, tmp, NULL);

#ifdef _OPENMP
    #pragma omp task shared(P1, Q1, T1) firstprivate(a, m, depth, max_depth)
#endif
    bs_parallel(P1, Q1, T1, a, m, depth + 1, max_depth);

    bs_parallel(P2, Q2, T2, m, b, depth + 1, max_depth);

#ifdef _OPENMP
    #pragma omp taskwait
#endif

    int do_merge_tasks = (depth <= g_merge_task_depth_limit && n >= g_merge_task_min_size);

#ifdef _OPENMP
    if (do_merge_tasks) {
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
    } else
#endif
    {
        mpz_mul(P, P1, P2);
        mpz_mul(Q, Q1, Q2);
        mpz_mul(T, T1, Q2);
        mpz_mul(tmp, P1, T2);
        mpz_add(T, T, tmp);
    }

    mpz_clears(P1, Q1, T1, P2, Q2, T2, tmp, NULL);
}

int main(int argc, char **argv) {
    long digits = 1000000;
    int max_depth = 10;

    if (argc >= 2) digits = atol(argv[1]);
    if (argc >= 3) max_depth = atoi(argv[2]);
    if (argc >= 4) g_leaf_size = atol(argv[3]);
    if (argc >= 5) g_task_min_size = atol(argv[4]);
    if (argc >= 6) g_merge_task_depth_limit = atoi(argv[5]);
    if (argc >= 7) g_merge_task_min_size = atol(argv[6]);

    if (digits <= 0 || max_depth < 0 || g_leaf_size <= 0 || g_task_min_size <= 0 ||
        g_merge_task_depth_limit < 0 || g_merge_task_min_size <= 0) {
        fprintf(stderr, "Usage: %s [digits] [task_depth] [leaf_size] [task_min_size] [merge_task_depth_limit] [merge_task_min_size]\n", argv[0]);
        return 1;
    }

#ifdef _OPENMP
    omp_set_dynamic(0);
#endif

    long terms = digits / 14 + 10;
    mpfr_prec_t prec = (mpfr_prec_t)(digits * 3.32192809488736234787) + 256;

    mpz_t P, Q, T;
    mpz_inits(P, Q, T, NULL);

    int thread_count = 1;
#ifdef _OPENMP
    thread_count = omp_get_max_threads();
#endif
    init_leaf_scratch(thread_count);

    double t0 = now_seconds();

#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single nowait
        bs_parallel(P, Q, T, 0, terms, 0, max_depth);
    }
#else
    bs_parallel(P, Q, T, 0, terms, 0, max_depth);
#endif

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

    printf("method=chudnovsky_exact_bs_openmp_yc_inspired\n");
    printf("digits=%ld\n", digits);
    printf("terms=%ld\n", terms);
    printf("task_depth=%d\n", max_depth);
    printf("leaf_size=%ld\n", g_leaf_size);
    printf("task_min_size=%ld\n", g_task_min_size);
    printf("merge_task_depth_limit=%d\n", g_merge_task_depth_limit);
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

    destroy_leaf_scratch();
    mpz_clears(P, Q, T, NULL);
    mpfr_clears(pi, sqrt10005, factor, qf, tf, (mpfr_ptr)0);
    mpfr_free_cache();
    return 0;
}

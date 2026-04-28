#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define A_CONST 13591409UL
#define B_CONST 545140134UL
#define C3_OVER_24_CONST 10939058860032000UL

static long g_leaf_size = 128;
static long g_task_min_size = 4096;
static mpz_t *tls_lhs;
static mpz_t *tls_rhs;

static double now_seconds(void){
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static void bs_leaf_fast(mpz_t P, mpz_t Q, mpz_t T, long a, long b){
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    mpz_ptr lhs = tls_lhs[tid];
    mpz_ptr rhs = tls_rhs[tid];

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

        mpz_set(lhs, T);
        mpz_mul_ui(lhs, lhs, f1);
        mpz_mul_ui(lhs, lhs, f2);
        mpz_mul_ui(lhs, lhs, f3);

        mpz_set(rhs, Q);
        mpz_mul_ui(rhs, rhs, linear);
        mpz_mul_ui(rhs, rhs, f1);
        mpz_mul_ui(rhs, rhs, f2);
        mpz_mul_ui(rhs, rhs, f3);
        if (uk & 1UL){
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

static void bs_serial(mpz_t P, mpz_t Q, mpz_t T, long a, long b){
    long n = b - a;
    if (n <= g_leaf_size){
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

static void bs_parallel(mpz_t P, mpz_t Q, mpz_t T, long a, long b, int depth, int max_depth){
    long n = b - a;
    if (n <= g_leaf_size){
        bs_leaf_fast(P, Q, T, a, b);
        return;
    }
    if (depth >= max_depth || n < g_task_min_size){
        bs_serial(P, Q, T, a, b);
        return;
    }

    long m = a + n / 2;
    mpz_t P1, Q1, T1, P2, Q2, T2, tmp;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, tmp, NULL);

    #pragma omp task shared(P1, Q1, T1) firstprivate(a, m, depth, max_depth)
    bs_parallel(P1, Q1, T1, a, m, depth + 1, max_depth);

    bs_parallel(P2, Q2, T2, m, b, depth + 1, max_depth);

    #pragma omp taskwait

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

    mpz_clears(P1, Q1, T1, P2, Q2, T2, tmp, NULL);
}

int main(int argc, char **argv){
    long digits = 100000000;
    int max_depth = 12;
    int threads = 1;

    if (argc >= 2) digits = atol(argv[1]);
    if (argc >= 3) max_depth = atoi(argv[2]);
    if (argc >= 4) g_leaf_size = atol(argv[3]);
    if (argc >= 5) g_task_min_size = atol(argv[4]);

#ifdef _OPENMP
    omp_set_dynamic(0);
    threads = omp_get_max_threads();
#endif

    tls_lhs = malloc((size_t)threads * sizeof(mpz_t));
    tls_rhs = malloc((size_t)threads * sizeof(mpz_t));
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        mpz_init(tls_lhs[tid]);
        mpz_init(tls_rhs[tid]);
    }

    long terms = (long)(digits / 14.181647462725477) + 10;
    mp_bitcnt_t prec = (mp_bitcnt_t)(digits * 3.32192809488736234787) + 64;

    mpz_t P, Q, T;
    mpz_inits(P, Q, T, NULL);

    double t0 = now_seconds();
    #pragma omp parallel
    {
        #pragma omp single nowait
        bs_parallel(P, Q, T, 0, terms, 0, max_depth);
    }
    double t1 = now_seconds();

    mpf_t qf, tf, sqrt10005, factor, pi;
    mpf_init2(qf, prec);
    mpf_init2(tf, prec);
    mpf_init2(sqrt10005, prec);
    mpf_init2(factor, prec);
    mpf_init2(pi, prec);

    mpf_set_z(qf, Q);
    mpf_set_z(tf, T);
    mpf_sqrt_ui(sqrt10005, 10005UL);
    mpf_mul_ui(factor, sqrt10005, 426880UL);
    mpf_mul(pi, factor, qf);
    mpf_div(pi, pi, tf);

    double t2 = now_seconds();

    printf("method=chudnovsky_extreme_mpf\n");
    printf("digits=%ld\n", digits);
    printf("terms=%ld\n", terms);
    printf("threads=%d\n", threads);
    printf("task_depth=%d\n", max_depth);
    printf("leaf_size=%ld\n", g_leaf_size);
    printf("task_min_size=%ld\n", g_task_min_size);
    printf("binary_splitting_time=%.6f s\n", t1 - t0);
    printf("final_eval_time=%.6f s\n", t2 - t1);
    printf("total_time=%.6f s\n", t2 - t0);
    gmp_printf("pi(first 50 digits)=%.50Ff\n", pi);

    mpz_clears(P, Q, T, NULL);
    mpf_clears(qf, tf, sqrt10005, factor, pi, NULL);
    for (int i = 0; i < threads; ++i){
        mpz_clear(tls_lhs[i]);
        mpz_clear(tls_rhs[i]);
    }
    free(tls_lhs);
    free(tls_rhs);
    return 0;
}

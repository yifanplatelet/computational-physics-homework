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

static long g_leaf_size = 256; // 增大叶子大小，减少递归开销
static long g_task_min_size = 4096;

// 线程专属工作区，彻底消除叶子节点内部的 malloc
mpz_t *tls_lhs;
mpz_t *tls_rhs;

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

// 极限优化的叶子节点
static void bs_leaf_fast(mpz_t P, mpz_t Q, mpz_t T, long a, long b) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num(); // 获取当前执行任务的线程 ID
#endif

    mpz_ptr lhs = tls_lhs[tid]; // 直接使用该线程的专属内存，不申请！
    mpz_ptr rhs = tls_rhs[tid];

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

        /* lhs = p(k) * T */
        mpz_set(lhs, T);
        mpz_mul_ui(lhs, lhs, f1);
        mpz_mul_ui(lhs, lhs, f2);
        mpz_mul_ui(lhs, lhs, f3);

        /* rhs = (A + Bk) * p(k) * Q */
        mpz_set(rhs, Q);
        mpz_mul_ui(rhs, rhs, linear);
        mpz_mul_ui(rhs, rhs, f1);
        mpz_mul_ui(rhs, rhs, f2);
        mpz_mul_ui(rhs, rhs, f3);
        
        if (uk & 1UL) {
            mpz_neg(rhs, rhs);
        }

        mpz_add(T, lhs, rhs);

        /* P *= p(k) */
        mpz_mul_ui(P, P, f1);
        mpz_mul_ui(P, P, f2);
        mpz_mul_ui(P, P, f3);

        /* Q *= q(k) */
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

    #pragma omp task shared(P1, Q1, T1) firstprivate(a, m, depth, max_depth)
    bs_parallel(P1, Q1, T1, a, m, depth + 1, max_depth);

    bs_parallel(P2, Q2, T2, m, b, depth + 1, max_depth);

    #pragma omp taskwait

    // 顶层 4 大乘法并行
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

int main(int argc, char **argv) {
    long digits = 200000000; // 默认挑战 2 亿位
    int max_depth = 12;      // 深任务树，喂饱 22 核

    if (argc >= 2) digits = atol(argv[1]);

    long terms = digits / 14 + 10;
    mpfr_prec_t prec = (mpfr_prec_t)(digits * 3.321928094887) + 256;

    int max_threads = 1;
#ifdef _OPENMP
    omp_set_dynamic(0);
    max_threads = omp_get_max_threads();
#endif

    // 初始化每个线程的专属工作区（消除叶子节点内存分配）
    tls_lhs = malloc(max_threads * sizeof(mpz_t));
    tls_rhs = malloc(max_threads * sizeof(mpz_t));
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        mpz_init(tls_lhs[tid]);
        mpz_init(tls_rhs[tid]);
    }

    mpz_t P, Q, T;
    mpz_inits(P, Q, T, NULL);

    double t0 = now_seconds();

    #pragma omp parallel
    {
        #pragma omp single nowait
        bs_parallel(P, Q, T, 0, terms, 0, max_depth);
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

    printf("method=chudnovsky_extreme_bs_22core\n");
    printf("digits=%ld\n", digits);
    printf("threads=%d\n", max_threads);
    printf("binary_splitting_time=%.6f s\n", t1 - t0);
    printf("final_eval_time=%.6f s\n", t2 - t1);
    printf("total_time=%.6f s\n", t2 - t0);
    mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);

    mpz_clears(P, Q, T, NULL);
    mpfr_clears(pi, sqrt10005, factor, qf, tf, (mpfr_ptr)0);
    
    // 清理专属工作区
    for (int i = 0; i < max_threads; i++) {
        mpz_clear(tls_lhs[i]);
        mpz_clear(tls_rhs[i]);
    }
    free(tls_lhs);
    free(tls_rhs);
    
    mpfr_free_cache();
    return 0;
}
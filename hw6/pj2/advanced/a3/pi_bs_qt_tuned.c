#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <gmp.h>
#include <mpfr.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define A_CONST 13591409UL
#define B_CONST 545140134UL
#define C3_OVER_24_CONST 10939058860032000UL

static long g_leaf_size = 256;
static long g_task_min_size = 16384;
static long g_merge_task_min_size = 262144;
static int g_merge_task_depth = 1;
static int g_print_pi = 1;

static double now_seconds(void){
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static unsigned long bitlen_u64(uint64_t x){
    if (x == 0){
        return 1;
    }
    return 64UL - (unsigned long)__builtin_clzll(x);
}

static mp_bitcnt_t ceil_log2_10_digits(mp_bitcnt_t digits){
    return (digits * 332193ULL + 99999ULL) / 100000ULL + 1ULL;
}

static void reserve_leaf_pqt(mpz_t P, mpz_t Q, mpz_t T, long a, long b){
    long n = b - a;
    uint64_t ub = b > 0 ? (uint64_t)b : 1ULL;
    unsigned long p_term_bits = 3UL * bitlen_u64(6ULL * ub + 1ULL);
    unsigned long q_term_bits = 3UL * bitlen_u64(ub) + bitlen_u64(C3_OVER_24_CONST);
    unsigned long t_term_bits = p_term_bits + q_term_bits + bitlen_u64(A_CONST + B_CONST * ub) + 2UL;

    mpz_realloc2(P, (mp_bitcnt_t)n * p_term_bits + 64);
    mpz_realloc2(Q, (mp_bitcnt_t)n * q_term_bits + 64);
    mpz_realloc2(T, (mp_bitcnt_t)n * t_term_bits + 64);
}

static void reserve_leaf_qt(mpz_t Q, mpz_t T, long a, long b){
    long n = b - a;
    uint64_t ub = b > 0 ? (uint64_t)b : 1ULL;
    unsigned long q_term_bits = 3UL * bitlen_u64(ub) + bitlen_u64(C3_OVER_24_CONST);
    unsigned long t_term_bits = 3UL * bitlen_u64(6ULL * ub + 1ULL) + q_term_bits + bitlen_u64(A_CONST + B_CONST * ub) + 2UL;

    mpz_realloc2(Q, (mp_bitcnt_t)n * q_term_bits + 64);
    mpz_realloc2(T, (mp_bitcnt_t)n * t_term_bits + 64);
}

static void reserve_mul_result(mpz_t dst, const mpz_t a, const mpz_t b){
    mp_size_t limbs = mpz_size(a) + mpz_size(b) + 1;
    if (limbs < 1){
        limbs = 1;
    }
    mpz_realloc2(dst, (mp_bitcnt_t)limbs * GMP_NUMB_BITS);
}

static void reserve_add_products_result(mpz_t dst,
                                        const mpz_t a1, const mpz_t b1,
                                        const mpz_t a2, const mpz_t b2){
    mp_size_t limbs1 = mpz_size(a1) + mpz_size(b1) + 1;
    mp_size_t limbs2 = mpz_size(a2) + mpz_size(b2) + 1;
    mp_size_t limbs = limbs1 > limbs2 ? limbs1 : limbs2;
    limbs += 1;
    if (limbs < 1){
        limbs = 1;
    }
    mpz_realloc2(dst, (mp_bitcnt_t)limbs * GMP_NUMB_BITS);
}

static void bs_leaf_pqt(mpz_t P, mpz_t Q, mpz_t T, long a, long b){
    reserve_leaf_pqt(P, Q, T, a, b);
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

static void bs_leaf_qt(mpz_t Q, mpz_t T, long a, long b){
    reserve_leaf_qt(Q, T, a, b);
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

    reserve_mul_result(P, P1, P2);
    reserve_mul_result(Q, Q1, Q2);
    reserve_add_products_result(T, T1, Q2, P1, T2);

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

    reserve_mul_result(Q, Q1, Q2);
    reserve_add_products_result(T, T1, Q2, P1, T2);

    mpz_mul(Q, Q1, Q2);
    mpz_mul(T, T1, Q2);
    mpz_addmul(T, P1, T2);

    mpz_clears(P1, Q1, T1, Q2, T2, NULL);
}

static void merge_pqt(mpz_t P, mpz_t Q, mpz_t T,
                      mpz_t P1, mpz_t Q1, mpz_t T1,
                      mpz_t P2, mpz_t Q2, mpz_t T2,
                      long n, int depth){
    int do_merge_tasks = 0;
#ifdef _OPENMP
    do_merge_tasks = (depth < g_merge_task_depth) && (n >= g_merge_task_min_size);
#endif

    reserve_mul_result(P, P1, P2);
    reserve_mul_result(Q, Q1, Q2);
    reserve_add_products_result(T, T1, Q2, P1, T2);

    if (!do_merge_tasks){
        mpz_mul(P, P1, P2);
        mpz_mul(Q, Q1, Q2);
        mpz_mul(T, T1, Q2);
        mpz_addmul(T, P1, T2);
        return;
    }

    mpz_t tmp;
    mpz_init2(tmp, (mp_bitcnt_t)(mpz_size(P1) + mpz_size(T2) + 1) * GMP_NUMB_BITS);

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
    int do_merge_tasks = 0;
#ifdef _OPENMP
    do_merge_tasks = (depth < g_merge_task_depth) && (n >= g_merge_task_min_size);
#endif

    reserve_mul_result(Q, Q1, Q2);
    reserve_add_products_result(T, T1, Q2, P1, T2);

    if (!do_merge_tasks){
        mpz_mul(Q, Q1, Q2);
        mpz_mul(T, T1, Q2);
        mpz_addmul(T, P1, T2);
        return;
    }

    mpz_t tmp;
    mpz_init2(tmp, (mp_bitcnt_t)(mpz_size(P1) + mpz_size(T2) + 1) * GMP_NUMB_BITS);

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

static void bs_parallel_pqt(mpz_t P, mpz_t Q, mpz_t T, long a, long b, int depth, int max_depth);
static void bs_parallel_qt(mpz_t Q, mpz_t T, long a, long b, int depth, int max_depth);

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

static int auto_task_depth(int threads){
    int depth = 0;
    int tasks = 1;
    if (threads < 1){
        threads = 1;
    }
    while (tasks < threads * 2){
        tasks <<= 1;
        depth++;
    }
    return depth;
}

static void compute_root_qt(mpz_t Q, mpz_t T, long terms, int max_depth){
    if (terms <= g_leaf_size){
        mpz_t Ptmp;
        mpz_init(Ptmp);
        bs_leaf_pqt(Ptmp, Q, T, 0, terms);
        mpz_clear(Ptmp);
        return;
    }

    long m = terms / 2;
    mpz_t P1, Q1, T1, Q2, T2;
    mpz_inits(P1, Q1, T1, Q2, T2, NULL);

    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            #pragma omp task shared(P1, Q1, T1) firstprivate(m, max_depth)
            bs_parallel_pqt(P1, Q1, T1, 0, m, 1, max_depth);

            bs_parallel_qt(Q2, T2, m, terms, 1, max_depth);

            #pragma omp taskwait
        }
    }

    reserve_mul_result(Q, Q1, Q2);
    reserve_add_products_result(T, T1, Q2, P1, T2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(T, T1, Q2);
    mpz_addmul(T, P1, T2);

    mpz_clears(P1, Q1, T1, Q2, T2, NULL);
}

static void final_eval_mpfr(mpfr_t pi, const mpz_t Q, const mpz_t T, mpfr_prec_t prec){
    mpfr_t sqrt10005, factor, qf, tf;
    mpfr_inits2(prec, sqrt10005, factor, qf, tf, (mpfr_ptr)0);

    mpfr_set_ui(sqrt10005, 10005, MPFR_RNDN);
    mpfr_sqrt(sqrt10005, sqrt10005, MPFR_RNDN);
    mpfr_mul_ui(factor, sqrt10005, 426880, MPFR_RNDN);
    mpfr_set_z(qf, Q, MPFR_RNDN);
    mpfr_set_z(tf, T, MPFR_RNDN);
    mpfr_mul(pi, factor, qf, MPFR_RNDN);
    mpfr_div(pi, pi, tf, MPFR_RNDN);

    mpfr_clears(sqrt10005, factor, qf, tf, (mpfr_ptr)0);
}

static void usage(const char *prog){
    fprintf(stderr,
        "Usage: %s [digits] [task_depth] [leaf_size] [task_min_size] [merge_task_depth] [merge_task_min_size] [print_pi]\n"
        "  task_depth=0 means auto\n",
        prog);
}

int main(int argc, char **argv){
    long digits = 100000000;
    int max_depth = 0;
    int threads = 1;

    if (argc >= 2) digits = atol(argv[1]);
    if (argc >= 3) max_depth = atoi(argv[2]);
    if (argc >= 4) g_leaf_size = atol(argv[3]);
    if (argc >= 5) g_task_min_size = atol(argv[4]);
    if (argc >= 6) g_merge_task_depth = atoi(argv[5]);
    if (argc >= 7) g_merge_task_min_size = atol(argv[6]);
    if (argc >= 8) g_print_pi = atoi(argv[7]);

    if (digits <= 0 || max_depth < 0 || g_leaf_size <= 0 || g_task_min_size <= 0 ||
        g_merge_task_depth < 0 || g_merge_task_min_size <= 0 ||
        (g_print_pi != 0 && g_print_pi != 1)){
        usage(argv[0]);
        return 1;
    }

#ifdef _OPENMP
    omp_set_dynamic(0);
    threads = omp_get_max_threads();
#endif
    if (max_depth == 0){
        max_depth = auto_task_depth(threads);
    }

    long terms = digits / 14 + 10;
    mpfr_prec_t prec = (mpfr_prec_t)(digits * 3.32192809488736234787) + 256;
    mp_bitcnt_t approx_bits = ceil_log2_10_digits((mp_bitcnt_t)digits) + 256;

    mpz_t Q, T;
    mpz_init2(Q, approx_bits);
    mpz_init2(T, approx_bits);

    double t0 = now_seconds();
    compute_root_qt(Q, T, terms, max_depth);
    double t1 = now_seconds();

    mpfr_t pi;
    mpfr_init2(pi, prec);
    final_eval_mpfr(pi, Q, T, prec);
    double t2 = now_seconds();

    printf("method=chudnovsky_qt_tuned_prealloc\n");
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
    printf("omp_threads=1\n");
#endif
    printf("binary_splitting_time=%.6f s\n", t1 - t0);
    printf("final_eval_time=%.6f s\n", t2 - t1);
    printf("total_time=%.6f s\n", t2 - t0);
    if (g_print_pi){
        mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);
    }

    mpfr_clear(pi);
    mpz_clears(Q, T, NULL);
    mpfr_free_cache();
    return 0;
}

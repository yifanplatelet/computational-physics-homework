#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <gmp.h>
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
static long g_guard_digits = 32;
static int g_print_pi = 1;
static int g_skew_num = 0;
static int g_skew_den = 1;
static int g_skew_depth = 0;
static long g_skew_min_size = 0;

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

static long choose_split(long a, long b, int depth){
    long n = b - a;
    if (g_skew_num > 0 && g_skew_den > g_skew_num &&
        depth < g_skew_depth && n >= g_skew_min_size){
        long m = a + (long)(((__int128)n * g_skew_num) / g_skew_den);
        if (m <= a){
            m = a + 1;
        }
        if (m >= b){
            m = b - 1;
        }
        return m;
    }
    return a + n / 2;
}

static void reserve_leaf_pqt(mpz_t P, mpz_t Q, mpz_t T, long a, long b){
    long n = b - a;
    uint64_t ub = (b > 0) ? (uint64_t)b : 1ULL;
    unsigned long p_term_bits = 3UL * bitlen_u64(6ULL * ub + 1ULL);
    unsigned long q_term_bits = 3UL * bitlen_u64(ub) + bitlen_u64(C3_OVER_24_CONST);
    unsigned long t_term_bits = p_term_bits + q_term_bits + bitlen_u64(A_CONST + B_CONST * ub) + 2UL;

    mpz_realloc2(P, (mp_bitcnt_t)n * p_term_bits + 64);
    mpz_realloc2(Q, (mp_bitcnt_t)n * q_term_bits + 64);
    mpz_realloc2(T, (mp_bitcnt_t)n * t_term_bits + 64);
}

static void reserve_leaf_qt(mpz_t Q, mpz_t T, long a, long b){
    long n = b - a;
    uint64_t ub = (b > 0) ? (uint64_t)b : 1ULL;
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

/*
 * Exact compact leaf:
 *   T <- p(k) * (T +/- (A + Bk) * Q)
 *   P <- P * p(k)
 *   Q <- Q * q(k)
 */
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

static void bs_serial_pqt(mpz_t P, mpz_t Q, mpz_t T, long a, long b, int depth);
static void bs_serial_qt(mpz_t Q, mpz_t T, long a, long b, int depth);

static void bs_serial_pqt(mpz_t P, mpz_t Q, mpz_t T, long a, long b, int depth){
    long n = b - a;
    if (n <= g_leaf_size){
        bs_leaf_pqt(P, Q, T, a, b);
        return;
    }

    long m = choose_split(a, b, depth);
    mpz_t P1, Q1, T1, P2, Q2, T2;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, NULL);

    bs_serial_pqt(P1, Q1, T1, a, m, depth + 1);
    bs_serial_pqt(P2, Q2, T2, m, b, depth + 1);

    reserve_mul_result(P, P1, P2);
    reserve_mul_result(Q, Q1, Q2);
    reserve_add_products_result(T, T1, Q2, P1, T2);

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(T, T1, Q2);
    mpz_addmul(T, P1, T2);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, NULL);
}

static void bs_serial_qt(mpz_t Q, mpz_t T, long a, long b, int depth){
    long n = b - a;
    if (n <= g_leaf_size){
        bs_leaf_qt(Q, T, a, b);
        return;
    }

    long m = choose_split(a, b, depth);
    mpz_t P1, Q1, T1, Q2, T2;
    mpz_inits(P1, Q1, T1, Q2, T2, NULL);

    bs_serial_pqt(P1, Q1, T1, a, m, depth + 1);
    bs_serial_qt(Q2, T2, m, b, depth + 1);

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
        bs_serial_pqt(P, Q, T, a, b, depth);
        return;
    }

    long m = choose_split(a, b, depth);
    mpz_t P1, Q1, T1, P2, Q2, T2;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, NULL);

#ifdef _OPENMP
    #pragma omp task shared(P1, Q1, T1) firstprivate(a, m, depth, max_depth)
#endif
    bs_parallel_pqt(P1, Q1, T1, a, m, depth + 1, max_depth);

    bs_parallel_pqt(P2, Q2, T2, m, b, depth + 1, max_depth);

#ifdef _OPENMP
    #pragma omp taskwait
#endif

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
        bs_serial_qt(Q, T, a, b, depth);
        return;
    }

    long m = choose_split(a, b, depth);
    mpz_t P1, Q1, T1, Q2, T2;
    mpz_inits(P1, Q1, T1, Q2, T2, NULL);

#ifdef _OPENMP
    #pragma omp task shared(P1, Q1, T1) firstprivate(a, m, depth, max_depth)
#endif
    bs_parallel_pqt(P1, Q1, T1, a, m, depth + 1, max_depth);

    bs_parallel_qt(Q2, T2, m, b, depth + 1, max_depth);

#ifdef _OPENMP
    #pragma omp taskwait
#endif

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
    return depth + 1;
}

/*
 * Root specialization:
 * we need only Q and T for final evaluation.
 * Left subtree needs full P/Q/T, right subtree can use the cheaper QT recursion.
 */
static void compute_root_qt(mpz_t Q, mpz_t T, long terms, int max_depth){
    if (terms <= g_leaf_size){
        mpz_t Ptmp;
        mpz_init(Ptmp);
        bs_leaf_pqt(Ptmp, Q, T, 0, terms);
        mpz_clear(Ptmp);
        return;
    }

    long m = choose_split(0, terms, 0);
    mpz_t P1, Q1, T1, Q2, T2;
    mpz_inits(P1, Q1, T1, Q2, T2, NULL);

#ifdef _OPENMP
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
#else
    bs_parallel_pqt(P1, Q1, T1, 0, m, 1, max_depth);
    bs_parallel_qt(Q2, T2, m, terms, 1, max_depth);
#endif

    reserve_mul_result(Q, Q1, Q2);
    reserve_add_products_result(T, T1, Q2, P1, T2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(T, T1, Q2);
    mpz_addmul(T, P1, T2);

    mpz_clears(P1, Q1, T1, Q2, T2, NULL);
}

static void compute_pi_scaled(mpz_t pi_scaled, const mpz_t Q, const mpz_t T, long digits, long guard_digits){
    mp_bitcnt_t scale_digits = (mp_bitcnt_t)(digits + guard_digits);
    mp_bitcnt_t scale_bits = ceil_log2_10_digits(scale_digits);
    mp_bitcnt_t sqrt_arg_bits = 2 * scale_bits + 16;
    mp_bitcnt_t sqrt_bits = scale_bits + 8;
    mp_bitcnt_t q_bits = (mp_bitcnt_t)mpz_size(Q) * GMP_NUMB_BITS;

    mpz_t sqrt_arg, sqrt_scaled, numerator;
    mpz_init2(sqrt_arg, sqrt_arg_bits);
    mpz_init2(sqrt_scaled, sqrt_bits);
    mpz_init2(numerator, q_bits + sqrt_bits + 32);

    mpz_ui_pow_ui(sqrt_arg, 10UL, (unsigned long)(2 * scale_digits));
    mpz_mul_ui(sqrt_arg, sqrt_arg, 10005UL);
    mpz_sqrt(sqrt_scaled, sqrt_arg);

    mpz_mul(numerator, Q, sqrt_scaled);
    mpz_mul_ui(numerator, numerator, 426880UL);
    mpz_fdiv_q(pi_scaled, numerator, T);

    mpz_clears(sqrt_arg, sqrt_scaled, numerator, NULL);
}

static void print_pi_prefix(const mpz_t pi_scaled, long digits, long guard_digits){
    const long frac_digits = 50;
    long total_scale = digits + guard_digits;
    long trim_digits = total_scale - frac_digits;
    mpz_t prefix, divisor, rounded, half_ulp;
    char *text;
    size_t len;

    mpz_init(prefix);
    mpz_init(divisor);
    mpz_init(rounded);
    mpz_init(half_ulp);

    if (trim_digits > 0){
        mpz_ui_pow_ui(divisor, 10UL, (unsigned long)trim_digits);
        mpz_set(rounded, pi_scaled);
        mpz_ui_pow_ui(half_ulp, 10UL, (unsigned long)(trim_digits - 1));
        mpz_mul_ui(half_ulp, half_ulp, 5UL);
        mpz_add(rounded, rounded, half_ulp);
        mpz_fdiv_q(prefix, rounded, divisor);
    }else{
        mpz_set(prefix, pi_scaled);
    }

    text = mpz_get_str(NULL, 10, prefix);
    len = strlen(text);

    if (len == 0){
        printf("pi(first 50 digits)=0.00000000000000000000000000000000000000000000000000\n");
        mpz_clears(prefix, divisor, rounded, half_ulp, NULL);
        return;
    }

    printf("pi(first 50 digits)=%c", text[0]);
    if (len == 1){
        putchar('.');
        for (long i = 0; i < frac_digits; ++i){
            putchar('0');
        }
        putchar('\n');
    }else{
        putchar('.');
        fputs(text + 1, stdout);
        for (size_t i = len - 1; i < (size_t)frac_digits; ++i){
            putchar('0');
        }
        putchar('\n');
    }

    void (*freefunc)(void *, size_t);
    mp_get_memory_functions(NULL, NULL, &freefunc);
    freefunc(text, len + 1);
    mpz_clears(prefix, divisor, rounded, half_ulp, NULL);
}

static void print_usage(const char *prog){
    fprintf(stderr,
        "Usage: %s [digits] [task_depth] [leaf_size] [task_min_size] [merge_task_depth] [merge_task_min_size] [guard_digits] [print_pi] [skew_num] [skew_den] [skew_depth] [skew_min_size]\n"
        "  digits               target decimal digits, default 100000000\n"
        "  task_depth           OpenMP recursion depth, 0 means auto\n"
        "  leaf_size            exact leaf block size, default 256\n"
        "  task_min_size        min range size before spawning tasks, default 4096\n"
        "  merge_task_depth     levels that parallelize merge multiplies, default 2\n"
        "  merge_task_min_size  min range size for merge-task fanout, default 65536\n"
        "  guard_digits         extra decimal guard digits for fixed-point finish, default 32\n"
        "  print_pi             1 prints first 50 digits, 0 skips it, default 1\n"
        "  skew_num/skew_den    optional top-level skewed split ratio, default 0/1\n"
        "  skew_depth           skewed split depth, default 0\n"
        "  skew_min_size        minimum range size for skew, default 0\n",
        prog);
}

int main(int argc, char **argv){
    long digits = 100000000;
    int max_depth = 0;
    int threads = 1;
    long terms;
    mp_bitcnt_t total_digits_bits;

    if (argc >= 2) digits = atol(argv[1]);
    if (argc >= 3) max_depth = atoi(argv[2]);
    if (argc >= 4) g_leaf_size = atol(argv[3]);
    if (argc >= 5) g_task_min_size = atol(argv[4]);
    if (argc >= 6) g_merge_task_depth = atoi(argv[5]);
    if (argc >= 7) g_merge_task_min_size = atol(argv[6]);
    if (argc >= 8) g_guard_digits = atol(argv[7]);
    if (argc >= 9) g_print_pi = atoi(argv[8]);
    if (argc >= 10) g_skew_num = atoi(argv[9]);
    if (argc >= 11) g_skew_den = atoi(argv[10]);
    if (argc >= 12) g_skew_depth = atoi(argv[11]);
    if (argc >= 13) g_skew_min_size = atol(argv[12]);

    if (digits <= 0 || max_depth < 0 || g_leaf_size <= 0 || g_task_min_size <= 0 ||
        g_merge_task_depth < 0 || g_merge_task_min_size <= 0 || g_guard_digits <= 0 ||
        (g_print_pi != 0 && g_print_pi != 1) || g_skew_num < 0 || g_skew_den <= 0 ||
        g_skew_depth < 0 || g_skew_min_size < 0 || g_skew_num >= g_skew_den){
        print_usage(argv[0]);
        return 1;
    }

#ifdef _OPENMP
    omp_set_dynamic(0);
    threads = omp_get_max_threads();
#endif

    if (max_depth == 0){
        max_depth = auto_task_depth(threads);
    }

    terms = digits / 14 + 10;
    total_digits_bits = ceil_log2_10_digits((mp_bitcnt_t)(digits + g_guard_digits));

    mpz_t Q, T, pi_scaled;
    mpz_init2(Q, total_digits_bits + 64);
    mpz_init2(T, total_digits_bits + 64);
    mpz_init2(pi_scaled, total_digits_bits + 64);

    double t0 = now_seconds();
    compute_root_qt(Q, T, terms, max_depth);
    double t1 = now_seconds();

    compute_pi_scaled(pi_scaled, Q, T, digits, g_guard_digits);
    double t2 = now_seconds();

    printf("method=chudnovsky_exact_bs_gmp_kernel_tuned\n");
    printf("digits=%ld\n", digits);
    printf("terms=%ld\n", terms);
    printf("task_depth=%d\n", max_depth);
    printf("leaf_size=%ld\n", g_leaf_size);
    printf("task_min_size=%ld\n", g_task_min_size);
    printf("merge_task_depth=%d\n", g_merge_task_depth);
    printf("merge_task_min_size=%ld\n", g_merge_task_min_size);
    printf("guard_digits=%ld\n", g_guard_digits);
    printf("skew=%d/%d\n", g_skew_num, g_skew_den);
    printf("skew_depth=%d\n", g_skew_depth);
    printf("skew_min_size=%ld\n", g_skew_min_size);
#ifdef _OPENMP
    printf("omp_threads=%d\n", omp_get_max_threads());
#else
    printf("omp_threads=1 (OpenMP disabled at compile time)\n");
#endif
    printf("binary_splitting_time=%.6f s\n", t1 - t0);
    printf("fixedpoint_finish_time=%.6f s\n", t2 - t1);
    printf("total_time=%.6f s\n", t2 - t0);

    if (g_print_pi){
        print_pi_prefix(pi_scaled, digits, g_guard_digits);
    }

    mpz_clears(Q, T, pi_scaled, NULL);
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
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

/*
 * Compact exact leaf recurrence:
 *   T <- p(k) * (T +/- (A + Bk) * Q)
 *   P <- P * p(k)
 *   Q <- Q * q(k)
 *
 * This removes the scratch lhs/rhs temporaries from the original 9/main.c.
 */
static void bs_leaf_compact(mpz_t P, mpz_t Q, mpz_t T, long a, long b){
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

static void bs_serial(mpz_t P, mpz_t Q, mpz_t T, long a, long b, int depth){
    long n = b - a;
    if (n <= g_leaf_size){
        bs_leaf_compact(P, Q, T, a, b);
        return;
    }

    long m = choose_split(a, b, depth);
    mpz_t P1, Q1, T1, P2, Q2, T2;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, NULL);

    bs_serial(P1, Q1, T1, a, m, depth + 1);
    bs_serial(P2, Q2, T2, m, b, depth + 1);

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(T, T1, Q2);
    mpz_addmul(T, P1, T2);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, NULL);
}

static void merge_node(mpz_t P, mpz_t Q, mpz_t T,
                       mpz_t P1, mpz_t Q1, mpz_t T1,
                       mpz_t P2, mpz_t Q2, mpz_t T2,
                       long n, int depth){
    int do_merge_tasks = 0;
#ifdef _OPENMP
    do_merge_tasks = (depth < g_merge_task_depth) && (n >= g_merge_task_min_size);
#endif

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

static void bs_parallel(mpz_t P, mpz_t Q, mpz_t T, long a, long b, int depth, int max_depth){
    long n = b - a;
    if (n <= g_leaf_size){
        bs_leaf_compact(P, Q, T, a, b);
        return;
    }
    if (depth >= max_depth || n < g_task_min_size){
        bs_serial(P, Q, T, a, b, depth);
        return;
    }

    long m = choose_split(a, b, depth);
    mpz_t P1, Q1, T1, P2, Q2, T2;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, NULL);

#ifdef _OPENMP
    #pragma omp task shared(P1, Q1, T1) firstprivate(a, m, depth, max_depth)
#endif
    bs_parallel(P1, Q1, T1, a, m, depth + 1, max_depth);

    bs_parallel(P2, Q2, T2, m, b, depth + 1, max_depth);

#ifdef _OPENMP
    #pragma omp taskwait
#endif

    merge_node(P, Q, T, P1, Q1, T1, P2, Q2, T2, n, depth);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, NULL);
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

static void compute_pi_scaled(mpz_t pi_scaled, const mpz_t Q, const mpz_t T, long digits, long guard_digits){
    unsigned long scale_digits = (unsigned long)(digits + guard_digits);
    mpz_t sqrt_arg, sqrt_scaled, numerator;
    mpz_inits(sqrt_arg, sqrt_scaled, numerator, NULL);

    mpz_ui_pow_ui(sqrt_arg, 10UL, 2UL * scale_digits);
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
        printf("pi(first 50 digits,truncated)=0.00000000000000000000000000000000000000000000000000\n");
        mpz_clear(prefix);
        mpz_clear(divisor);
        mpz_clear(rounded);
        mpz_clear(half_ulp);
        return;
    }

    printf("pi(first 50 digits,truncated)=%c", text[0]);
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
    mpz_clear(prefix);
    mpz_clear(divisor);
    mpz_clear(rounded);
    mpz_clear(half_ulp);
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
        "  skew_num/skew_den    optional top-level skew ratio, default 0/1\n"
        "  skew_depth           skewed split depth, default 0\n"
        "  skew_min_size        minimum range size for skew, default 0\n",
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

    long terms = digits / 14 + 10;
    mpz_t P, Q, T, pi_scaled;
    mpz_inits(P, Q, T, pi_scaled, NULL);

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

    compute_pi_scaled(pi_scaled, Q, T, digits, g_guard_digits);

    double t2 = now_seconds();

    printf("method=chudnovsky_exact_bs_fixedpoint_finish\n");
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

    mpz_clears(P, Q, T, pi_scaled, NULL);
    return 0;
}

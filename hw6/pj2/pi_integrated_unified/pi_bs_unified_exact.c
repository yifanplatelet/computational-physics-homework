#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define A 13591409UL
#define B 545140134UL
#define C3_OVER_24 10939058860032000UL

typedef struct {
    long digits;
    int task_depth;
    long leaf_size;
    long task_min_size;
    int skew_num;
    int skew_den;
    int skew_depth;
    long skew_min_size;
    int endgame_mode;   /* 0 = mpfr_div, 1 = newton reciprocal */
} Config;

static Config g_cfg = {1000000, 6, 256, 4096, 0, 1, 0, 0, 0};

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static inline void term_pqt(mpz_t p, mpz_t q, mpz_t t, unsigned long k) {
    mpz_set_ui(p, 6 * k - 5);
    mpz_mul_ui(p, p, 2 * k - 1);
    mpz_mul_ui(p, p, 6 * k - 1);

    mpz_set_ui(q, k);
    mpz_mul(q, q, q);
    mpz_mul_ui(q, q, k);
    mpz_mul_ui(q, q, C3_OVER_24);

    mpz_set_ui(t, B);
    mpz_mul_ui(t, t, k);
    mpz_add_ui(t, t, A);
    mpz_mul(t, t, p);
    if (k & 1UL) {
        mpz_neg(t, t);
    }
}

static inline long choose_split(long a, long b, int depth) {
    long n = b - a;
    if (g_cfg.skew_num > 0 && g_cfg.skew_den > g_cfg.skew_num &&
        depth < g_cfg.skew_depth && n >= g_cfg.skew_min_size) {
        long m = a + (long)(((__int128)n * g_cfg.skew_num) / g_cfg.skew_den);
        if (m <= a) m = a + 1;
        if (m >= b) m = b - 1;
        return m;
    }
    return a + n / 2;
}

static void bs_leaf_reverse(mpz_t P, mpz_t Q, mpz_t T, long a, long b) {
    mpz_t p, q, t, tmp;
    mpz_inits(p, q, t, tmp, NULL);

    mpz_set_ui(P, 1);
    mpz_set_ui(Q, 1);
    mpz_set_ui(T, 0);

    for (long k = b - 1; k >= a; --k) {
        if (k == 0) {
            mpz_addmul_ui(T, Q, A);
            continue;
        }

        term_pqt(p, q, t, (unsigned long)k);

        /* T <- t*Q + p*T */
        mpz_mul(tmp, p, T);
        mpz_mul(T, t, Q);
        mpz_add(T, T, tmp);

        mpz_mul(P, P, p);
        mpz_mul(Q, Q, q);
    }

    mpz_clears(p, q, t, tmp, NULL);
}

static void bs_serial(mpz_t P, mpz_t Q, mpz_t T, long a, long b, int depth) {
    long n = b - a;
    if (n <= g_cfg.leaf_size) {
        bs_leaf_reverse(P, Q, T, a, b);
        return;
    }

    long m = choose_split(a, b, depth);
    mpz_t P1, Q1, T1, P2, Q2, T2;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, NULL);

    bs_serial(P1, Q1, T1, a, m, depth + 1);
    bs_serial(P2, Q2, T2, m, b, depth + 1);

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(T, P1, T2);
    mpz_addmul(T, T1, Q2);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, NULL);
}

static void bs_parallel(mpz_t P, mpz_t Q, mpz_t T, long a, long b, int depth) {
    long n = b - a;
    if (n <= g_cfg.leaf_size) {
        bs_leaf_reverse(P, Q, T, a, b);
        return;
    }
    if (depth >= g_cfg.task_depth || n < g_cfg.task_min_size) {
        bs_serial(P, Q, T, a, b, depth);
        return;
    }

    long m = choose_split(a, b, depth);
    mpz_t P1, Q1, T1, P2, Q2, T2;
    mpz_inits(P1, Q1, T1, P2, Q2, T2, NULL);

    #pragma omp task shared(P1, Q1, T1)
    bs_parallel(P1, Q1, T1, a, m, depth + 1);

    #pragma omp task shared(P2, Q2, T2)
    bs_parallel(P2, Q2, T2, m, b, depth + 1);

    #pragma omp taskwait

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(T, P1, T2);
    mpz_addmul(T, T1, Q2);

    mpz_clears(P1, Q1, T1, P2, Q2, T2, NULL);
}

static void final_eval_direct(mpfr_t pi, mpz_t Q, mpz_t T, mpfr_prec_t prec) {
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

static void final_eval_newton(mpfr_t pi, mpz_t Q, mpz_t T, mpfr_prec_t prec) {
    mpfr_t sqrt10005, factor, qf, tf, inv, tmp;
    mpfr_prec_t start_prec = 256;
    if (start_prec > prec) start_prec = prec;

    mpfr_inits2(prec, sqrt10005, factor, qf, tf, inv, tmp, (mpfr_ptr)0);

    mpfr_set_ui(sqrt10005, 10005, MPFR_RNDN);
    mpfr_sqrt(sqrt10005, sqrt10005, MPFR_RNDN);
    mpfr_mul_ui(factor, sqrt10005, 426880, MPFR_RNDN);
    mpfr_set_z(qf, Q, MPFR_RNDN);
    mpfr_set_z(tf, T, MPFR_RNDN);

    /* seed reciprocal */
    mpfr_prec_round(inv, start_prec, MPFR_RNDN);
    mpfr_prec_round(tmp, start_prec, MPFR_RNDN);
    mpfr_prec_round(tf, start_prec, MPFR_RNDN);
    mpfr_ui_div(inv, 1, tf, MPFR_RNDN);

    while (start_prec < prec) {
        mpfr_prec_t next_prec = start_prec * 2;
        if (next_prec > prec) next_prec = prec;
        mpfr_prec_round(inv, next_prec, MPFR_RNDN);
        mpfr_prec_round(tmp, next_prec, MPFR_RNDN);
        mpfr_prec_round(tf, next_prec, MPFR_RNDN);
        mpfr_set_z(tf, T, MPFR_RNDN);
        /* inv <- inv * (2 - tf*inv) */
        mpfr_mul(tmp, tf, inv, MPFR_RNDN);
        mpfr_ui_sub(tmp, 2, tmp, MPFR_RNDN);
        mpfr_mul(inv, inv, tmp, MPFR_RNDN);
        start_prec = next_prec;
    }

    mpfr_prec_round(qf, prec, MPFR_RNDN);
    mpfr_set_z(qf, Q, MPFR_RNDN);
    mpfr_mul(pi, factor, qf, MPFR_RNDN);
    mpfr_mul(pi, pi, inv, MPFR_RNDN);

    mpfr_clears(sqrt10005, factor, qf, tf, inv, tmp, (mpfr_ptr)0);
}

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [digits] [task_depth] [leaf_size] [task_min_size] [skew_num] [skew_den] [skew_depth] [skew_min_size] [endgame_mode]\n"
        "  endgame_mode: 0 = MPFR direct divide, 1 = Newton reciprocal\n"
        "Example baseline: %s 100000000 6 256 4096 0 1 0 0 0\n"
        "Example top-skew: %s 100000000 6 256 4096 2 3 2 131072 0\n",
        prog, prog, prog);
}

int main(int argc, char **argv) {
    if (argc >= 2) g_cfg.digits = atol(argv[1]);
    if (argc >= 3) g_cfg.task_depth = atoi(argv[2]);
    if (argc >= 4) g_cfg.leaf_size = atol(argv[3]);
    if (argc >= 5) g_cfg.task_min_size = atol(argv[4]);
    if (argc >= 6) g_cfg.skew_num = atoi(argv[5]);
    if (argc >= 7) g_cfg.skew_den = atoi(argv[6]);
    if (argc >= 8) g_cfg.skew_depth = atoi(argv[7]);
    if (argc >= 9) g_cfg.skew_min_size = atol(argv[8]);
    if (argc >= 10) g_cfg.endgame_mode = atoi(argv[9]);

    if (g_cfg.digits <= 0 || g_cfg.task_depth < 0 || g_cfg.leaf_size <= 0 || g_cfg.task_min_size <= 0 ||
        g_cfg.skew_den <= 0 || g_cfg.skew_num < 0 || g_cfg.skew_depth < 0 || g_cfg.skew_min_size < 0 ||
        (g_cfg.endgame_mode != 0 && g_cfg.endgame_mode != 1)) {
        usage(argv[0]);
        return 1;
    }

    long terms = g_cfg.digits / 14 + 10;
    mpfr_prec_t prec = (mpfr_prec_t)(g_cfg.digits * 3.32192809488736234787) + 256;

    mpz_t P, Q, T;
    mpz_inits(P, Q, T, NULL);

    double t0 = now_seconds();

    #pragma omp parallel
    {
        #pragma omp single nowait
        bs_parallel(P, Q, T, 0, terms, 0);
    }

    double t1 = now_seconds();

    mpfr_t pi;
    mpfr_init2(pi, prec);
    if (g_cfg.endgame_mode == 0) {
        final_eval_direct(pi, Q, T, prec);
    } else {
        final_eval_newton(pi, Q, T, prec);
    }

    double t2 = now_seconds();

    printf("method=unified_exact_bs\n");
    printf("digits=%ld\n", g_cfg.digits);
    printf("terms=%ld\n", terms);
    printf("task_depth=%d\n", g_cfg.task_depth);
    printf("leaf_size=%ld\n", g_cfg.leaf_size);
    printf("task_min_size=%ld\n", g_cfg.task_min_size);
    printf("skew=%d/%d\n", g_cfg.skew_num, g_cfg.skew_den);
    printf("skew_depth=%d\n", g_cfg.skew_depth);
    printf("skew_min_size=%ld\n", g_cfg.skew_min_size);
    printf("endgame_mode=%d\n", g_cfg.endgame_mode);
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
    mpfr_clear(pi);
    mpfr_free_cache();
    return 0;
}

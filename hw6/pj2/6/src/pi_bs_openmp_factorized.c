#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define A 13591409UL
#define B 545140134UL
#define C3_OVER_24 10939058860032000UL

static long g_leaf_size = 128;
static long g_task_min_size = 4096;
static int g_cancel_level = 4;
static long g_cancel_min_size = 8192;

static uint32_t *g_spf_odd = NULL;
static uint32_t g_spf_max_n = 0;
static unsigned long long g_cancel_count = 0;
static unsigned long long g_cancel_prime_hits = 0;

typedef struct {
    uint32_t *prime;
    uint32_t *exp;
    size_t size;
    size_t cap;
} factor_list_t;

typedef struct {
    mpz_t P, Q, T;
    factor_list_t fp, fq;
} bs_node_t;

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static void fl_init(factor_list_t *fl) {
    fl->prime = NULL;
    fl->exp = NULL;
    fl->size = 0;
    fl->cap = 0;
}

static void fl_free(factor_list_t *fl) {
    free(fl->prime);
    free(fl->exp);
    fl->prime = NULL;
    fl->exp = NULL;
    fl->size = 0;
    fl->cap = 0;
}

static void fl_clear(factor_list_t *fl) {
    fl->size = 0;
}

static void fl_reserve(factor_list_t *fl, size_t need) {
    if (need <= fl->cap) return;
    size_t new_cap = fl->cap ? fl->cap : 8;
    while (new_cap < need) new_cap *= 2;
    fl->prime = (uint32_t *)realloc(fl->prime, new_cap * sizeof(uint32_t));
    fl->exp = (uint32_t *)realloc(fl->exp, new_cap * sizeof(uint32_t));
    if (!fl->prime || !fl->exp) {
        fprintf(stderr, "Out of memory in fl_reserve\n");
        exit(2);
    }
    fl->cap = new_cap;
}

static void fl_push_back(factor_list_t *fl, uint32_t p, uint32_t e) {
    if (e == 0) return;
    if (fl->size > 0 && fl->prime[fl->size - 1] == p) {
        fl->exp[fl->size - 1] += e;
        return;
    }
    fl_reserve(fl, fl->size + 1);
    fl->prime[fl->size] = p;
    fl->exp[fl->size] = e;
    fl->size++;
}

static void fl_copy(factor_list_t *dst, const factor_list_t *src) {
    fl_reserve(dst, src->size);
    if (src->size) {
        memcpy(dst->prime, src->prime, src->size * sizeof(uint32_t));
        memcpy(dst->exp, src->exp, src->size * sizeof(uint32_t));
    }
    dst->size = src->size;
}

static void fl_merge_sum_into(factor_list_t *dst, const factor_list_t *a, const factor_list_t *b) {
    size_t ia = 0, ib = 0, k = 0;
    fl_reserve(dst, a->size + b->size);
    while (ia < a->size && ib < b->size) {
        if (a->prime[ia] == b->prime[ib]) {
            dst->prime[k] = a->prime[ia];
            dst->exp[k] = a->exp[ia] + b->exp[ib];
            ia++; ib++; k++;
        } else if (a->prime[ia] < b->prime[ib]) {
            dst->prime[k] = a->prime[ia];
            dst->exp[k] = a->exp[ia];
            ia++; k++;
        } else {
            dst->prime[k] = b->prime[ib];
            dst->exp[k] = b->exp[ib];
            ib++; k++;
        }
    }
    while (ia < a->size) {
        dst->prime[k] = a->prime[ia];
        dst->exp[k] = a->exp[ia];
        ia++; k++;
    }
    while (ib < b->size) {
        dst->prime[k] = b->prime[ib];
        dst->exp[k] = b->exp[ib];
        ib++; k++;
    }
    dst->size = k;
}

static void fl_compact(factor_list_t *fl) {
    size_t j = 0;
    for (size_t i = 0; i < fl->size; ++i) {
        if (fl->exp[i] != 0) {
            if (j != i) {
                fl->prime[j] = fl->prime[i];
                fl->exp[j] = fl->exp[i];
            }
            j++;
        }
    }
    fl->size = j;
}

static void bs_node_init(bs_node_t *node) {
    mpz_inits(node->P, node->Q, node->T, NULL);
    fl_init(&node->fp);
    fl_init(&node->fq);
}

static void bs_node_clear(bs_node_t *node) {
    mpz_clears(node->P, node->Q, node->T, NULL);
    fl_free(&node->fp);
    fl_free(&node->fq);
}

static void build_spf_odd(uint32_t max_n) {
    if (max_n <= g_spf_max_n && g_spf_odd != NULL) return;
    free(g_spf_odd);
    g_spf_max_n = max_n;
    size_t len = ((size_t)max_n >> 1) + 1;
    g_spf_odd = (uint32_t *)calloc(len, sizeof(uint32_t));
    if (!g_spf_odd) {
        fprintf(stderr, "Failed to allocate SPF table\n");
        exit(2);
    }
    uint32_t limit = (uint32_t)(sqrt((double)max_n) + 1.0);
    for (uint32_t p = 3; p <= limit; p += 2) {
        size_t idx = (size_t)p >> 1;
        if (g_spf_odd[idx] != 0) continue;
        uint64_t start = (uint64_t)p * (uint64_t)p;
        for (uint64_t x = start; x <= max_n; x += 2ULL * p) {
            size_t xi = (size_t)(x >> 1);
            if (g_spf_odd[xi] == 0) g_spf_odd[xi] = p;
        }
    }
}

static void factorize_uint_sorted_append(factor_list_t *fl, uint32_t n, uint32_t mult) {
    if (n == 0 || mult == 0) return;
    if ((n & 1U) == 0U) {
        uint32_t e = 0;
        while ((n & 1U) == 0U) {
            n >>= 1;
            e++;
        }
        fl_push_back(fl, 2U, e * mult);
    }
    while (n > 1U) {
        uint32_t p = g_spf_odd[n >> 1];
        if (p == 0) p = n;
        uint32_t e = 0;
        do {
            n /= p;
            e++;
        } while (n > 1U && n % p == 0U);
        fl_push_back(fl, p, e * mult);
    }
}

static void build_term_factors_p(factor_list_t *dst, uint32_t k) {
    fl_clear(dst);
    factorize_uint_sorted_append(dst, 2U * k - 1U, 1U);
    factorize_uint_sorted_append(dst, 6U * k - 5U, 1U);
    factorize_uint_sorted_append(dst, 6U * k - 1U, 1U);
}

static void build_term_factors_q(factor_list_t *dst, uint32_t k) {
    fl_clear(dst);
    factorize_uint_sorted_append(dst, k, 3U);
    fl_push_back(dst, 2U, 15U);
    fl_push_back(dst, 3U, 2U);
    fl_push_back(dst, 5U, 3U);
    fl_push_back(dst, 23U, 3U);
    fl_push_back(dst, 29U, 3U);
}

static void common_factor_build_and_cancel(mpz_t leftP, factor_list_t *leftFp,
                                           mpz_t rightQ, factor_list_t *rightFq) {
    factor_list_t common;
    fl_init(&common);
    fl_reserve(&common, (leftFp->size < rightFq->size) ? leftFp->size : rightFq->size);

    size_t i = 0, j = 0;
    while (i < leftFp->size && j < rightFq->size) {
        uint32_t pa = leftFp->prime[i];
        uint32_t pb = rightFq->prime[j];
        if (pa == pb) {
            uint32_t ea = leftFp->exp[i];
            uint32_t eb = rightFq->exp[j];
            uint32_t m = ea < eb ? ea : eb;
            if (m) {
                fl_push_back(&common, pa, m);
                leftFp->exp[i] -= m;
                rightFq->exp[j] -= m;
            }
            i++; j++;
        } else if (pa < pb) {
            i++;
        } else {
            j++;
        }
    }

    if (common.size != 0) {
        mpz_t g, pe;
        mpz_inits(g, pe, NULL);
        mpz_set_ui(g, 1U);
        for (size_t k = 0; k < common.size; ++k) {
            mpz_ui_pow_ui(pe, common.prime[k], common.exp[k]);
            mpz_mul(g, g, pe);
        }
        mpz_divexact(leftP, leftP, g);
        mpz_divexact(rightQ, rightQ, g);
        mpz_clears(g, pe, NULL);
        fl_compact(leftFp);
        fl_compact(rightFq);
#ifdef _OPENMP
        #pragma omp atomic update
#endif
        g_cancel_count += 1ULL;
#ifdef _OPENMP
        #pragma omp atomic update
#endif
        g_cancel_prime_hits += (unsigned long long)common.size;
    }

    fl_free(&common);
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

static void bs_leaf(bs_node_t *node, long a, long b) {
    mpz_t p, q, t, tmp;
    mpz_inits(p, q, t, tmp, NULL);
    factor_list_t term_fp, term_fq, tmp_fl;
    fl_init(&term_fp);
    fl_init(&term_fq);
    fl_init(&tmp_fl);

    mpz_set_ui(node->P, 1U);
    mpz_set_ui(node->Q, 1U);
    mpz_set_ui(node->T, 0U);
    fl_clear(&node->fp);
    fl_clear(&node->fq);

    for (long kk = b - 1; kk >= a; --kk) {
        if (kk == 0) {
            mpz_addmul_ui(node->T, node->Q, A);
            continue;
        }
        uint32_t k = (uint32_t)kk;
        term_pqt(p, q, t, k);
        build_term_factors_p(&term_fp, k);
        build_term_factors_q(&term_fq, k);

        mpz_mul(tmp, p, node->T);
        mpz_mul(node->T, t, node->Q);
        mpz_add(node->T, node->T, tmp);

        mpz_mul(node->P, node->P, p);
        mpz_mul(node->Q, node->Q, q);

        fl_merge_sum_into(&tmp_fl, &node->fp, &term_fp);
        fl_copy(&node->fp, &tmp_fl);
        fl_merge_sum_into(&tmp_fl, &node->fq, &term_fq);
        fl_copy(&node->fq, &tmp_fl);
    }

    fl_free(&term_fp);
    fl_free(&term_fq);
    fl_free(&tmp_fl);
    mpz_clears(p, q, t, tmp, NULL);
}

static void combine_nodes(bs_node_t *out, bs_node_t *left, bs_node_t *right, int level, long interval_size) {
    if (g_cancel_level >= 0 && level >= g_cancel_level && interval_size >= g_cancel_min_size) {
        common_factor_build_and_cancel(left->P, &left->fp, right->Q, &right->fq);
    }

    mpz_mul(out->P, left->P, right->P);
    mpz_mul(out->Q, left->Q, right->Q);
    mpz_mul(out->T, left->P, right->T);
    mpz_addmul(out->T, left->T, right->Q);

    fl_merge_sum_into(&out->fp, &left->fp, &right->fp);
    fl_merge_sum_into(&out->fq, &left->fq, &right->fq);
}

static void bs_serial(bs_node_t *out, long a, long b, int level) {
    long n = b - a;
    if (n <= g_leaf_size) {
        bs_leaf(out, a, b);
        return;
    }
    long m = a + n / 2;
    bs_node_t left, right;
    bs_node_init(&left);
    bs_node_init(&right);

    bs_serial(&left, a, m, level + 1);
    bs_serial(&right, m, b, level + 1);
    combine_nodes(out, &left, &right, level, n);

    bs_node_clear(&left);
    bs_node_clear(&right);
}

static void bs_parallel(bs_node_t *out, long a, long b, int level, int max_depth) {
    long n = b - a;
    if (n <= g_leaf_size) {
        bs_leaf(out, a, b);
        return;
    }
    if (level >= max_depth || n < g_task_min_size) {
        bs_serial(out, a, b, level);
        return;
    }
    long m = a + n / 2;
    bs_node_t left, right;
    bs_node_init(&left);
    bs_node_init(&right);

    #pragma omp task shared(left)
    bs_parallel(&left, a, m, level + 1, max_depth);

    #pragma omp task shared(right)
    bs_parallel(&right, m, b, level + 1, max_depth);

    #pragma omp taskwait

    combine_nodes(out, &left, &right, level, n);

    bs_node_clear(&left);
    bs_node_clear(&right);
}

int main(int argc, char **argv) {
    long digits = 10000000;
    int max_depth = 5;
    if (argc >= 2) digits = atol(argv[1]);
    if (argc >= 3) max_depth = atoi(argv[2]);
    if (argc >= 4) g_leaf_size = atol(argv[3]);
    if (argc >= 5) g_task_min_size = atol(argv[4]);
    if (argc >= 6) g_cancel_level = atoi(argv[5]);
    if (argc >= 7) g_cancel_min_size = atol(argv[6]);
    if (digits <= 0 || max_depth < 0 || g_leaf_size <= 0 || g_task_min_size <= 0 || g_cancel_min_size < 0) {
        fprintf(stderr, "Usage: %s [digits] [task_depth] [leaf_size] [task_min_size] [cancel_level] [cancel_min_size]\n", argv[0]);
        return 1;
    }

    long terms = digits / 14 + 10;
    uint32_t sieve_max = (uint32_t)(6L * terms + 16L);
    build_spf_odd(sieve_max);

    mpfr_prec_t prec = (mpfr_prec_t)(digits * 3.32192809488736234787) + 256;
    bs_node_t root;
    bs_node_init(&root);

    double t0 = now_seconds();
    #pragma omp parallel
    {
        #pragma omp single nowait
        bs_parallel(&root, 0, terms, 0, max_depth);
    }
    double t1 = now_seconds();

    mpfr_t pi, sqrt10005, factor, qf, tf;
    mpfr_inits2(prec, pi, sqrt10005, factor, qf, tf, (mpfr_ptr)0);
    mpfr_set_ui(sqrt10005, 10005U, MPFR_RNDN);
    mpfr_sqrt(sqrt10005, sqrt10005, MPFR_RNDN);
    mpfr_mul_ui(factor, sqrt10005, 426880U, MPFR_RNDN);
    mpfr_set_z(qf, root.Q, MPFR_RNDN);
    mpfr_set_z(tf, root.T, MPFR_RNDN);
    mpfr_mul(pi, factor, qf, MPFR_RNDN);
    mpfr_div(pi, pi, tf, MPFR_RNDN);
    double t2 = now_seconds();

    printf("method=factorized_openmp_exact_bs\n");
    printf("digits=%ld\n", digits);
    printf("terms=%ld\n", terms);
    printf("task_depth=%d\n", max_depth);
    printf("leaf_size=%ld\n", g_leaf_size);
    printf("task_min_size=%ld\n", g_task_min_size);
    printf("cancel_level=%d\n", g_cancel_level);
    printf("cancel_min_size=%ld\n", g_cancel_min_size);
    printf("sieve_max=%u\n", sieve_max);
#ifdef _OPENMP
    printf("omp_threads=%d\n", omp_get_max_threads());
#else
    printf("omp_threads=1\n");
#endif
    printf("cancel_count=%llu\n", g_cancel_count);
    printf("cancel_prime_hits=%llu\n", g_cancel_prime_hits);
    printf("binary_splitting_time=%.6f s\n", t1 - t0);
    printf("final_eval_time=%.6f s\n", t2 - t1);
    printf("total_time=%.6f s\n", t2 - t0);
    mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);

    bs_node_clear(&root);
    free(g_spf_odd);
    g_spf_odd = NULL;
    mpfr_clears(pi, sqrt10005, factor, qf, tf, (mpfr_ptr)0);
    mpfr_free_cache();
    return 0;
}

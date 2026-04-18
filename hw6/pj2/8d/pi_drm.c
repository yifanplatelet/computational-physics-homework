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

// 定义稀疏素数因子结构
typedef struct {
    int p;
    int e;
} fac_t;

typedef struct {
    fac_t *arr;
    int size;
    int cap;
} fac_list_t;

int *spf; // 最小素因子数组 (Sieve)

// 初始化 Sieve 筛法
void init_sieve(long max_val) {
    spf = malloc((max_val + 1) * sizeof(int));
    for (long i = 0; i <= max_val; i++) spf[i] = i;
    for (long i = 2; i * i <= max_val; i++) {
        if (spf[i] == i) {
            for (long j = i * i; j <= max_val; j += i) {
                if (spf[j] == j) spf[j] = i;
            }
        }
    }
}

void fac_init(fac_list_t *fl, int cap) {
    fl->arr = malloc(cap * sizeof(fac_t));
    fl->size = 0;
    fl->cap = cap;
}

void fac_free(fac_list_t *fl) {
    free(fl->arr);
}

void fac_push(fac_list_t *fl, int p, int e) {
    if (fl->size >= fl->cap) {
        fl->cap *= 2;
        fl->arr = realloc(fl->arr, fl->cap * sizeof(fac_t));
    }
    fl->arr[fl->size++] = (fac_t){p, e};
}

// 质因数分解，结果天然有序
void factorize(long n, fac_list_t *fl) {
    fac_init(fl, 10);
    while (n > 1) {
        int p = spf[n];
        int e = 0;
        while (n % p == 0) { e++; n /= p; }
        fac_push(fl, p, e);
    }
}

// 获取 Q 的常数部分因子 (C^3 / 24)
void get_C3_24_factors(fac_list_t *fl) {
    fac_init(fl, 5);
    fac_push(fl, 2, 15);
    fac_push(fl, 3, 2);
    fac_push(fl, 5, 3);
    fac_push(fl, 23, 3);
    fac_push(fl, 29, 3);
}

// 两个有序因子列表相乘 (指数相加)
void fac_mul(fac_list_t *out, fac_list_t *a, fac_list_t *b) {
    fac_init(out, a->size + b->size);
    int i = 0, j = 0;
    while (i < a->size && j < b->size) {
        if (a->arr[i].p < b->arr[j].p) {
            fac_push(out, a->arr[i].p, a->arr[i].e); i++;
        } else if (a->arr[i].p > b->arr[j].p) {
            fac_push(out, b->arr[j].p, b->arr[j].e); j++;
        } else {
            fac_push(out, a->arr[i].p, a->arr[i].e + b->arr[j].e);
            i++; j++;
        }
    }
    while (i < a->size) { fac_push(out, a->arr[i].p, a->arr[i].e); i++; }
    while (j < b->size) { fac_push(out, b->arr[j].p, b->arr[j].e); j++; }
}

// 核心：提取公共素数因子并剔除 (GCD Extraction)
void fac_remove_gcd(fac_list_t *p1, fac_list_t *q2, fac_list_t *out_p1, fac_list_t *out_q2) {
    fac_init(out_p1, p1->size);
    fac_init(out_q2, q2->size);
    int i = 0, j = 0;
    while (i < p1->size && j < q2->size) {
        if (p1->arr[i].p < q2->arr[j].p) {
            fac_push(out_p1, p1->arr[i].p, p1->arr[i].e); i++;
        } else if (p1->arr[i].p > q2->arr[j].p) {
            fac_push(out_q2, q2->arr[j].p, q2->arr[j].e); j++;
        } else {
            int min_e = p1->arr[i].e < q2->arr[j].e ? p1->arr[i].e : q2->arr[j].e;
            if (p1->arr[i].e > min_e) fac_push(out_p1, p1->arr[i].p, p1->arr[i].e - min_e);
            if (q2->arr[j].e > min_e) fac_push(out_q2, q2->arr[j].p, q2->arr[j].e - min_e);
            i++; j++;
        }
    }
    while (i < p1->size) { fac_push(out_p1, p1->arr[i].p, p1->arr[i].e); i++; }
    while (j < q2->size) { fac_push(out_q2, q2->arr[j].p, q2->arr[j].e); j++; }
}

// 快速素数连乘树 (将稀疏数组转化为 GMP 整数)
void build_product_tree(mpz_t res, fac_list_t *fl) {
    if (fl->size == 0) { mpz_set_ui(res, 1); return; }
    mpz_t *tree = malloc(fl->size * sizeof(mpz_t));
    for (int i = 0; i < fl->size; i++) {
        mpz_init(tree[i]);
        mpz_ui_pow_ui(tree[i], fl->arr[i].p, fl->arr[i].e);
    }
    int step = 1;
    while (step < fl->size) {
        for (int i = 0; i + step < fl->size; i += 2 * step) {
            mpz_mul(tree[i], tree[i], tree[i + step]);
        }
        step *= 2;
    }
    mpz_set(res, tree[0]);
    for (int i = 0; i < fl->size; i++) mpz_clear(tree[i]);
    free(tree);
}

// 叶子节点计算
void bs_leaf(long a, fac_list_t *P, fac_list_t *Q, mpz_t T) {
    if (a == 0) {
        fac_init(P, 1); fac_init(Q, 1);
        mpz_set_ui(T, A_CONST);
        return;
    }
    fac_list_t f1, f2, f3, p_tmp;
    factorize(6*a-5, &f1); factorize(2*a-1, &f2); factorize(6*a-1, &f3);
    fac_mul(&p_tmp, &f1, &f2); fac_mul(P, &p_tmp, &f3);
    fac_free(&f1); fac_free(&f2); fac_free(&f3); fac_free(&p_tmp);

    fac_list_t fk, fk2, fk3, fc;
    factorize(a, &fk); fac_mul(&fk2, &fk, &fk); fac_mul(&fk3, &fk2, &fk);
    get_C3_24_factors(&fc); fac_mul(Q, &fk3, &fc);
    fac_free(&fk); fac_free(&fk2); fac_free(&fk3); fac_free(&fc);

    mpz_t mpz_P; mpz_init(mpz_P);
    build_product_tree(mpz_P, P);
    unsigned long linear = A_CONST + B_CONST * a;
    mpz_mul_ui(T, mpz_P, linear);
    if (a & 1) mpz_neg(T, T);
    mpz_clear(mpz_P);
}

// 递归分解与提取 (混合多线程)
void bs_parallel(long a, long b, fac_list_t *P, fac_list_t *Q, mpz_t T, int depth) {
    if (b - a == 1) {
        bs_leaf(a, P, Q, T);
        return;
    }
    long m = a + (b - a) / 2;
    fac_list_t P1, Q1, P2, Q2;
    mpz_t T1, T2; mpz_init(T1); mpz_init(T2);

    if (depth < 6) {
        #pragma omp task shared(P1, Q1, T1)
        bs_parallel(a, m, &P1, &Q1, T1, depth + 1);
        bs_parallel(m, b, &P2, &Q2, T2, depth + 1);
        #pragma omp taskwait
    } else {
        bs_parallel(a, m, &P1, &Q1, T1, depth + 1);
        bs_parallel(m, b, &P2, &Q2, T2, depth + 1);
    }

    fac_mul(P, &P1, &P2);
    fac_mul(Q, &Q1, &Q2);

    // 绝杀环节：去除 P1 和 Q2 的公共因子
    fac_list_t P1_opt, Q2_opt;
    fac_remove_gcd(&P1, &Q2, &P1_opt, &Q2_opt);

    mpz_t mpz_P1, mpz_Q2;
    mpz_init(mpz_P1); mpz_init(mpz_Q2);
    
    // 如果是在浅层，将多线程并发推到极限
    if (depth < 6) {
        #pragma omp task shared(mpz_P1, P1_opt)
        build_product_tree(mpz_P1, &P1_opt);
        #pragma omp task shared(mpz_Q2, Q2_opt)
        build_product_tree(mpz_Q2, &Q2_opt);
        #pragma omp taskwait

        mpz_t tmp_T1, tmp_T2;
        mpz_init(tmp_T1); mpz_init(tmp_T2);

        #pragma omp task shared(tmp_T1, T1, mpz_Q2)
        mpz_mul(tmp_T1, T1, mpz_Q2);
        #pragma omp task shared(tmp_T2, T2, mpz_P1)
        mpz_mul(tmp_T2, T2, mpz_P1);
        #pragma omp taskwait

        mpz_add(T, tmp_T1, tmp_T2);
        mpz_clears(tmp_T1, tmp_T2, NULL);
    } else {
        build_product_tree(mpz_P1, &P1_opt);
        build_product_tree(mpz_Q2, &Q2_opt);
        mpz_mul(T1, T1, mpz_Q2);
        mpz_mul(T2, T2, mpz_P1);
        mpz_add(T, T1, T2);
    }

    mpz_clears(T1, T2, mpz_P1, mpz_Q2, NULL);
    fac_free(&P1); fac_free(&Q1); fac_free(&P2); fac_free(&Q2);
    fac_free(&P1_opt); fac_free(&Q2_opt);
}

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

int main(int argc, char **argv) {
    long digits = 100000000;
    if (argc >= 2) digits = atol(argv[1]);

    long terms = digits / 14 + 10;
    
    // 初始化筛法 (提取最大因数界限)
    printf("Initializing Prime Sieve...\n");
    init_sieve(6 * terms);

    mpz_t T; mpz_init(T);
    fac_list_t root_P, root_Q;

    double t0 = now_seconds();

    #pragma omp parallel
    {
        #pragma omp single nowait
        bs_parallel(0, terms, &root_P, &root_Q, T, 0);
    }

    // 最后才将 Q 的素数列表还原为 GMP 整数
    mpz_t Q; mpz_init(Q);
    build_product_tree(Q, &root_Q);

    double t1 = now_seconds();

    mpfr_prec_t prec = (mpfr_prec_t)(digits * 3.321928094887) + 256;
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

    printf("Algorithm: DRM (Prime Factorization + GCD Extraction)\n");
    printf("Digits: %ld\n", digits);
    printf("Terms: %ld\n", terms);
    printf("Binary Splitting Time: %.6f s\n", t1 - t0);
    printf("Final Eval Time: %.6f s\n", t2 - t1);
    printf("Total Time: %.6f s\n", t2 - t0);
    mpfr_printf("pi(first 50 digits)=%.50RNf\n", pi);

    mpz_clears(T, Q, NULL);
    fac_free(&root_P); fac_free(&root_Q); free(spf);
    mpfr_clears(pi, sqrt10005, factor, qf, tf, (mpfr_ptr)0);
    mpfr_free_cache();
    return 0;
}
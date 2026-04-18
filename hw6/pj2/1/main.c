#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>

#define DIGITS 100000

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

int main(void) {
    // 十进制位数 -> 二进制精度，额外留保护位
    mpfr_prec_t prec = (mpfr_prec_t)(DIGITS * 3.3219280948873626) + 128;

    mpfr_t sum, term, pi_calc, pi_ref, diff, sqrt10005, C;
    mpfr_t numer_f, denom_f;
    mpz_t sixk_fact, threek_fact, k_fact, k_fact_cubed;
    mpz_t numer_z, denom_z, pow_z, tmp_z, linear_z;

    mpfr_inits2(prec, sum, term, pi_calc, pi_ref, diff, sqrt10005, C, numer_f, denom_f, (mpfr_ptr) 0);
    mpz_inits(sixk_fact, threek_fact, k_fact, k_fact_cubed, numer_z, denom_z, pow_z, tmp_z, linear_z, (mpz_ptr) 0);

    mpfr_set_zero(sum, 0);

    // 每项约增加 14 位，略多取几项
    int terms = DIGITS / 14 + 5;

    double t0 = now_seconds();

    for (int k = 0; k < terms; ++k) {
        // (6k)!, (3k)!, k!
        mpz_fac_ui(sixk_fact, 6u * (unsigned long)k);
        mpz_fac_ui(threek_fact, 3u * (unsigned long)k);
        mpz_fac_ui(k_fact, (unsigned long)k);

        // (k!)^3
        mpz_mul(k_fact_cubed, k_fact, k_fact);
        mpz_mul(k_fact_cubed, k_fact_cubed, k_fact);

        // linear_z = 13591409 + 545140134*k
        mpz_set_ui(linear_z, 545140134u);
        mpz_mul_ui(linear_z, linear_z, (unsigned long)k);
        mpz_add_ui(linear_z, linear_z, 13591409u);

        // numer_z = (6k)! * (13591409 + 545140134k)
        mpz_mul(numer_z, sixk_fact, linear_z);

        // denom_z = (3k)! * (k!)^3 * 640320^(3k)
        mpz_ui_pow_ui(pow_z, 640320u, 3u * (unsigned long)k);
        mpz_mul(denom_z, threek_fact, k_fact_cubed);
        mpz_mul(denom_z, denom_z, pow_z);

        // 转成 mpfr 做高精度除法
        mpfr_set_z(numer_f, numer_z, MPFR_RNDN);
        mpfr_set_z(denom_f, denom_z, MPFR_RNDN);
        mpfr_div(term, numer_f, denom_f, MPFR_RNDN);

        if (k & 1) {
            mpfr_neg(term, term, MPFR_RNDN);
        }

        mpfr_add(sum, sum, term, MPFR_RNDN);
    }

    // C = 426880 * sqrt(10005)
    mpfr_set_ui(sqrt10005, 10005u, MPFR_RNDN);
    mpfr_sqrt(sqrt10005, sqrt10005, MPFR_RNDN);
    mpfr_set_ui(C, 426880u, MPFR_RNDN);
    mpfr_mul(C, C, sqrt10005, MPFR_RNDN);

    // pi_calc = C / sum
    mpfr_div(pi_calc, C, sum, MPFR_RNDN);

    double t1 = now_seconds();

    // MPFR 直接给出高精度 pi 作为参考值
    mpfr_const_pi(pi_ref, MPFR_RNDN);

    // diff = |pi_calc - pi_ref|
    mpfr_sub(diff, pi_calc, pi_ref, MPFR_RNDN);
    mpfr_abs(diff, diff, MPFR_RNDN);

    printf("目标小数位数: %d\n", DIGITS);
    printf("Chudnovsky 项数: %d\n", terms);
    printf("计算耗时: %.6f 秒\n", t1 - t0);

    // 输出前 50 位看看
    mpfr_printf("pi_calc (前50位) = %.50RNf\n", pi_calc);
    mpfr_printf("pi_ref  (前50位) = %.50RNf\n", pi_ref);

    // 误差输出为科学计数法
    mpfr_printf("|pi_calc - pi_ref| = %.20Re\n", diff);

    // 还可以估计误差对应的大致正确位数
    if (mpfr_zero_p(diff)) {
        printf("误差为 0（在当前设置精度下完全一致）\n");
    } else {
        mpfr_t log10diff;
        mpfr_init2(log10diff, prec);
        mpfr_log10(log10diff, diff, MPFR_RNDN);
        mpfr_printf("log10(误差) = %.10Rf\n", log10diff);
        mpfr_clear(log10diff);
    }

    mpfr_clears(sum, term, pi_calc, pi_ref, diff, sqrt10005, C, numer_f, denom_f, (mpfr_ptr) 0);
    mpz_clears(sixk_fact, threek_fact, k_fact, k_fact_cubed, numer_z, denom_z, pow_z, tmp_z, linear_z, (mpz_ptr) 0);

    mpfr_free_cache();
    return 0;
}
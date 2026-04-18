#include <stdio.h>
#include <time.h>
#include <math.h>

/* Gregory-Leibniz 正序求和 */
long double pi_forward_gl(long long n) {
    long double sum = 0.0L;
    long long k;
    for (k = 1; k <= n; k++) {
        long double term = 1.0L / (2.0L * k - 1.0L);
        if (k % 2 == 0)
            sum -= term;
        else
            sum += term;
    }
    return 4.0L * sum;
}

/* Gregory-Leibniz 倒序求和 */
long double pi_backward_gl(long long n) {
    long double sum = 0.0L;
    long long k;
    for (k = n; k >= 1; k--) {
        long double term = 1.0L / (2.0L * k - 1.0L);
        if (k % 2 == 0)
            sum -= term;
        else
            sum += term;
    }
    return 4.0L * sum;
}

int main(void) {
    /* 用一个较小的 n 来测试时间 */
    long long test_n = 5000000LL;              /* 500万，可自行改成 1000000 或 10000000 */
    long long target_n = 40000000000LL;        /* 约 4e10，达到 10 位精度所需项数 */

    long double pi_f, pi_b, pi_true;
    double time_f, time_b;
    double est_time_f, est_time_b;

    clock_t start, end;

    pi_true = acosl(-1.0L);

    printf("Gregory-Leibniz series for pi\n");
    printf("============================================\n");
    printf("Test n   = %lld\n", test_n);
    printf("Target n = %lld (estimated for 10-digit accuracy)\n\n", target_n);

    /* 正序计时 */
    start = clock();
    pi_f = pi_forward_gl(test_n);
    end = clock();
    time_f = (double)(end - start) / CLOCKS_PER_SEC;

    /* 倒序计时 */
    start = clock();
    pi_b = pi_backward_gl(test_n);
    end = clock();
    time_b = (double)(end - start) / CLOCKS_PER_SEC;

    /* 按线性比例估算 4e10 项所需时间 */
    est_time_f = time_f * (double)target_n / (double)test_n;
    est_time_b = time_b * (double)target_n / (double)test_n;

    printf("Forward result  = %.15Lf\n", pi_f);
    printf("Backward result = %.15Lf\n", pi_b);
    printf("True pi         = %.15Lf\n\n", pi_true);

    printf("Forward error   = %.15Le\n", fabsl(pi_f - pi_true));
    printf("Backward error  = %.15Le\n\n", fabsl(pi_b - pi_true));

    printf("Measured time with n = %lld:\n", test_n);
    printf("  Forward time  = %.6f seconds\n", time_f);
    printf("  Backward time = %.6f seconds\n\n", time_b);

    printf("Estimated time for n = %lld:\n", target_n);
    printf("  Forward time  = %.2f seconds (%.2f minutes, %.2f hours)\n",
           est_time_f, est_time_f / 60.0, est_time_f / 3600.0);
    printf("  Backward time = %.2f seconds (%.2f minutes, %.2f hours)\n\n",
           est_time_b, est_time_b / 60.0, est_time_b / 3600.0);

    printf("Conclusion:\n");
    printf("1. Backward summation is usually slightly more accurate than forward summation\n");
    printf("   because it reduces floating-point rounding error.\n");
    printf("2. However, the Gregory-Leibniz series converges very slowly.\n");
    printf("3. To obtain about 10 decimal digits of accuracy, we need roughly n ~ 4e10 terms.\n");
    printf("4. Based on the measured time for a smaller n, directly computing 4e10 terms\n");
    printf("   would take a very long time, so it is generally impractical on an ordinary computer.\n");

    return 0;
}
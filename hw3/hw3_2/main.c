#include <stdio.h>
#include <quadmath.h>

// 用反正切级数计算 arctan(x)
__float128 arctan_series(__float128 x)
{
    __float128 sum = 0.0Q;
    __float128 term = x;
    __float128 x2 = x * x;
    int n = 0;

    while (fabsq(term) > 1e-40Q) {
        sum += term / (2 * n + 1);
        term *= -x2;
        n++;
    }

    return sum;
}

int main()
{
    // Machin公式计算pi
    __float128 a = arctan_series(1.0Q / 5.0Q);
    __float128 b = arctan_series(1.0Q / 239.0Q);
    __float128 pi_calc = 4.0Q * (4.0Q * a - b);

    // pi真值：利用 acos(-1)
    __float128 pi_true = acosq(-1.0Q);

    // 误差
    __float128 abs_error = fabsq(pi_calc - pi_true);

    char buf1[128], buf2[128], buf3[128];

    quadmath_snprintf(buf1, sizeof(buf1), "%.30Qf", pi_calc);
    quadmath_snprintf(buf2, sizeof(buf2), "%.30Qf", pi_true);
    quadmath_snprintf(buf3, sizeof(buf3), "%.40Qe", abs_error);

    printf("Machin公式计算得到的 pi = %s\n", buf1);
    printf("pi 真值(acos(-1))    = %s\n", buf2);
    printf("绝对误差             = %s\n", buf3);

    return 0;
}
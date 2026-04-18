#include <stdio.h>
#include <math.h>
#include <float.h>

int main(void) {
    /* IEEE 754 binary16 理论参数 */
    const int p = 11;        // 有效二进制位数（含隐藏位）
    const int emin = -14;    // 最小正规数指数
    const int emax = 15;     // 最大正规数指数

    double eps = ldexp(1.0, 1 - p);                 // 2^(1-p) = 2^-10
    double u   = ldexp(1.0, -p);                    // 2^-11
    double xmin_normal = ldexp(1.0, emin);          // 2^-14
    double xmin_sub    = ldexp(1.0, emin - (p - 1)); // 2^-24
    double xmax = (2.0 - ldexp(1.0, 1 - p)) * ldexp(1.0, emax);

    printf("IEEE 754 half precision (binary16)\n");
    printf("---------------------------------\n");
    printf("effective precision p      = %d bits\n", p);
    printf("machine epsilon eps        = 2^(1-p) = %.17g\n", eps);
    printf("unit roundoff u            = 2^(-p)  = %.17g\n", u);
    printf("min positive normal        = %.17g\n", xmin_normal);
    printf("min positive subnormal     = %.17g\n", xmin_sub);
    printf("max finite number          = %.17g\n", xmax);
    printf("approx decimal digits      = %.4f\n", p * log10(2.0));

#if defined(FLT16_MANT_DIG)
    printf("\nCompiler reports _Float16 support:\n");
    printf("FLT16_EPSILON              = %.17g\n", (double)FLT16_EPSILON);
    printf("FLT16_MIN                  = %.17g\n", (double)FLT16_MIN);
    printf("FLT16_TRUE_MIN             = %.17g\n", (double)FLT16_TRUE_MIN);
    printf("FLT16_MAX                  = %.17g\n", (double)FLT16_MAX);

    /* 用 _Float16 实际演示舍入误差 */
    _Float16 a = (_Float16)10000.0;
    _Float16 b = (_Float16)1.0;
    _Float16 c = (_Float16)((a + b) - a);

    printf("\nRoundoff example:\n");
    printf("half   : (10000 + 1) - 10000 = %.17g\n", (double)c);
    printf("double : (10000 + 1) - 10000 = %.17g\n", (10000.0 + 1.0) - 10000.0);

    /* 再看 1 附近的间隔 */
    _Float16 one = (_Float16)1.0;
    _Float16 half_eps = (_Float16)(eps / 2.0);

    printf("\nCheck spacing around 1:\n");
    printf("1 + eps/2 in half = %.17g\n", (double)(_Float16)(one + half_eps));
    printf("1 + eps   in half = %.17g\n", (double)(_Float16)(one + (_Float16)eps));

#else
    printf("\nThis compiler does not provide _Float16 / FLT16_* macros.\n");
    printf("But the theoretical values above are the exact IEEE 754 binary16 values.\n");
#endif

    return 0;
}
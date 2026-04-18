#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <quadmath.h>

using namespace std;

// ============================================================
// __float128 转字符串
// ============================================================
string qtos(__float128 x, const char* fmt = "%.36Qg") {
    char buf[256];
    quadmath_snprintf(buf, sizeof(buf), fmt, x);
    return string(buf);
}

// ============================================================
// 直接 Taylor 部分和
// S(x,N) = sum_{n=0}^{N} x^n / n!
// ============================================================
template <typename T>
T S_direct(T x, int N) {
    T sum = 1;
    T term = 1;
    for (int n = 1; n <= N; ++n) {
        term = term * x / static_cast<T>(n);
        sum += term;
    }
    return sum;
}

// ============================================================
// 改进方法：x<0 时使用 e^x = 1 / S(-x,N)
// ============================================================
template <typename T>
T S_improved(T x, int N) {
    if (x < 0) {
        return static_cast<T>(1) / S_direct(-x, N);
    } else {
        return S_direct(x, N);
    }
}

// ============================================================
// 精确 exp
// ============================================================
float exact_exp(float x) { return expf(x); }
double exact_exp(double x) { return exp(x); }
__float128 exact_exp(__float128 x) { return expq(x); }

// ============================================================
// 误差
// ============================================================
template <typename T>
T abs_error(T approx, T exact) {
    return fabs(approx - exact);
}

__float128 abs_error(__float128 approx, __float128 exact) {
    return fabsq(approx - exact);
}

template <typename T>
T rel_error(T approx, T exact) {
    if (exact == static_cast<T>(0)) return static_cast<T>(0);
    return fabs((approx - exact) / exact);
}

__float128 rel_error(__float128 approx, __float128 exact) {
    if (exact == 0) return 0;
    return fabsq((approx - exact) / exact);
}

// ============================================================
// 普通打印
// ============================================================
template <typename T>
void print_value(const string& label, T value) {
    cout << left << setw(32) << label
         << setprecision(18) << value << "\n";
}

void print_value_q(const string& label, __float128 value) {
    cout << left << setw(32) << label
         << qtos(value, "%.36Qg") << "\n";
}

// ============================================================
// 每个测试点的结果打印（float）
// ============================================================
void print_case_float(float x, int N) {
    float direct = S_direct<float>(x, N);
    float improved = S_improved<float>(x, N);
    float exact = exact_exp(x);

    cout << "------------------------------------------------------------\n";
    cout << "x = " << x << "\n";
    print_value("Direct Taylor S(x,N) =", direct);
    print_value("Improved value       =", improved);
    print_value("std::exp(x)          =", exact);
    print_value("Abs error of direct  =", abs_error(direct, exact));
    print_value("Rel error of direct  =", rel_error(direct, exact));
    print_value("Abs error improved   =", abs_error(improved, exact));
    print_value("Rel error improved   =", rel_error(improved, exact));
    cout << "\n";
}

// ============================================================
// 每个测试点的结果打印（double）
// ============================================================
void print_case_double(double x, int N) {
    double direct = S_direct<double>(x, N);
    double improved = S_improved<double>(x, N);
    double exact = exact_exp(x);

    cout << "------------------------------------------------------------\n";
    cout << "x = " << x << "\n";
    print_value("Direct Taylor S(x,N) =", direct);
    print_value("Improved value       =", improved);
    print_value("std::exp(x)          =", exact);
    print_value("Abs error of direct  =", abs_error(direct, exact));
    print_value("Rel error of direct  =", rel_error(direct, exact));
    print_value("Abs error improved   =", abs_error(improved, exact));
    print_value("Rel error improved   =", rel_error(improved, exact));
    cout << "\n";
}

// ============================================================
// 每个测试点的结果打印（quadruple）
// ============================================================
void print_case_quad(__float128 x, int N) {
    __float128 direct = S_direct<__float128>(x, N);
    __float128 improved = S_improved<__float128>(x, N);
    __float128 exact = exact_exp(x);

    cout << "------------------------------------------------------------\n";
    cout << "x = " << qtos(x, "%.20Qg") << "\n";
    print_value_q("Direct Taylor S(x,N) =", direct);
    print_value_q("Improved value       =", improved);
    print_value_q("expq(x)              =", exact);
    print_value_q("Abs error of direct  =", abs_error(direct, exact));
    print_value_q("Rel error of direct  =", rel_error(direct, exact));
    print_value_q("Abs error improved   =", abs_error(improved, exact));
    print_value_q("Rel error improved   =", rel_error(improved, exact));
    cout << "\n";
}

// ============================================================
// 每组算完后的分析
// ============================================================
void explain_single_group() {
    cout << "================ SINGLE PRECISION ANALYSIS =================\n";
    cout << "1) In single precision, the number of significant digits is limited.\n";
    cout << "2) For x > 0, the Taylor series terms are all positive, so the summation is relatively stable.\n";
    cout << "3) For x < 0, the series becomes alternating, which means large positive and negative terms\n";
    cout << "   cancel each other.\n";
    cout << "4) This cancellation causes severe round-off loss in float arithmetic.\n";
    cout << "5) Therefore, direct S(x,N) is usually poor for negative x, especially for x = -10.\n";
    cout << "6) The improved formula e^x = 1 / S(-x,N) avoids this cancellation and is much more stable.\n";
    cout << "============================================================\n\n";
}

void explain_double_group() {
    cout << "================ DOUBLE PRECISION ANALYSIS =================\n";
    cout << "1) Double precision keeps more significant digits than single precision.\n";
    cout << "2) Therefore, round-off error is smaller and the results are generally better.\n";
    cout << "3) However, for x < 0, cancellation still exists in the direct Taylor sum.\n";
    cout << "4) So even though double is better than float, direct summation may still lose accuracy\n";
    cout << "   when x is sufficiently negative.\n";
    cout << "5) The improved reciprocal method remains more reliable for negative x.\n";
    cout << "============================================================\n\n";
}

void explain_quad_group() {
    cout << "============== QUADRUPLE PRECISION ANALYSIS ================\n";
    cout << "1) Quadruple precision preserves far more digits than float or double.\n";
    cout << "2) Because of that, the direct Taylor result is improved significantly.\n";
    cout << "3) Even for negative x, cancellation becomes less destructive than before.\n";
    cout << "4) However, the direct formula is still not the best numerical algorithm for x < 0.\n";
    cout << "5) The identity e^x = 1 / S(-x,N) is still the more stable and preferred approach.\n";
    cout << "6) This shows that increasing precision helps, but using a stable algorithm helps even more.\n";
    cout << "============================================================\n\n";
}

// ============================================================
// 整体总分析
// ============================================================
void final_explanation() {
    cout << "==================== FINAL DISCUSSION ======================\n";
    cout << "Part (a): Why does direct S(x,N) fail for x < 0?\n";
    cout << "Because when x is negative, the Taylor series becomes alternating:\n";
    cout << "1 - |x| + |x|^2/2! - |x|^3/3! + ...\n";
    cout << "The intermediate terms can be large, but the final answer e^x is small.\n";
    cout << "So the final result is obtained by subtracting nearly equal large numbers.\n";
    cout << "This is called catastrophic cancellation, and it destroys significant digits.\n\n";

    cout << "Part (b): Why does the improved method work?\n";
    cout << "For x < 0, write e^x = 1 / e^{-x}.\n";
    cout << "Now -x is positive, so S(-x,N) is a sum of positive terms.\n";
    cout << "That summation is numerically stable, and then taking the reciprocal gives a much better result.\n\n";

    cout << "What happens in single, double, and quadruple precision?\n";
    cout << "1) Single precision is affected the most by cancellation.\n";
    cout << "2) Double precision reduces round-off and improves accuracy.\n";
    cout << "3) Quadruple precision improves further and can preserve many more digits.\n";
    cout << "4) But the main lesson is that a numerically stable method is more important than only using higher precision.\n\n";

    cout << "Overall conclusion:\n";
    cout << "For x > 0, direct Taylor summation is fine.\n";
    cout << "For x < 0, the stable choice is e^x = 1 / S(-x,N).\n";
    cout << "============================================================\n";
}

// ============================================================
// 主程序
// ============================================================
int main() {
    cout << fixed;
    int N = 20;

    // ======================== SINGLE =========================
    cout << "######################## FLOAT TESTS ########################\n\n";
    vector<float> xs_float = {10.0f, 2.0f, -2.0f, -10.0f};
    for (float x : xs_float) {
        print_case_float(x, N);
    }
    explain_single_group();

    // ======================== DOUBLE =========================
    cout << "####################### DOUBLE TESTS ########################\n\n";
    vector<double> xs_double = {10.0, 2.0, -2.0, -10.0};
    for (double x : xs_double) {
        print_case_double(x, N);
    }
    explain_double_group();

    // ====================== QUADRUPLE ========================
    cout << "##################### QUADRUPLE TESTS #######################\n\n";
    vector<__float128> xs_quad = {10.0Q, 2.0Q, -2.0Q, -10.0Q};
    for (__float128 x : xs_quad) {
        print_case_quad(x, N);
    }
    explain_quad_group();

    // ======================== FINAL ==========================
    final_explanation();

    return 0;
}
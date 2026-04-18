#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>

using namespace std;

// ============================================================
// 理论值：a_n = ∫_0^1 x^n e^(x-1) dx
// 用 Simpson 积分
// ============================================================
long double exactValueIntegral(int n) {
    const int M = 200000;
    long double h = 1.0L / M;
    long double sum = 0.0L;

    for (int i = 0; i <= M; i++) {
        long double x = i * h;
        long double fx = powl(x, n) * expl(x - 1.0L);

        if (i == 0 || i == M)
            sum += fx;
        else if (i % 2 == 0)
            sum += 2.0L * fx;
        else
            sum += 4.0L * fx;
    }

    return sum * h / 3.0L;
}

// ============================================================
// 原程序真实递推（double）
// ============================================================
vector<double> computeRealRecursion(int N) {
    vector<double> a(N + 1);
    a[0] = (exp(1.0) - 1.0) / exp(1.0);

    for (int n = 0; n < N; ++n) {
        a[n + 1] = 1.0 - (n + 1) * a[n];
    }
    return a;
}

// ============================================================
// 非负大整数，base = 1e9
// ============================================================
struct BigUInt {
    static const uint32_t BASE = 1000000000U;
    vector<uint32_t> d; // little-endian

    void trim() {
        while (!d.empty() && d.back() == 0) d.pop_back();
    }

    bool isZero() const {
        return d.empty();
    }

    static BigUInt fromUint64(unsigned long long x) {
        BigUInt a;
        if (x == 0) return a;
        while (x) {
            a.d.push_back((uint32_t)(x % BASE));
            x /= BASE;
        }
        return a;
    }

    void mulSmall(uint32_t m) {
        if (isZero() || m == 0) {
            d.clear();
            return;
        }
        if (m == 1) return;

        unsigned long long carry = 0;
        for (size_t i = 0; i < d.size(); ++i) {
            unsigned long long cur = 1ULL * d[i] * m + carry;
            d[i] = (uint32_t)(cur % BASE);
            carry = cur / BASE;
        }
        while (carry) {
            d.push_back((uint32_t)(carry % BASE));
            carry /= BASE;
        }
    }

    uint32_t divSmall(uint32_t m) {
        unsigned long long rem = 0;
        for (int i = (int)d.size() - 1; i >= 0; --i) {
            unsigned long long cur = rem * BASE + d[i];
            d[i] = (uint32_t)(cur / m);
            rem = cur % m;
        }
        trim();
        return (uint32_t)rem;
    }

    string toDecimalString() const {
        if (isZero()) return "0";
        BigUInt tmp = *this;
        vector<uint32_t> parts;

        while (!tmp.isZero()) {
            parts.push_back(tmp.divSmall(BASE));
        }

        string s = to_string(parts.back());
        for (int i = (int)parts.size() - 2; i >= 0; --i) {
            string t = to_string(parts[i]);
            s += string(9 - t.size(), '0') + t;
        }
        return s;
    }

    // 仅用于显示近似值
    long double toLongDoubleApprox() const {
        if (isZero()) return 0.0L;

        string s = toDecimalString();
        int take = min((int)s.size(), 18);

        long double mant = 0.0L;
        for (int i = 0; i < take; ++i) {
            mant = mant * 10.0L + (s[i] - '0');
        }

        int exp10 = (int)s.size() - take;
        return mant * powl(10.0L, exp10);
    }
};

// ============================================================
// 有符号大整数
// ============================================================
struct SignedBigInt {
    int sign;   // -1, 0, +1
    BigUInt mag;

    SignedBigInt() : sign(0) {}

    bool isZero() const {
        return sign == 0 || mag.isZero();
    }

    static SignedBigInt fromScaledLongDouble(long double x, int SHIFT) {
        SignedBigInt a;
        if (x == 0.0L) return a;

        long double scaled = ldexpl(x, SHIFT);  // x * 2^SHIFT
        long double rounded = roundl(scaled);

        if (rounded == 0.0L) return a;

        a.sign = (rounded > 0.0L ? 1 : -1);
        unsigned long long v = (unsigned long long) llroundl(fabsl(rounded));
        a.mag = BigUInt::fromUint64(v);
        return a;
    }

    void mulSmallSigned(int k) {
        if (isZero() || k == 0) {
            sign = 0;
            mag.d.clear();
            return;
        }
        if (k < 0) {
            sign = -sign;
            k = -k;
        }
        mag.mulSmall((uint32_t)k);
        if (mag.isZero()) sign = 0;
    }

    // 仅用于显示近似值：先转 long double 再除 2^SHIFT
    long double toApproxError(int SHIFT) const {
        if (isZero()) return 0.0L;
        long double v = mag.toLongDoubleApprox();
        if (sign < 0) v = -v;
        return ldexpl(v, -SHIFT);
    }

    string toScientificString(int SHIFT, int sigDigits = 12) const {
        long double v = toApproxError(SHIFT);
        ostringstream oss;
        oss << scientific << setprecision(sigDigits) << v;
        return oss.str();
    }
};

int main() {
    int N;
    cout << "请输入 N: ";
    cin >> N;

    // 这个值可调：
    // M_n = round(e_n * 2^SHIFT)
    // SHIFT 越大，初值误差保留越细
    const int SHIFT = 80;

    // -----------------------------
    // 理论值
    // -----------------------------
    vector<long double> a_theory(N + 1);
    a_theory[0] = 1.0L - 1.0L / expl(1.0L);
    for (int n = 1; n <= N; ++n) {
        a_theory[n] = exactValueIntegral(n);
    }

    // -----------------------------
    // 原程序实际结果
    // -----------------------------
    vector<double> a_real = computeRealRecursion(N);

    vector<long double> real_err(N + 1);
    for (int n = 0; n <= N; ++n) {
        real_err[n] = (long double)a_real[n] - a_theory[n];
    }

    // -----------------------------
    // 误差模型：
    // e_n ≈ (-1)^n n! e0
    //
    // 这里 e0 取“原程序实际初值”与“理论初值”的差：
    // e0 = a0_real - a0_theory
    //
    // 内部不直接存 e_n，而存
    // M_n = round(e_n * 2^SHIFT)
    // 于是：
    // M_{n+1} = -(n+1) M_n
    // -----------------------------
    long double a0_theory = a_theory[0];
    double a0_real = a_real[0];
    long double e0 = (long double)a0_real - a0_theory;

    vector<SignedBigInt> pred(N + 1);
    pred[0] = SignedBigInt::fromScaledLongDouble(e0, SHIFT);

    for (int n = 0; n < N; ++n) {
        pred[n + 1] = pred[n];
        pred[n + 1].mulSmallSigned(-(n + 1));
    }

    // -----------------------------
    // Part 1
    // -----------------------------
    cout << "\n==================== Part 1 ====================\n";
    cout << "原程序递推值、理论值、真实误差\n\n";

    cout << scientific << setprecision(10);
    cout << setw(5)  << "n"
         << setw(20) << "a_n(real)"
         << setw(24) << "a_n(theory)"
         << setw(24) << "real_error"
         << "\n";

    for (int n = 0; n <= N; ++n) {
        cout << setw(5)  << n
             << setw(20) << a_real[n]
             << setw(24) << a_theory[n]
             << setw(24) << real_err[n]
             << "\n";
    }

    // -----------------------------
    // Part 2
    // -----------------------------
    cout << "\n==================== Part 2 ====================\n";
    cout << "误差模型：e_n ≈ (-1)^n n! e0\n";
    cout << "其中 e0 = a0_real - a0_theory\n";
    cout << "为了避免浮点爆炸，内部采用：\n";
    cout << "M_n = round(e_n * 2^SHIFT), SHIFT = " << SHIFT << "\n";
    cout << "然后精确递推：M_(n+1) = -(n+1) M_n\n";
    cout << "输出时再除以 2^SHIFT。\n\n";

    cout << scientific << setprecision(18);
    cout << "a0_theory = " << a0_theory << "\n";
    cout << "a0_real   = " << a0_real << "\n";
    cout << "e0        = " << e0 << "\n";
    cout << "pred[0]   = " << pred[0].toScientificString(SHIFT, 12) << "\n\n";

    cout << scientific << setprecision(10);
    cout << setw(5)  << "n"
         << setw(24) << "predicted_error"
         << "\n";

    for (int n = 0; n <= N; ++n) {
        cout << setw(5)  << n
             << setw(24) << pred[n].toApproxError(SHIFT)
             << "\n";
    }

    // -----------------------------
    // Part 3
    // -----------------------------
    cout << "\n==================== Part 3 ====================\n";
    cout << "将预测误差与真实误差比较\n\n";

    cout << setw(5)  << "n"
         << setw(24) << "real_error"
         << setw(24) << "predicted_error"
         << setw(24) << "difference"
         << "\n";

    for (int n = 0; n <= N; ++n) {
        long double est = pred[n].toApproxError(SHIFT);
        long double diff = real_err[n] - est;

        cout << setw(5)  << n
             << setw(24) << real_err[n]
             << setw(24) << est
             << setw(24) << diff
             << "\n";
    }

    return 0;
}
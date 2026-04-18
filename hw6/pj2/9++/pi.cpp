// pi_chudnovsky_gmp.cpp
#include <gmpxx.h>
#include <omp.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <string>

static const mpz_class A("13591409");
static const mpz_class B("545140134");
static const mpz_class C3_OVER_24("10939058860032000"); // 640320^3 / 24

struct BSResult {
    mpz_class P;
    mpz_class Q;
    mpz_class T;
};

static mpz_class pow_u64(mpz_class base, uint64_t exp) {
    mpz_class result = 1;
    while (exp > 0) {
        if (exp & 1) result *= base;
        exp >>= 1;
        if (exp) base *= base;
    }
    return result;
}

static void bs_serial(uint64_t a, uint64_t b, BSResult& r) {
    r.P = 1;
    r.Q = 1;
    r.T = 0;

    for (uint64_t k = a; k < b; ++k) {
        mpz_class pk, qk, tk;

        if (k == 0) {
            pk = 1;
            qk = 1;
            tk = A;
        } else {
            pk = mpz_class(6 * k - 5);
            pk *= mpz_class(2 * k - 1);
            pk *= mpz_class(6 * k - 1);

            qk = mpz_class(k);
            qk *= mpz_class(k);
            qk *= mpz_class(k);
            qk *= C3_OVER_24;

            tk = pk * (A + B * mpz_class(k));
            if (k & 1) tk = -tk;
        }

        mpz_class newT = r.T * qk + r.P * tk;
        r.P *= pk;
        r.Q *= qk;
        r.T.swap(newT);
    }
}

static void bs_range(uint64_t a, uint64_t b, BSResult& r,
                     uint64_t leaf_size, uint64_t task_cutoff) {
    uint64_t len = b - a;
    if (len <= leaf_size) {
        bs_serial(a, b, r);
        return;
    }

    uint64_t m = a + len / 2;
    BSResult left, right;

    if (len >= task_cutoff) {
        #pragma omp task shared(left)
        bs_range(a, m, left, leaf_size, task_cutoff);

        #pragma omp task shared(right)
        bs_range(m, b, right, leaf_size, task_cutoff);

        #pragma omp taskwait
    } else {
        bs_range(a, m, left, leaf_size, task_cutoff);
        bs_range(m, b, right, leaf_size, task_cutoff);
    }

    r.P = left.P * right.P;
    r.Q = left.Q * right.Q;
    r.T = left.T * right.Q + left.P * right.T;
}

static uint64_t estimate_terms(uint64_t digits) {
    // Chudnovsky 每项大约增加 14.1816474627 位
    return static_cast<uint64_t>(digits / 14.181647462725477) + 10;
}

static void write_pi_decimal(const mpz_class& pi_scaled, uint64_t digits,
                             const std::string& out_path) {
    mpz_class scale = pow_u64(10, digits);
    mpz_class int_part = pi_scaled / scale;
    mpz_class frac_part = pi_scaled % scale;

    FILE* fp = std::fopen(out_path.c_str(), "wb");
    if (!fp) {
        std::perror("fopen");
        std::exit(1);
    }

    mpz_out_str(fp, 10, int_part.get_mpz_t());
    std::fputc('.', fp);

    uint64_t frac_len = 1;
    if (frac_part != 0) {
        frac_len = mpz_sizeinbase(frac_part.get_mpz_t(), 10);
    }

    for (uint64_t i = frac_len; i < digits; ++i) {
        std::fputc('0', fp);
    }
    mpz_out_str(fp, 10, frac_part.get_mpz_t());
    std::fputc('\n', fp);

    std::fclose(fp);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <digits> <output.txt> [leaf_size] [task_cutoff]\n";
        return 1;
    }

    const uint64_t digits = std::strtoull(argv[1], nullptr, 10);
    const std::string out_path = argv[2];
    const uint64_t leaf_size = (argc >= 4) ? std::strtoull(argv[3], nullptr, 10) : 64;
    const uint64_t task_cutoff = (argc >= 5) ? std::strtoull(argv[4], nullptr, 10) : (1ull << 15);

    // 为了稳妥，留一些 guard digits
    const uint64_t guard = 32;
    const uint64_t work_digits = digits + guard;
    const uint64_t n_terms = estimate_terms(work_digits);

    std::cerr << "digits      = " << digits << "\n";
    std::cerr << "work_digits = " << work_digits << "\n";
    std::cerr << "n_terms     = " << n_terms << "\n";
    std::cerr << "threads     = " << omp_get_max_threads() << "\n";

    auto t0 = std::chrono::steady_clock::now();

    BSResult r;
    #pragma omp parallel
    {
        #pragma omp single nowait
        bs_range(0, n_terms, r, leaf_size, task_cutoff);
    }

    auto t1 = std::chrono::steady_clock::now();
    std::cerr << "[1/3] binary splitting done.\n";

    // sqrt(10005) * 10^work_digits
    mpz_class scale2 = pow_u64(10, 2 * work_digits);
    scale2 *= 10005;

    mpz_class sqrt_c;
    mpz_sqrt(sqrt_c.get_mpz_t(), scale2.get_mpz_t());

    auto t2 = std::chrono::steady_clock::now();
    std::cerr << "[2/3] sqrt done.\n";

    mpz_class pi_scaled = (r.Q * 426880 * sqrt_c) / r.T;

    // 去掉 guard digits
    mpz_class div_guard = pow_u64(10, guard);
    pi_scaled /= div_guard;

    auto t3 = std::chrono::steady_clock::now();
    std::cerr << "[3/3] division done.\n";

    write_pi_decimal(pi_scaled, digits, out_path);

    auto t4 = std::chrono::steady_clock::now();

    auto ms_bs   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    auto ms_sqrt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    auto ms_div  = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
    auto ms_out  = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();

    std::cerr << "binary split : " << ms_bs   / 1000.0 << " s\n";
    std::cerr << "sqrt         : " << ms_sqrt / 1000.0 << " s\n";
    std::cerr << "final divide : " << ms_div  / 1000.0 << " s\n";
    std::cerr << "write output : " << ms_out  / 1000.0 << " s\n";

    return 0;
}
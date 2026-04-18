#include <gmpxx.h>
#include <omp.h>
#include <sys/resource.h>
#include <sys/time.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

static const mpz_class A("13591409");
static const mpz_class B("545140134");
static const mpz_class C3_OVER_24("10939058860032000"); // 640320^3 / 24

struct BSResult {
    mpz_class P;
    mpz_class Q;
    mpz_class T;
};

struct Timer {
    using clock = std::chrono::steady_clock;
    clock::time_point start_tp;
    Timer() : start_tp(clock::now()) {}
    void reset() { start_tp = clock::now(); }
    double elapsed_sec() const {
        return std::chrono::duration<double>(clock::now() - start_tp).count();
    }
};

enum class OutputMode {
    NONE,
    TXT,
    RAW
};

static void die(const std::string& msg) {
    std::cerr << "Error: " << msg << "\n";
    std::exit(1);
}

static OutputMode parse_output_mode(const std::string& s) {
    if (s == "none") return OutputMode::NONE;
    if (s == "txt")  return OutputMode::TXT;
    if (s == "raw")  return OutputMode::RAW;
    die("invalid output mode: " + s + " (expected none|txt|raw)");
    return OutputMode::NONE;
}

static const char* output_mode_name(OutputMode mode) {
    switch (mode) {
        case OutputMode::NONE: return "none";
        case OutputMode::TXT:  return "txt";
        case OutputMode::RAW:  return "raw";
    }
    return "unknown";
}

static mpz_class pow10_ui(unsigned long exp) {
    mpz_class x;
    mpz_ui_pow_ui(x.get_mpz_t(), 10UL, exp);
    return x;
}

static uint64_t estimate_terms(uint64_t digits) {
    return static_cast<uint64_t>(digits / 14.181647462725477) + 10;
}

static size_t file_size_bytes(const std::string& path) {
    FILE* fp = std::fopen(path.c_str(), "rb");
    if (!fp) return 0;
    if (std::fseek(fp, 0, SEEK_END) != 0) {
        std::fclose(fp);
        return 0;
    }
    long long pos = std::ftell(fp);
    std::fclose(fp);
    return pos > 0 ? static_cast<size_t>(pos) : 0;
}

static double current_rss_mib() {
    FILE* fp = std::fopen("/proc/self/status", "r");
    if (!fp) return -1.0;

    char line[512];
    long long kb = -1;
    while (std::fgets(line, sizeof(line), fp)) {
        if (std::strncmp(line, "VmRSS:", 6) == 0) {
            std::sscanf(line + 6, "%lld", &kb);
            break;
        }
    }
    std::fclose(fp);
    if (kb < 0) return -1.0;
    return static_cast<double>(kb) / 1024.0;
}

static double peak_rss_mib() {
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0) return -1.0;
#if defined(__linux__)
    return static_cast<double>(ru.ru_maxrss) / 1024.0; // Linux returns KB
#else
    return static_cast<double>(ru.ru_maxrss) / (1024.0 * 1024.0);
#endif
}

static void print_stage(const std::string& name, double sec) {
    std::cerr << std::fixed << std::setprecision(3)
              << "[TIME] " << std::setw(18) << std::left << name
              << " : " << sec << " s\n";
}

static void print_mem(const std::string& name) {
    const double rss = current_rss_mib();
    const double peak = peak_rss_mib();
    std::cerr << std::fixed << std::setprecision(1)
              << "[MEM ] " << std::setw(18) << std::left << name
              << " : rss=";
    if (rss >= 0.0) std::cerr << rss << " MiB";
    else std::cerr << "n/a";
    std::cerr << ", peak=";
    if (peak >= 0.0) std::cerr << peak << " MiB";
    else std::cerr << "n/a";
    std::cerr << "\n";
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
            if (k & 1ULL) tk = -tk;
        }

        mpz_class newT = r.T * qk + r.P * tk;
        r.P *= pk;
        r.Q *= qk;
        r.T.swap(newT);
    }
}

static void bs_range(uint64_t a, uint64_t b, BSResult& r,
                     uint64_t leaf_size, uint64_t task_cutoff) {
    const uint64_t len = b - a;
    if (len <= leaf_size) {
        bs_serial(a, b, r);
        return;
    }

    const uint64_t m = a + len / 2;
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

static void write_pi_txt_fast(const mpz_class& pi_scaled, uint64_t digits,
                              const std::string& out_path) {
    // Fast path for moderate outputs: convert once, then stream with large stdio buffer.
    // Avoids extra huge division/modulo used only for formatting.
    std::string s = pi_scaled.get_str(10);
    if (s.empty()) die("internal error: empty decimal string");

    FILE* fp = std::fopen(out_path.c_str(), "wb");
    if (!fp) {
        std::perror("fopen");
        std::exit(1);
    }

    static const size_t BUF_SIZE = 1u << 20;
    std::unique_ptr<char[]> buf(new char[BUF_SIZE]);
    std::setvbuf(fp, buf.get(), _IOFBF, BUF_SIZE);

    std::fputc(s[0], fp);
    std::fputc('.', fp);

    const size_t frac_len = (s.size() >= 2) ? (s.size() - 1) : 0;
    if (digits > frac_len) {
        size_t pad = static_cast<size_t>(digits - frac_len);
        std::string zeros(std::min<size_t>(pad, 1u << 20), '0');
        while (pad > 0) {
            const size_t chunk = std::min(pad, zeros.size());
            std::fwrite(zeros.data(), 1, chunk, fp);
            pad -= chunk;
        }
    }

    if (frac_len > 0) {
        std::fwrite(s.data() + 1, 1, frac_len, fp);
    }
    std::fputc('\n', fp);
    std::fclose(fp);
}

static void write_pi_raw(const mpz_class& pi_scaled, const std::string& out_path) {
    FILE* fp = std::fopen(out_path.c_str(), "wb");
    if (!fp) {
        std::perror("fopen");
        std::exit(1);
    }
    if (mpz_out_raw(fp, pi_scaled.get_mpz_t()) == 0) {
        std::fclose(fp);
        die("mpz_out_raw failed");
    }
    std::fclose(fp);
}

static void print_usage(const char* argv0) {
    std::cerr
        << "Usage:\n"
        << "  " << argv0 << " <digits> [output_path] [leaf_size] [task_cutoff] [output_mode]\n\n"
        << "Arguments:\n"
        << "  digits       target decimal digits after decimal point\n"
        << "  output_path  output file path (omit or use '-' when output_mode=none)\n"
        << "  leaf_size    serial leaf block size, default 256\n"
        << "  task_cutoff  recurse in parallel only when len >= cutoff, default 262144\n"
        << "  output_mode  none | txt | raw, default txt\n\n"
        << "Examples:\n"
        << "  " << argv0 << " 100000000 pi_100m.txt\n"
        << "  " << argv0 << " 1000000000 - 256 262144 none\n"
        << "  " << argv0 << " 1000000000 pi_1b.raw 256 262144 raw\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    const uint64_t digits = std::strtoull(argv[1], nullptr, 10);
    if (digits == 0) die("digits must be > 0");

    const std::string out_path = (argc >= 3) ? argv[2] : "pi.txt";
    const uint64_t leaf_size = (argc >= 4) ? std::strtoull(argv[3], nullptr, 10) : 256ULL;
    const uint64_t task_cutoff = (argc >= 5) ? std::strtoull(argv[4], nullptr, 10) : (1ULL << 18);
    const OutputMode out_mode = (argc >= 6) ? parse_output_mode(argv[5]) : OutputMode::TXT;

    if (leaf_size == 0) die("leaf_size must be > 0");
    if (task_cutoff == 0) die("task_cutoff must be > 0");
    if (out_mode != OutputMode::NONE && out_path == "-") die("output_path cannot be '-' unless output_mode=none");

    const uint64_t guard_digits = 32;
    const uint64_t work_digits = digits + guard_digits;
    const uint64_t n_terms = estimate_terms(work_digits);

    std::cerr << "digits       = " << digits << "\n";
    std::cerr << "work_digits  = " << work_digits << "\n";
    std::cerr << "guard_digits = " << guard_digits << "\n";
    std::cerr << "n_terms      = " << n_terms << "\n";
    std::cerr << "threads      = " << omp_get_max_threads() << "\n";
    std::cerr << "leaf_size    = " << leaf_size << "\n";
    std::cerr << "task_cutoff  = " << task_cutoff << "\n";
    std::cerr << "output_mode  = " << output_mode_name(out_mode) << "\n";
    if (out_mode != OutputMode::NONE) {
        std::cerr << "output_path  = " << out_path << "\n";
    }
    std::cerr << "\n";

    Timer wall_timer;
    Timer compute_timer;

    // Stage 1: binary splitting
    Timer timer_bs;
    BSResult r;
    #pragma omp parallel
    {
        #pragma omp single nowait
        bs_range(0, n_terms, r, leaf_size, task_cutoff);
    }
    const double sec_bs = timer_bs.elapsed_sec();
    print_stage("binary splitting", sec_bs);
    print_mem("after binary split");

    // Stage 2: sqrt(10005) * 10^work_digits
    Timer timer_sqrt;
    mpz_class scale2 = pow10_ui(static_cast<unsigned long>(2 * work_digits));
    scale2 *= 10005;

    mpz_class sqrt_c;
    mpz_sqrt(sqrt_c.get_mpz_t(), scale2.get_mpz_t());
    const double sec_sqrt = timer_sqrt.elapsed_sec();
    print_stage("sqrt", sec_sqrt);
    print_mem("after sqrt");

    // Stage 3: final division
    Timer timer_div;
    mpz_class pi_scaled = (r.Q * 426880 * sqrt_c) / r.T;
    if (guard_digits > 0) {
        mpz_class div_guard = pow10_ui(static_cast<unsigned long>(guard_digits));
        pi_scaled /= div_guard;
    }
    const double sec_div = timer_div.elapsed_sec();
    print_stage("final division", sec_div);
    print_mem("after division");

    const double sec_compute = compute_timer.elapsed_sec();
    print_stage("COMPUTE TOTAL", sec_compute);

    // Optional output stage, intentionally excluded from compute total
    double sec_out = 0.0;
    size_t out_bytes = 0;
    if (out_mode != OutputMode::NONE) {
        Timer timer_out;
        if (out_mode == OutputMode::TXT) {
            write_pi_txt_fast(pi_scaled, digits, out_path);
        } else {
            write_pi_raw(pi_scaled, out_path);
        }
        sec_out = timer_out.elapsed_sec();
        out_bytes = file_size_bytes(out_path);
        print_stage("write output", sec_out);
        print_mem("after write");
    }

    const double sec_wall = wall_timer.elapsed_sec();
    print_stage("WALL TOTAL", sec_wall);

    std::cerr << "\n";
    std::cerr << std::fixed << std::setprecision(3);
    const double terms_per_sec = (sec_bs > 0.0) ? (double)n_terms / sec_bs : 0.0;
    const double digits_per_sec_compute = (sec_compute > 0.0) ? (double)digits / sec_compute : 0.0;
    const double digits_per_sec_wall = (sec_wall > 0.0) ? (double)digits / sec_wall : 0.0;

    std::cerr << "[STAT] terms/sec             : " << terms_per_sec << "\n";
    std::cerr << "[STAT] digits/sec(compute)   : " << digits_per_sec_compute << "\n";
    std::cerr << "[STAT] digits/sec(wall)      : " << digits_per_sec_wall << "\n";
    if (out_mode != OutputMode::NONE) {
        const double mib = static_cast<double>(out_bytes) / (1024.0 * 1024.0);
        const double write_mib_per_sec = (sec_out > 0.0) ? mib / sec_out : 0.0;
        std::cerr << "[STAT] output size(MiB)      : " << mib << "\n";
        std::cerr << "[STAT] write speed(MiB/s)    : " << write_mib_per_sec << "\n";
    }
    std::cerr << "[STAT] peak RSS(MiB)         : " << peak_rss_mib() << "\n";

    return 0;
}

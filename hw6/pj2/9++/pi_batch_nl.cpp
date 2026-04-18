#include <gmpxx.h>
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/resource.h>
#include <unistd.h>
#include <vector>

static const mpz_class A("13591409");
static const mpz_class B("545140134");
static const mpz_class C3_OVER_24("10939058860032000"); // 640320^3 / 24

struct Node {
    mpz_class P;
    mpz_class Q;
    mpz_class T;
};

struct Timer {
    using clock = std::chrono::steady_clock;
    clock::time_point t0;
    Timer() : t0(clock::now()) {}
    void reset() { t0 = clock::now(); }
    double sec() const {
        return std::chrono::duration<double>(clock::now() - t0).count();
    }
};

enum class OutputMode {
    TXT,
    RAW,
    NONE,
};

static void print_stage(const std::string& name, double sec) {
    std::cerr << std::fixed << std::setprecision(3)
              << "[TIME] " << std::setw(18) << std::left << name
              << " : " << sec << " s\n";
}

static void die_usage(const char* argv0) {
    std::cerr << "Usage: " << argv0
              << " <digits> <output_path|-> [leaf_block] [node_parallel_threshold]"
              << " [big_limb_threshold] [mode=txt|raw|none] [inner_merge_threads]\n"
              << "\n"
              << "  digits                  digits after decimal point\n"
              << "  output_path             output path, or '-' when mode=none\n"
              << "  leaf_block              OpenMP schedule chunk for leaf generation (default 2048)\n"
              << "  node_parallel_threshold use N-dim parallel when pair_count >= threshold * threads (default 2)\n"
              << "  big_limb_threshold      use L-dim merge when operand limbs are large enough (default 2048)\n"
              << "  mode                    txt | raw | none (default txt)\n"
              << "  inner_merge_threads     threads used inside one big merge in L-dim mode (default 4)\n";
    std::exit(1);
}

static OutputMode parse_mode(const std::string& s) {
    if (s == "txt") return OutputMode::TXT;
    if (s == "raw") return OutputMode::RAW;
    if (s == "none") return OutputMode::NONE;
    std::cerr << "Unknown mode: " << s << "\n";
    std::exit(1);
}

static mpz_class pow10_ui(uint64_t exp) {
    mpz_class x;
    mpz_ui_pow_ui(x.get_mpz_t(), 10ul, static_cast<unsigned long>(exp));
    return x;
}

static uint64_t estimate_terms(uint64_t digits) {
    return static_cast<uint64_t>(digits / 14.181647462725477) + 10;
}

static size_t limb_estimate(const Node& n) {
    size_t p = mpz_size(n.P.get_mpz_t());
    size_t q = mpz_size(n.Q.get_mpz_t());
    size_t t = mpz_size(n.T.get_mpz_t());
    return std::max({p, q, t});
}

static double current_rss_mib() {
    long rss_pages = 0;
    FILE* fp = std::fopen("/proc/self/statm", "r");
    if (!fp) return 0.0;
    long dummy = 0;
    if (std::fscanf(fp, "%ld %ld", &dummy, &rss_pages) != 2) {
        std::fclose(fp);
        return 0.0;
    }
    std::fclose(fp);
    const long page_size = sysconf(_SC_PAGESIZE);
    return static_cast<double>(rss_pages) * static_cast<double>(page_size) / (1024.0 * 1024.0);
}

static double peak_rss_mib() {
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0) return 0.0;
#if defined(__APPLE__)
    return static_cast<double>(ru.ru_maxrss) / (1024.0 * 1024.0);
#else
    return static_cast<double>(ru.ru_maxrss) / 1024.0;
#endif
}

static Node make_leaf(uint64_t k) {
    Node n;
    if (k == 0) {
        n.P = 1;
        n.Q = 1;
        n.T = A;
        return n;
    }

    n.P = mpz_class(6 * k - 5);
    n.P *= mpz_class(2 * k - 1);
    n.P *= mpz_class(6 * k - 1);

    n.Q = mpz_class(k);
    n.Q *= mpz_class(k);
    n.Q *= mpz_class(k);
    n.Q *= C3_OVER_24;

    n.T = n.P * (A + B * mpz_class(k));
    if (k & 1) n.T = -n.T;

    return n;
}

static void build_leaves(std::vector<Node>& leaves, uint64_t n_terms, int leaf_block) {
    leaves.resize(n_terms);
#pragma omp parallel for schedule(static, leaf_block)
    for (long long i = 0; i < static_cast<long long>(n_terms); ++i) {
        leaves[static_cast<size_t>(i)] = make_leaf(static_cast<uint64_t>(i));
    }
}

static void merge_serial(const Node& a, const Node& b, Node& out) {
    out.P = a.P * b.P;
    out.Q = a.Q * b.Q;
    out.T = a.T * b.Q + a.P * b.T;
}

static void merge_internal_parallel(const Node& a, const Node& b, Node& out, int inner_threads) {
    mpz_class p, q, t1, t2;
#pragma omp parallel sections num_threads(inner_threads)
    {
#pragma omp section
        { p = a.P * b.P; }
#pragma omp section
        { q = a.Q * b.Q; }
#pragma omp section
        { t1 = a.T * b.Q; }
#pragma omp section
        { t2 = a.P * b.T; }
    }
    out.P.swap(p);
    out.Q.swap(q);
    out.T = t1 + t2;
}

struct ReduceStats {
    double sec = 0.0;
    int levels = 0;
};

static ReduceStats reduce_tree_dynamic(std::vector<Node>& current,
                                       int max_threads,
                                       int node_parallel_threshold,
                                       size_t big_limb_threshold,
                                       int inner_merge_threads,
                                       bool verbose_levels) {
    ReduceStats stats;
    Timer timer;
    int level = 0;

    while (current.size() > 1) {
        const size_t node_count = current.size();
        const size_t pair_count = node_count / 2;
        const size_t next_count = (node_count + 1) / 2;
        const size_t approx_limbs = limb_estimate(current[0]);

        std::vector<Node> next(next_count);

        enum class Strategy { N_DIM, L_DIM, SERIAL };
        Strategy strategy = Strategy::SERIAL;

        if (pair_count >= static_cast<size_t>(node_parallel_threshold) * static_cast<size_t>(max_threads)) {
            strategy = Strategy::N_DIM;
        } else if (approx_limbs >= big_limb_threshold && max_threads > 1) {
            strategy = Strategy::L_DIM;
        }

        if (strategy == Strategy::N_DIM) {
#pragma omp parallel for schedule(static)
            for (long long i = 0; i < static_cast<long long>(pair_count); ++i) {
                const size_t idx = static_cast<size_t>(i);
                merge_serial(current[2 * idx], current[2 * idx + 1], next[idx]);
            }
        } else if (strategy == Strategy::L_DIM) {
            const int inner = std::max(2, std::min(inner_merge_threads, max_threads));
            for (size_t i = 0; i < pair_count; ++i) {
                merge_internal_parallel(current[2 * i], current[2 * i + 1], next[i], inner);
            }
        } else {
            for (size_t i = 0; i < pair_count; ++i) {
                merge_serial(current[2 * i], current[2 * i + 1], next[i]);
            }
        }

        if (node_count & 1) {
            next.back() = std::move(current.back());
        }

        if (verbose_levels) {
            const char* s = (strategy == Strategy::N_DIM)
                                ? "N-dim"
                                : (strategy == Strategy::L_DIM ? "L-dim" : "serial");
            std::cerr << "[LEVEL " << level << "] "
                      << "nodes=" << node_count
                      << " pairs=" << pair_count
                      << " limbs~=" << approx_limbs
                      << " strategy=" << s
                      << " rss=" << std::fixed << std::setprecision(1) << current_rss_mib() << " MiB\n";
        }

        current.swap(next);
        ++level;
    }

    stats.sec = timer.sec();
    stats.levels = level;
    return stats;
}

static void write_pi_decimal_fast(const mpz_class& pi_scaled,
                                  uint64_t digits,
                                  const std::string& out_path) {
    std::string s = pi_scaled.get_str(10);

    FILE* fp = std::fopen(out_path.c_str(), "wb");
    if (!fp) {
        std::perror("fopen");
        std::exit(1);
    }

    static const size_t BUF_SIZE = 1 << 20;
    std::setvbuf(fp, nullptr, _IOFBF, BUF_SIZE);

    if (s.empty()) {
        std::fclose(fp);
        std::cerr << "Unexpected empty pi string\n";
        std::exit(1);
    }

    std::fputc(s[0], fp);
    std::fputc('.', fp);

    size_t frac_len = (s.size() >= 2) ? (s.size() - 1) : 0;
    if (digits > frac_len) {
        size_t pad = static_cast<size_t>(digits - frac_len);
        std::string zeros(std::min<size_t>(pad, 1 << 20), '0');
        while (pad > 0) {
            size_t chunk = std::min(pad, zeros.size());
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
        std::cerr << "mpz_out_raw failed\n";
        std::exit(1);
    }
    std::fclose(fp);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        die_usage(argv[0]);
    }

    const uint64_t digits = std::strtoull(argv[1], nullptr, 10);
    const std::string out_path = argv[2];
    const int leaf_block = (argc >= 4) ? std::atoi(argv[3]) : 2048;
    const int node_parallel_threshold = (argc >= 5) ? std::atoi(argv[4]) : 2;
    const size_t big_limb_threshold = (argc >= 6) ? static_cast<size_t>(std::strtoull(argv[5], nullptr, 10)) : 2048;
    const OutputMode mode = (argc >= 7) ? parse_mode(argv[6]) : OutputMode::TXT;
    const int inner_merge_threads = (argc >= 8) ? std::atoi(argv[7]) : 4;

    if (leaf_block <= 0 || node_parallel_threshold <= 0 || inner_merge_threads <= 0) {
        std::cerr << "leaf_block, node_parallel_threshold, and inner_merge_threads must be positive.\n";
        return 1;
    }

    omp_set_dynamic(0);
    omp_set_max_active_levels(2);

    const int max_threads = omp_get_max_threads();
    const uint64_t guard = 32;
    const uint64_t work_digits = digits + guard;
    const uint64_t n_terms = estimate_terms(work_digits);

    std::cerr << "digits                  = " << digits << "\n";
    std::cerr << "work_digits             = " << work_digits << "\n";
    std::cerr << "n_terms                 = " << n_terms << "\n";
    std::cerr << "threads                 = " << max_threads << "\n";
    std::cerr << "leaf_block              = " << leaf_block << "\n";
    std::cerr << "node_parallel_threshold = " << node_parallel_threshold << "\n";
    std::cerr << "big_limb_threshold      = " << big_limb_threshold << "\n";
    std::cerr << "inner_merge_threads     = " << inner_merge_threads << "\n";
    std::cerr << "mode                    = "
              << (mode == OutputMode::TXT ? "txt" : mode == OutputMode::RAW ? "raw" : "none") << "\n\n";

    Timer wall_timer;
    Timer compute_timer;

    Timer t_leaf;
    std::vector<Node> current;
    build_leaves(current, n_terms, leaf_block);
    const double sec_leaf = t_leaf.sec();
    print_stage("leaf batch", sec_leaf);

    ReduceStats red = reduce_tree_dynamic(current,
                                          max_threads,
                                          node_parallel_threshold,
                                          big_limb_threshold,
                                          inner_merge_threads,
                                          true);
    print_stage("tree reduce", red.sec);

    Node& root = current[0];

    Timer t_sqrt;
    mpz_class scale2 = pow10_ui(2 * work_digits);
    scale2 *= 10005;

    mpz_class sqrt_c;
    mpz_sqrt(sqrt_c.get_mpz_t(), scale2.get_mpz_t());
    const double sec_sqrt = t_sqrt.sec();
    print_stage("sqrt", sec_sqrt);

    Timer t_div;
    mpz_class pi_scaled = (root.Q * 426880 * sqrt_c) / root.T;
    mpz_class div_guard = pow10_ui(guard);
    pi_scaled /= div_guard;
    const double sec_div = t_div.sec();
    print_stage("final division", sec_div);

    const double sec_compute = compute_timer.sec();
    print_stage("COMPUTE TOTAL", sec_compute);

    double sec_out = 0.0;
    if (mode != OutputMode::NONE) {
        Timer t_out;
        if (mode == OutputMode::TXT) {
            write_pi_decimal_fast(pi_scaled, digits, out_path);
        } else {
            write_pi_raw(pi_scaled, out_path);
        }
        sec_out = t_out.sec();
    }
    print_stage("write output", sec_out);

    const double sec_wall = wall_timer.sec();
    print_stage("WALL TOTAL", sec_wall);

    std::cerr << std::fixed << std::setprecision(3);
    std::cerr << "[STAT] reduce_levels       : " << red.levels << "\n";
    std::cerr << "[STAT] terms/sec           : " << (sec_compute > 0.0 ? static_cast<double>(n_terms) / sec_compute : 0.0) << "\n";
    std::cerr << "[STAT] digits/sec(compute) : " << (sec_compute > 0.0 ? static_cast<double>(digits) / sec_compute : 0.0) << "\n";
    std::cerr << "[STAT] rss_now(MiB)        : " << current_rss_mib() << "\n";
    std::cerr << "[STAT] rss_peak(MiB)       : " << peak_rss_mib() << "\n";

    return 0;
}

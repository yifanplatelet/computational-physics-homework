#define _POSIX_C_SOURCE 200809L
#define _FILE_OFFSET_BITS 64

#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/statvfs.h>
#include <time.h>

#include "cpu_backend.h"
#include "cuda_backend.h"

static double elapsed_seconds(struct timespec a, struct timespec b) {
    return (double)(b.tv_sec - a.tv_sec) +
           (double)(b.tv_nsec - a.tv_nsec) / 1e9;
}

/* 检查 width*height*3 是否溢出，并返回像素字节数 */
static int compute_image_bytes(int w, int h, size_t *out_bytes) {
    if (w <= 0 || h <= 0) return 0;

    uint64_t pixels = (uint64_t)(uint32_t)w * (uint64_t)(uint32_t)h;
    uint64_t bytes64 = pixels * 3ull;

    if (bytes64 > (uint64_t)SIZE_MAX) {
        return 0;
    }

    *out_bytes = (size_t)bytes64;
    return 1;
}

/* 打印文件系统剩余空间，方便排查磁盘不足 */
static void print_disk_free_space(const char *path) {
    struct statvfs vfs;
    if (statvfs(path, &vfs) == 0) {
        unsigned long long free_bytes =
            (unsigned long long)vfs.f_bavail * (unsigned long long)vfs.f_frsize;
        fprintf(stderr, "Disk free space (approx): %.2f GB\n",
                (double)free_bytes / (1024.0 * 1024.0 * 1024.0));
    } else {
        perror("statvfs");
    }
}

/* 写完后检查文件实际大小 */
static int verify_file_size(const char *filename, size_t header_bytes, size_t image_bytes) {
    struct stat st;
    if (stat(filename, &st) != 0) {
        perror("stat");
        return 0;
    }

    uint64_t expected = (uint64_t)header_bytes + (uint64_t)image_bytes;
    uint64_t actual = (uint64_t)st.st_size;

    fprintf(stderr, "Expected file size: %" PRIu64 " bytes\n", expected);
    fprintf(stderr, "Actual file size  : %" PRIu64 " bytes\n", actual);

    if (actual != expected) {
        fprintf(stderr, "ERROR: output file size mismatch.\n");
        return 0;
    }
    return 1;
}

/* 分块写，避免单次 fwrite 太大；并且严格检查写入结果 */
static int write_ppm_checked(const char *filename, const unsigned char *img, int w, int h) {
    size_t image_bytes = 0;
    if (!compute_image_bytes(w, h, &image_bytes)) {
        fprintf(stderr, "ERROR: image size overflow: %d x %d\n", w, h);
        return 0;
    }

    fprintf(stderr, "About to write PPM: %dx%d, pixel bytes = %" PRIu64 "\n",
            w, h, (uint64_t)image_bytes);
    print_disk_free_space(".");

    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        perror("fopen");
        return 0;
    }

    int header_len = fprintf(fp, "P6\n%d %d\n255\n", w, h);
    if (header_len < 0) {
        perror("fprintf");
        fclose(fp);
        return 0;
    }

    size_t header_bytes = (size_t)header_len;

    const size_t CHUNK = 64u * 1024u * 1024u; /* 64 MB */
    size_t total_written = 0;

    while (total_written < image_bytes) {
        size_t remain = image_bytes - total_written;
        size_t this_chunk = remain < CHUNK ? remain : CHUNK;

        size_t written = fwrite(img + total_written, 1, this_chunk, fp);
        if (written != this_chunk) {
            fprintf(stderr,
                    "ERROR: fwrite failed at offset %" PRIu64
                    ", wrote %zu / %zu bytes\n",
                    (uint64_t)total_written, written, this_chunk);
            if (ferror(fp)) {
                perror("fwrite");
            }
            fclose(fp);
            return 0;
        }

        total_written += written;
    }

    if (fflush(fp) != 0) {
        perror("fflush");
        fclose(fp);
        return 0;
    }

    if (fclose(fp) != 0) {
        perror("fclose");
        return 0;
    }

    if (!verify_file_size(filename, header_bytes, image_bytes)) {
        return 0;
    }

    fprintf(stderr, "PPM write OK: %s\n", filename);
    return 1;
}

static void print_usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Options:\n"
        "  --backend=cpu|gpu       choose backend (default: cpu)\n"
        "  --method=newton|secant  choose method  (default: newton)\n"
        "  --output=FILE           output ppm file\n"
        "  --width=N               image width    (default: 2000)\n"
        "  --height=N              image height   (default: 2000)\n"
        "  --max_iter=N            max iterations (default: 50)\n"
        "  --tol=VALUE             tolerance      (default: 1e-6)\n"
        "  --threads=N             CPU threads; 0 means max (default: 0)\n"
        "  --xmin=VALUE            viewport xmin  (default: -2.0)\n"
        "  --xmax=VALUE            viewport xmax  (default:  2.0)\n"
        "  --ymin=VALUE            viewport ymin  (default: -2.0)\n"
        "  --ymax=VALUE            viewport ymax  (default:  2.0)\n",
        prog
    );
}

int main(int argc, char *argv[]) {
    FractalConfig cfg;

    cfg.width = 2000;
    cfg.height = 2000;
    cfg.max_iter = 50;
    cfg.tol = 1e-6;
    cfg.xmin = -2.0;
    cfg.xmax =  2.0;
    cfg.ymin = -2.0;
    cfg.ymax =  2.0;
    cfg.method = METHOD_NEWTON;
    cfg.threads = 0;

    const char *backend = "cpu";
    const char *outfile = "fractal.ppm";

    for (int i = 1; i < argc; ++i) {
        if (strncmp(argv[i], "--backend=", 10) == 0) {
            backend = argv[i] + 10;
        } else if (strncmp(argv[i], "--method=", 9) == 0) {
            const char *m = argv[i] + 9;
            if (strcmp(m, "newton") == 0) cfg.method = METHOD_NEWTON;
            else if (strcmp(m, "secant") == 0) cfg.method = METHOD_SECANT;
            else {
                fprintf(stderr, "Unknown method: %s\n", m);
                return 1;
            }
        } else if (strncmp(argv[i], "--output=", 9) == 0) {
            outfile = argv[i] + 9;
        } else if (strncmp(argv[i], "--width=", 8) == 0) {
            cfg.width = atoi(argv[i] + 8);
        } else if (strncmp(argv[i], "--height=", 9) == 0) {
            cfg.height = atoi(argv[i] + 9);
        } else if (strncmp(argv[i], "--max_iter=", 11) == 0) {
            cfg.max_iter = atoi(argv[i] + 11);
        } else if (strncmp(argv[i], "--tol=", 6) == 0) {
            cfg.tol = atof(argv[i] + 6);
        } else if (strncmp(argv[i], "--threads=", 10) == 0) {
            cfg.threads = atoi(argv[i] + 10);
        } else if (strncmp(argv[i], "--xmin=", 7) == 0) {
            cfg.xmin = atof(argv[i] + 7);
        } else if (strncmp(argv[i], "--xmax=", 7) == 0) {
            cfg.xmax = atof(argv[i] + 7);
        } else if (strncmp(argv[i], "--ymin=", 7) == 0) {
            cfg.ymin = atof(argv[i] + 7);
        } else if (strncmp(argv[i], "--ymax=", 7) == 0) {
            cfg.ymax = atof(argv[i] + 7);
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
    }

    if (cfg.width <= 0 || cfg.height <= 0 || cfg.max_iter <= 0 || cfg.tol <= 0.0) {
        fprintf(stderr, "Invalid numeric arguments.\n");
        return 1;
    }

    size_t bytes = 0;
    if (!compute_image_bytes(cfg.width, cfg.height, &bytes)) {
        fprintf(stderr, "ERROR: image buffer size overflow for %d x %d\n",
                cfg.width, cfg.height);
        return 1;
    }

    fprintf(stderr, "Allocating image buffer: %" PRIu64 " bytes (%.2f GB)\n",
            (uint64_t)bytes, (double)bytes / (1024.0 * 1024.0 * 1024.0));

    unsigned char *img = (unsigned char *)malloc(bytes);
    if (!img) {
        fprintf(stderr, "Failed to allocate image buffer (%zu bytes)\n", bytes);
        return 1;
    }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int used_gpu = 0;

    if (strcmp(backend, "gpu") == 0) {
        int rc = render_gpu(img, &cfg);
        if (rc == 0) {
            used_gpu = 1;
        } else {
            fprintf(stderr, "[warn] GPU backend unavailable or failed, fallback to CPU.\n");
            render_cpu(img, &cfg);
        }
    } else {
        render_cpu(img, &cfg);
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);

    if (!write_ppm_checked(outfile, img, cfg.width, cfg.height)) {
        fprintf(stderr, "ERROR: failed to write complete PPM file.\n");
        free(img);
        return 1;
    }

    double sec = elapsed_seconds(t0, t1);
    fprintf(stderr, "Output      : %s\n", outfile);
    fprintf(stderr, "Backend     : %s\n", used_gpu ? "gpu" : "cpu");
    fprintf(stderr, "Method      : %s\n", cfg.method == METHOD_NEWTON ? "newton" : "secant");
    fprintf(stderr, "Resolution  : %d x %d\n", cfg.width, cfg.height);
    fprintf(stderr, "Max iter    : %d\n", cfg.max_iter);
    fprintf(stderr, "Tolerance   : %.3e\n", cfg.tol);
    fprintf(stderr, "Compute time: %.6f seconds\n", sec);

    free(img);
    return 0;
}
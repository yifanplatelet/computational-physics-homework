#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "cpu_backend.h"
#include "cuda_backend.h"

static double elapsed_seconds(struct timespec a, struct timespec b) {
    return (double)(b.tv_sec - a.tv_sec) +
           (double)(b.tv_nsec - a.tv_nsec) / 1e9;
}

static void write_ppm(const char *filename, const unsigned char *img, int w, int h) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        perror("fopen");
        exit(1);
    }

    fprintf(fp, "P6\n%d %d\n255\n", w, h);
    fwrite(img, 1, (size_t)(3 * w * h), fp);
    fclose(fp);
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
        "  --symmetry=off|auto|on  conjugate-symmetry optimization\n"
        "                         (default: off; use auto/on only for symmetric windows)\n"
        "  --xmin=VALUE            viewport xmin  (default: -2.0)\n"
        "  --xmax=VALUE            viewport xmax  (default:  2.0)\n"
        "  --ymin=VALUE            viewport ymin  (default: -2.0)\n"
        "  --ymax=VALUE            viewport ymax  (default:  2.0)\n"
        "\n"
        "Examples:\n"
        "  %s --backend=cpu --method=newton --output=newton.ppm\n"
        "  %s --backend=gpu --method=secant --output=secant_gpu.ppm\n",
        prog, prog, prog
    );
}

int main(int argc, char *argv[]) {
    FractalConfig cfg;

    cfg.width = 8000;
    cfg.height = 8000;
    cfg.max_iter = 50;
    cfg.tol = 1e-6;
    cfg.xmin = -2.0;
    cfg.xmax =  2.0;
    cfg.ymin = -2.0;
    cfg.ymax =  2.0;
    cfg.method = METHOD_NEWTON;
    cfg.threads = 0;
    cfg.symmetry_mode = SYMMETRY_OFF;

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
        } else if (strncmp(argv[i], "--symmetry=", 11) == 0) {
            const char *s = argv[i] + 11;
            if (strcmp(s, "off") == 0) cfg.symmetry_mode = SYMMETRY_OFF;
            else if (strcmp(s, "auto") == 0) cfg.symmetry_mode = SYMMETRY_AUTO;
            else if (strcmp(s, "on") == 0) cfg.symmetry_mode = SYMMETRY_ON;
            else {
                fprintf(stderr, "Unknown symmetry mode: %s\n", s);
                return 1;
            }
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

    size_t bytes = (size_t)3 * (size_t)cfg.width * (size_t)cfg.height;
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

    write_ppm(outfile, img, cfg.width, cfg.height);

    double sec = elapsed_seconds(t0, t1);
    fprintf(stderr, "Output      : %s\n", outfile);
    fprintf(stderr, "Backend     : %s\n", used_gpu ? "gpu" : "cpu");
    fprintf(stderr, "Method      : %s\n", cfg.method == METHOD_NEWTON ? "newton" : "secant");
    fprintf(stderr, "Resolution  : %d x %d\n", cfg.width, cfg.height);
    fprintf(stderr, "Max iter    : %d\n", cfg.max_iter);
    fprintf(stderr, "Tolerance   : %.3e\n", cfg.tol);
    fprintf(stderr, "Symmetry    : %s\n",
            cfg.symmetry_mode == SYMMETRY_OFF ? "off" :
            (cfg.symmetry_mode == SYMMETRY_AUTO ? "auto" : "on"));
    fprintf(stderr, "Elapsed     : %.6f seconds\n", sec);

    free(img);
    return 0;
}

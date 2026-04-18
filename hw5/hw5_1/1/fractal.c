#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <time.h>

typedef struct {
    unsigned char r, g, b;
} Color;

typedef struct {
    int converged;   // 1 = converged, 0 = failed
    int root_id;     // 0 -> root 1, 1/2 -> complex roots
    int iter;
} Result;

static inline double complex f(double complex z) {
    return z*z*z - 1.0;
}

static inline double complex fp(double complex z) {
    return 3.0*z*z;
}

static int nearest_root(double complex z) {
    const double complex roots[3] = {
        1.0 + 0.0*I,
        -0.5 + (sqrt(3.0)/2.0)*I,
        -0.5 - (sqrt(3.0)/2.0)*I
    };

    double d0 = cabs(z - roots[0]);
    double d1 = cabs(z - roots[1]);
    double d2 = cabs(z - roots[2]);

    if (d0 <= d1 && d0 <= d2) return 0;
    if (d1 <= d0 && d1 <= d2) return 1;
    return 2;
}

Result iterate_newton(double complex z0, int max_iter, double tol) {
    Result res;
    res.converged = 0;
    res.root_id = -1;
    res.iter = max_iter;

    double complex z = z0;

    for (int k = 0; k < max_iter; ++k) {
        double complex fv = f(z);
        double complex df = fp(z);

        if (cabs(df) < 1e-14) {
            return res;
        }

        double complex z_next = z - fv / df;

        if (cabs(z_next - z) < tol || cabs(f(z_next)) < tol) {
            res.converged = 1;
            res.root_id = nearest_root(z_next);
            res.iter = k + 1;
            return res;
        }

        z = z_next;
    }

    return res;
}

Result iterate_secant(double complex z0, int max_iter, double tol, double complex delta) {
    Result res;
    res.converged = 0;
    res.root_id = -1;
    res.iter = max_iter;

    double complex z_prev = z0;
    double complex z_cur  = z0 + delta;

    for (int k = 0; k < max_iter; ++k) {
        double complex f_prev = f(z_prev);
        double complex f_cur  = f(z_cur);
        double complex denom  = f_cur - f_prev;

        if (cabs(denom) < 1e-14) {
            return res;
        }

        double complex z_next = z_cur - f_cur * (z_cur - z_prev) / denom;

        if (cabs(z_next - z_cur) < tol || cabs(f(z_next)) < tol) {
            res.converged = 1;
            res.root_id = nearest_root(z_next);
            res.iter = k + 1;
            return res;
        }

        z_prev = z_cur;
        z_cur  = z_next;
    }

    return res;
}

static inline unsigned char clamp_u8(double x) {
    if (x < 0.0) return 0;
    if (x > 255.0) return 255;
    return (unsigned char)(x + 0.5);
}

Color shade_two_color(int root_id, int iter, int max_iter, int converged) {
    Color A = {255, 210, 60};
    Color B = {70, 180, 255};

    if (!converged) {
        Color black = {0, 0, 0};
        return black;
    }

    double t = (double)iter / (double)max_iter;
    double brightness = 1.0 - 0.82 * t;

    Color base = (root_id == 0) ? A : B;

    Color out;
    out.r = clamp_u8(base.r * brightness);
    out.g = clamp_u8(base.g * brightness);
    out.b = clamp_u8(base.b * brightness);
    return out;
}

static double elapsed_seconds(struct timespec a, struct timespec b) {
    return (double)(b.tv_sec - a.tv_sec) +
           (double)(b.tv_nsec - a.tv_nsec) / 1e9;
}

int main(int argc, char *argv[]) {
    int width = 20000;
    int height = 20000;
    int max_iter = 50;
    double tol = 1e-5;

    double xmin = -2.0, xmax = 2.0;
    double ymin = -2.0, ymax = 2.0;

    const char *method = "newton";
    if (argc >= 2) {
        method = argv[1];
    }

    const char *outfile = NULL;
    if (argc >= 3) {
        outfile = argv[2];
    } else {
        outfile = (strcmp(method, "secant") == 0) ? "secant_2color.ppm" : "newton_2color.ppm";
    }

    FILE *fp_out = fopen(outfile, "wb");
    if (!fp_out) {
        perror("fopen");
        return 1;
    }

    fprintf(fp_out, "P6\n%d %d\n255\n", width, height);

    double complex delta = 1e-3 + 1e-3*I;

    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    for (int j = 0; j < height; ++j) {
        double y = ymax - (ymax - ymin) * j / (height - 1.0);

        for (int i = 0; i < width; ++i) {
            double x = xmin + (xmax - xmin) * i / (width - 1.0);
            double complex z0 = x + y * I;

            Result res;
            if (strcmp(method, "secant") == 0) {
                res = iterate_secant(z0, max_iter, tol, delta);
            } else {
                res = iterate_newton(z0, max_iter, tol);
            }

            Color c = shade_two_color(res.root_id, res.iter, max_iter, res.converged);
            fputc(c.r, fp_out);
            fputc(c.g, fp_out);
            fputc(c.b, fp_out);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = elapsed_seconds(t_start, t_end);

    fclose(fp_out);

    fprintf(stderr, "Done. Output: %s\n", outfile);
    fprintf(stderr, "Method=%s, size=%dx%d, max_iter=%d, tol=%g\n",
            method, width, height, max_iter, tol);
    fprintf(stderr, "Elapsed time: %.6f seconds\n", elapsed);

    return 0;
}
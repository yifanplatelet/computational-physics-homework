#include "cpu_backend.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    int converged;
    int root_id;
    int iter;
} Result;

/* 三个根 */
static const double ROOT0_X =  1.0;
static const double ROOT0_Y =  0.0;

static const double ROOT1_X = -0.5;
static const double ROOT1_Y =  0.86602540378443864676; /* sqrt(3)/2 */

static const double ROOT2_X = -0.5;
static const double ROOT2_Y = -0.86602540378443864676; /* -sqrt(3)/2 */

static inline unsigned char clamp_u8(double x) {
    if (x < 0.0) return 0;
    if (x > 255.0) return 255;
    return (unsigned char)(x + 0.5);
}

static inline double dist2(double x1, double y1, double x2, double y2) {
    double dx = x1 - x2;
    double dy = y1 - y2;
    return dx * dx + dy * dy;
}

/* 若已接近某根，直接返回 root_id；否则返回 -1 */
static inline int detect_root_if_close(double x, double y, double tol2) {
    if (dist2(x, y, ROOT0_X, ROOT0_Y) < tol2) return 0;
    if (dist2(x, y, ROOT1_X, ROOT1_Y) < tol2) return 1;
    if (dist2(x, y, ROOT2_X, ROOT2_Y) < tol2) return 2;
    return -1;
}

/* 仅作兜底分类 */
static inline int nearest_root_xy(double x, double y) {
    double d0 = dist2(x, y, ROOT0_X, ROOT0_Y);
    double d1 = dist2(x, y, ROOT1_X, ROOT1_Y);
    double d2 = dist2(x, y, ROOT2_X, ROOT2_Y);

    if (d0 <= d1 && d0 <= d2) return 0;
    if (d1 <= d0 && d1 <= d2) return 1;
    return 2;
}

/*
 * 优化版 Newton：
 * 对 z^3 - 1 = 0，Newton 可化简为
 *
 *   z_{n+1} = (2/3) z + 1 / (3 z^2)
 *
 * 这样避免每次显式计算 f(z)、f'(z) 和一次完整复除法。
 */
static inline Result iterate_newton_xy_optimized(double x0, double y0,
                                                 int max_iter, double tol) {
    Result res;
    res.converged = 0;
    res.root_id = -1;
    res.iter = max_iter;

    double x = x0;
    double y = y0;

    const double tol2 = tol * tol;
    const double two_thirds = 2.0 / 3.0;
    const double one_third  = 1.0 / 3.0;

    for (int k = 0; k < max_iter; ++k) {
        /* 先判断当前点是否已经很接近某个根 */
        int rid = detect_root_if_close(x, y, tol2);
        if (rid >= 0) {
            res.converged = 1;
            res.root_id = rid;
            res.iter = k;
            return res;
        }

        /* z^2 = (a + bi) */
        double a = x * x - y * y;
        double b = 2.0 * x * y;

        /* 1 / z^2 = (a - bi) / (a^2 + b^2) */
        double denom = a * a + b * b;
        if (denom < 1e-28) {
            return res;
        }

        double inv_r =  a / denom;
        double inv_i = -b / denom;

        /* z_next = (2/3) z + (1/3) * (1 / z^2) */
        double xn = two_thirds * x + one_third * inv_r;
        double yn = two_thirds * y + one_third * inv_i;

        /* 再判断新点是否接近某根 */
        rid = detect_root_if_close(xn, yn, tol2);
        if (rid >= 0) {
            res.converged = 1;
            res.root_id = rid;
            res.iter = k + 1;
            return res;
        }

        /* 步长判据 */
        double dx = xn - x;
        double dy = yn - y;
        if (dx * dx + dy * dy < tol2) {
            res.converged = 1;
            res.root_id = nearest_root_xy(xn, yn);
            res.iter = k + 1;
            return res;
        }

        x = xn;
        y = yn;
    }

    return res;
}

/*
 * Secant 仍保留原来的专用实部/虚部写法。
 * 这里也加了接近根提前退出。
 */
static inline Result iterate_secant_xy_optimized(double x0, double y0,
                                                 int max_iter, double tol) {
    Result res;
    res.converged = 0;
    res.root_id = -1;
    res.iter = max_iter;

    double x_prev = x0;
    double y_prev = y0;
    double x_cur  = x0 + 1e-3;
    double y_cur  = y0 + 1e-3;

    const double tol2 = tol * tol;

    for (int k = 0; k < max_iter; ++k) {
        int rid = detect_root_if_close(x_cur, y_cur, tol2);
        if (rid >= 0) {
            res.converged = 1;
            res.root_id = rid;
            res.iter = k;
            return res;
        }

        /* f(z_prev) = z^3 - 1 */
        double z2r = x_prev * x_prev - y_prev * y_prev;
        double z2i = 2.0 * x_prev * y_prev;
        double fpr = z2r * x_prev - z2i * y_prev - 1.0;
        double fpi = z2r * y_prev + z2i * x_prev;

        /* f(z_cur) = z^3 - 1 */
        z2r = x_cur * x_cur - y_cur * y_cur;
        z2i = 2.0 * x_cur * y_cur;
        double fcr = z2r * x_cur - z2i * y_cur - 1.0;
        double fci = z2r * y_cur + z2i * x_cur;

        /* denom = f_cur - f_prev */
        double dr = fcr - fpr;
        double di = fci - fpi;
        double denom = dr * dr + di * di;
        if (denom < 1e-28) {
            return res;
        }

        /* (z_cur - z_prev) / (f_cur - f_prev) */
        double zr = x_cur - x_prev;
        double zi = y_cur - y_prev;
        double rr = (zr * dr + zi * di) / denom;
        double ri = (zi * dr - zr * di) / denom;

        /* f_cur * ratio */
        double mr = fcr * rr - fci * ri;
        double mi = fcr * ri + fci * rr;

        double x_next = x_cur - mr;
        double y_next = y_cur - mi;

        rid = detect_root_if_close(x_next, y_next, tol2);
        if (rid >= 0) {
            res.converged = 1;
            res.root_id = rid;
            res.iter = k + 1;
            return res;
        }

        double dx = x_next - x_cur;
        double dy = y_next - y_cur;
        if (dx * dx + dy * dy < tol2) {
            res.converged = 1;
            res.root_id = nearest_root_xy(x_next, y_next);
            res.iter = k + 1;
            return res;
        }

        x_prev = x_cur;
        y_prev = y_cur;
        x_cur = x_next;
        y_cur = y_next;
    }

    return res;
}

static inline void shade_two_color(unsigned char *r,
                                   unsigned char *g,
                                   unsigned char *b,
                                   int root_id,
                                   int iter,
                                   int max_iter,
                                   int converged) {
    if (!converged) {
        *r = 0;
        *g = 0;
        *b = 0;
        return;
    }

    const double Ar = 255.0, Ag = 210.0, Ab =  60.0;  /* 实根 */
    const double Br =  70.0, Bg = 180.0, Bb = 255.0;  /* 两个复根 */

    double t = (double)iter / (double)max_iter;
    double brightness = 1.0 - 0.82 * t;

    double base_r = (root_id == 0) ? Ar : Br;
    double base_g = (root_id == 0) ? Ag : Bg;
    double base_b = (root_id == 0) ? Ab : Bb;

    *r = clamp_u8(base_r * brightness);
    *g = clamp_u8(base_g * brightness);
    *b = clamp_u8(base_b * brightness);
}

void render_cpu(unsigned char *img, const FractalConfig *cfg) {
    const int width = cfg->width;
    const int height = cfg->height;

    const double xmin = cfg->xmin;
    const double xmax = cfg->xmax;
    const double ymin = cfg->ymin;
    const double ymax = cfg->ymax;

    const double dx = (xmax - xmin) / (double)(width - 1);
    const double dy = (ymax - ymin) / (double)(height - 1);

#ifdef _OPENMP
    int threads = cfg->threads > 0 ? cfg->threads : omp_get_max_threads();
    omp_set_num_threads(threads);
    fprintf(stderr, "[cpu] OpenMP threads = %d\n", threads);

    #pragma omp parallel for schedule(guided)
#endif
    for (int j = 0; j < height; ++j) {
        double y = ymax - j * dy;

        for (int i = 0; i < width; ++i) {
            double x = xmin + i * dx;

            Result res;
            if (cfg->method == METHOD_SECANT) {
                res = iterate_secant_xy_optimized(x, y, cfg->max_iter, cfg->tol);
            } else {
                res = iterate_newton_xy_optimized(x, y, cfg->max_iter, cfg->tol);
            }

            unsigned char r, g, b;
            shade_two_color(&r, &g, &b,
                            res.root_id, res.iter, cfg->max_iter, res.converged);

            size_t idx = (size_t)3 * ((size_t)j * (size_t)width + (size_t)i);
            img[idx + 0] = r;
            img[idx + 1] = g;
            img[idx + 2] = b;
        }
    }
}
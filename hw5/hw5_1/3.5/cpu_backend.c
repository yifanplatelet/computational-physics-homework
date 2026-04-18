#include "cpu_backend.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __AVX__
#include <immintrin.h>
#endif

typedef struct {
    int converged;
    int root_id;
    int iter;
} Result;

static inline unsigned char clamp_u8(double x) {
    if (x < 0.0) return 0;
    if (x > 255.0) return 255;
    return (unsigned char)(x + 0.5);
}

static inline int nearest_root_xy(double x, double y) {
    const double s = 0.86602540378443864676; /* sqrt(3) / 2 */

    double d0 = (x - 1.0) * (x - 1.0) + y * y;
    double xp = x + 0.5;
    double d1 = xp * xp + (y - s) * (y - s);
    double d2 = xp * xp + (y + s) * (y + s);

    if (d0 < d1) return (d0 < d2) ? 0 : 2;
    return (d1 < d2) ? 1 : 2;
}

static inline Result iterate_newton_xy(double x0, double y0, int max_iter, double tol) {
    const double one_third = 1.0 / 3.0;
    const double singular_r2 = 1e-14;

    Result res;
    res.converged = 0;
    res.root_id = -1;
    res.iter = max_iter;

    double x = x0;
    double y = y0;
    double tol2 = tol * tol;

    for (int k = 0; k < max_iter; ++k) {
        double xx = x * x;
        double yy = y * y;
        double xy = x * y;
        double r2 = xx + yy;

        if (r2 < singular_r2) {
            return res;
        }

        double inv = one_third / (r2 * r2);
        double dx = -one_third * x + (xx - yy) * inv;
        double dy = -one_third * y - (xy + xy) * inv;

        x += dx;
        y += dy;

        if (dx * dx + dy * dy < tol2) {
            res.converged = 1;
            res.root_id = nearest_root_xy(x, y);
            res.iter = k + 1;
            return res;
        }
    }

    return res;
}

static inline Result iterate_secant_xy(double x0, double y0, int max_iter, double tol) {
    Result res;
    res.converged = 0;
    res.root_id = -1;
    res.iter = max_iter;

    /* Use a purely real perturbation so conjugate symmetry is preserved. */
    double x_prev = x0;
    double y_prev = y0;
    double x_cur  = x0 + 1e-3;
    double y_cur  = y0;

    double tol2 = tol * tol;

    for (int k = 0; k < max_iter; ++k) {
        double z2r = x_prev * x_prev - y_prev * y_prev;
        double z2i = 2.0 * x_prev * y_prev;
        double fpr = z2r * x_prev - z2i * y_prev - 1.0;
        double fpi = z2r * y_prev + z2i * x_prev;

        z2r = x_cur * x_cur - y_cur * y_cur;
        z2i = 2.0 * x_cur * y_cur;
        double fcr = z2r * x_cur - z2i * y_cur - 1.0;
        double fci = z2r * y_cur + z2i * x_cur;

        double dr = fcr - fpr;
        double di = fci - fpi;
        double denom = dr * dr + di * di;

        if (denom < 1e-28) {
            return res;
        }

        double zr = x_cur - x_prev;
        double zi = y_cur - y_prev;
        double rr = (zr * dr + zi * di) / denom;
        double ri = (zi * dr - zr * di) / denom;

        double mr = fcr * rr - fci * ri;
        double mi = fcr * ri + fci * rr;

        double dx = -mr;
        double dy = -mi;

        double x_next = x_cur + dx;
        double y_next = y_cur + dy;

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
        *r = 0; *g = 0; *b = 0;
        return;
    }

    const double Ar = 255.0, Ag = 210.0, Ab =  60.0;
    const double Br =  70.0, Bg = 180.0, Bb = 255.0;

    double t = (double)iter / (double)max_iter;
    double brightness = 1.0 - 0.82 * t;

    double base_r = (root_id == 0) ? Ar : Br;
    double base_g = (root_id == 0) ? Ag : Bg;
    double base_b = (root_id == 0) ? Ab : Bb;

    *r = clamp_u8(base_r * brightness);
    *g = clamp_u8(base_g * brightness);
    *b = clamp_u8(base_b * brightness);
}

static inline void shade_result(unsigned char *dst, const Result *res, int max_iter) {
    shade_two_color(&dst[0], &dst[1], &dst[2],
                    res->root_id, res->iter, max_iter, res->converged);
}

static inline void render_row_newton_scalar(unsigned char *row,
                                            const double *xcoords,
                                            double y,
                                            int width,
                                            int max_iter,
                                            double tol) {
    for (int i = 0; i < width; ++i) {
        Result res = iterate_newton_xy(xcoords[i], y, max_iter, tol);
        shade_result(row + (size_t)3 * (size_t)i, &res, max_iter);
    }
}

static inline void render_row_secant_scalar(unsigned char *row,
                                            const double *xcoords,
                                            double y,
                                            int width,
                                            int max_iter,
                                            double tol) {
    for (int i = 0; i < width; ++i) {
        Result res = iterate_secant_xy(xcoords[i], y, max_iter, tol);
        shade_result(row + (size_t)3 * (size_t)i, &res, max_iter);
    }
}

#ifdef __AVX__
static inline __m256d lane_mask_pd(int mask) {
    return _mm256_setr_pd((mask & 0x1) ? -1.0 : 0.0,
                          (mask & 0x2) ? -1.0 : 0.0,
                          (mask & 0x4) ? -1.0 : 0.0,
                          (mask & 0x8) ? -1.0 : 0.0);
}

static inline void render_row_newton_simd(unsigned char *row,
                                          const double *xcoords,
                                          double y,
                                          int width,
                                          int max_iter,
                                          double tol) {
    const __m256d one_third = _mm256_set1_pd(1.0 / 3.0);
    const __m256d tol2v = _mm256_set1_pd(tol * tol);
    const __m256d singular_r2v = _mm256_set1_pd(1e-14);
    const __m256d onev = _mm256_set1_pd(1.0);
    const __m256d y0v = _mm256_set1_pd(y);

    int i = 0;
    for (; i + 4 <= width; i += 4) {
        __m256d x = _mm256_loadu_pd(xcoords + i);
        __m256d yv = y0v;
        int converged[4] = {0, 0, 0, 0};
        int iters[4] = {max_iter, max_iter, max_iter, max_iter};
        int active_mask = 0xF;

        for (int k = 0; k < max_iter && active_mask; ++k) {
            __m256d xx = _mm256_mul_pd(x, x);
            __m256d yy = _mm256_mul_pd(yv, yv);
            __m256d xy = _mm256_mul_pd(x, yv);
            __m256d r2 = _mm256_add_pd(xx, yy);

            int singular_mask = _mm256_movemask_pd(
                _mm256_cmp_pd(r2, singular_r2v, _CMP_LT_OQ)) & active_mask;

            __m256d r2_safe = _mm256_blendv_pd(r2, onev,
                                               _mm256_cmp_pd(r2, singular_r2v, _CMP_LT_OQ));
            __m256d inv = _mm256_div_pd(one_third, _mm256_mul_pd(r2_safe, r2_safe));
            __m256d dx = _mm256_add_pd(
                _mm256_mul_pd(_mm256_sub_pd(xx, yy), inv),
                _mm256_mul_pd(_mm256_set1_pd(-1.0 / 3.0), x));
            __m256d dy = _mm256_sub_pd(
                _mm256_mul_pd(_mm256_set1_pd(-1.0 / 3.0), yv),
                _mm256_mul_pd(_mm256_add_pd(xy, xy), inv));

            __m256d x_next = _mm256_add_pd(x, dx);
            __m256d y_next = _mm256_add_pd(yv, dy);
            __m256d step2 = _mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy));

            int done_mask = _mm256_movemask_pd(
                _mm256_cmp_pd(step2, tol2v, _CMP_LT_OQ)) & active_mask & (~singular_mask);

            __m256d active_v = lane_mask_pd(active_mask);
            x = _mm256_blendv_pd(x, x_next, active_v);
            yv = _mm256_blendv_pd(yv, y_next, active_v);

            int finished_mask = done_mask | singular_mask;
            int newly_converged = done_mask;
            for (int lane = 0; lane < 4; ++lane) {
                if (newly_converged & (1 << lane)) {
                    converged[lane] = 1;
                    iters[lane] = k + 1;
                }
            }

            active_mask &= ~finished_mask;
        }

        double xf[4];
        double yf[4];
        _mm256_storeu_pd(xf, x);
        _mm256_storeu_pd(yf, yv);

        for (int lane = 0; lane < 4; ++lane) {
            Result res;
            res.converged = converged[lane];
            res.root_id = converged[lane] ? nearest_root_xy(xf[lane], yf[lane]) : -1;
            res.iter = iters[lane];
            shade_result(row + (size_t)3 * (size_t)(i + lane), &res, max_iter);
        }
    }

    for (; i < width; ++i) {
        Result res = iterate_newton_xy(xcoords[i], y, max_iter, tol);
        shade_result(row + (size_t)3 * (size_t)i, &res, max_iter);
    }
}
#endif

static inline int has_real_axis_symmetry(const FractalConfig *cfg) {
    double scale = fmax(1.0, fmax(fabs(cfg->ymin), fabs(cfg->ymax)));
    return fabs(cfg->ymin + cfg->ymax) <= 1e-12 * scale;
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

    double *xcoords = (double *)malloc((size_t)width * sizeof(double));
    if (!xcoords) {
        fprintf(stderr, "[cpu] failed to allocate x coordinate buffer\n");
        return;
    }

    for (int i = 0; i < width; ++i) {
        xcoords[i] = xmin + (double)i * dx;
    }

    const int mirror_rows = has_real_axis_symmetry(cfg);
    const int rows_to_compute = mirror_rows ? (height + 1) / 2 : height;

#ifdef _OPENMP
    int threads = cfg->threads > 0 ? cfg->threads : omp_get_max_threads();
    omp_set_num_threads(threads);
    fprintf(stderr, "[cpu] OpenMP threads = %d\n", threads);
    fprintf(stderr, "[cpu] real-axis symmetry = %s\n", mirror_rows ? "on" : "off");
    fprintf(stderr, "[cpu] Newton SIMD (AVX) = %s\n",
#if defined(__AVX__)
            (cfg->method == METHOD_NEWTON) ? "on" : "off"
#else
            "off"
#endif
    );

    #pragma omp parallel for schedule(guided)
#endif
    for (int j = 0; j < rows_to_compute; ++j) {
        double y = ymax - (double)j * dy;
        unsigned char *row_top = img + (size_t)3 * (size_t)j * (size_t)width;

        if (cfg->method == METHOD_SECANT) {
            render_row_secant_scalar(row_top, xcoords, y, width, cfg->max_iter, cfg->tol);
        } else {
#ifdef __AVX__
            render_row_newton_simd(row_top, xcoords, y, width, cfg->max_iter, cfg->tol);
#else
            render_row_newton_scalar(row_top, xcoords, y, width, cfg->max_iter, cfg->tol);
#endif
        }

        if (mirror_rows) {
            int mj = height - 1 - j;
            if (mj != j) {
                unsigned char *row_bottom = img + (size_t)3 * (size_t)mj * (size_t)width;
                memcpy(row_bottom, row_top, (size_t)3 * (size_t)width);
            }
        }
    }

    free(xcoords);
}

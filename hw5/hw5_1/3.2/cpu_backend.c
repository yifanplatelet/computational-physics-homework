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

static inline unsigned char clamp_u8(double x) {
    if (x < 0.0) return 0;
    if (x > 255.0) return 255;
    return (unsigned char)(x + 0.5);
}

static inline int nearest_root_xy(double x, double y) {
    const double rx0 = 1.0,  ry0 = 0.0;
    const double rx1 = -0.5, ry1 =  sqrt(3.0) / 2.0;
    const double rx2 = -0.5, ry2 = -sqrt(3.0) / 2.0;

    double d0 = (x - rx0)*(x - rx0) + (y - ry0)*(y - ry0);
    double d1 = (x - rx1)*(x - rx1) + (y - ry1)*(y - ry1);
    double d2 = (x - rx2)*(x - rx2) + (y - ry2)*(y - ry2);

    if (d0 <= d1 && d0 <= d2) return 0;
    if (d1 <= d0 && d1 <= d2) return 1;
    return 2;
}

static inline Result iterate_newton_xy(double x0, double y0, int max_iter, double tol) {
    Result res;
    res.converged = 0;
    res.root_id = -1;
    res.iter = max_iter;

    double x = x0, y = y0;
    double tol2 = tol * tol;

    for (int k = 0; k < max_iter; ++k) {
        /* z^2 */
        double z2r = x*x - y*y;
        double z2i = 2.0*x*y;

        /* z^3 */
        double z3r = z2r*x - z2i*y;
        double z3i = z2r*y + z2i*x;

        /* f(z)=z^3-1 */
        double fr = z3r - 1.0;
        double fi = z3i;

        /* f'(z)=3z^2 */
        double dfr = 3.0 * z2r;
        double dfi = 3.0 * z2i;

        double denom = dfr*dfr + dfi*dfi;
        if (denom < 1e-28) {
            return res;
        }

        /* f/df */
        double qr = (fr*dfr + fi*dfi) / denom;
        double qi = (fi*dfr - fr*dfi) / denom;

        double xn = x - qr;
        double yn = y - qi;

        double dx = xn - x;
        double dy = yn - y;

        /* 也可用 f(z_next) 判断，这里先用步长 */
        if (dx*dx + dy*dy < tol2) {
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

static inline Result iterate_secant_xy(double x0, double y0, int max_iter, double tol) {
    Result res;
    res.converged = 0;
    res.root_id = -1;
    res.iter = max_iter;

    /* z_prev = z0, z_cur = z0 + delta */
    double x_prev = x0;
    double y_prev = y0;
    double x_cur  = x0 + 1e-3;
    double y_cur  = y0 + 1e-3;

    double tol2 = tol * tol;

    for (int k = 0; k < max_iter; ++k) {
        /* f(z_prev) = z^3 - 1 */
        double z2r = x_prev*x_prev - y_prev*y_prev;
        double z2i = 2.0*x_prev*y_prev;
        double fpr = z2r*x_prev - z2i*y_prev - 1.0;
        double fpi = z2r*y_prev + z2i*x_prev;

        /* f(z_cur) = z^3 - 1 */
        z2r = x_cur*x_cur - y_cur*y_cur;
        z2i = 2.0*x_cur*y_cur;
        double fcr = z2r*x_cur - z2i*y_cur - 1.0;
        double fci = z2r*y_cur + z2i*x_cur;

        /* denom = f_cur - f_prev */
        double dr = fcr - fpr;
        double di = fci - fpi;
        double denom = dr*dr + di*di;

        if (denom < 1e-28) {
            return res;
        }

        /* (z_cur - z_prev) / (f_cur - f_prev) */
        double zr = x_cur - x_prev;
        double zi = y_cur - y_prev;
        double rr = (zr*dr + zi*di) / denom;
        double ri = (zi*dr - zr*di) / denom;

        /* f_cur * ratio */
        double mr = fcr*rr - fci*ri;
        double mi = fcr*ri + fci*rr;

        double x_next = x_cur - mr;
        double y_next = y_cur - mi;

        double dx = x_next - x_cur;
        double dy = y_next - y_cur;

        if (dx*dx + dy*dy < tol2) {
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

    const double Ar = 255.0, Ag = 210.0, Ab =  60.0;  /* 实根 */
    const double Br =  70.0, Bg = 180.0, Bb = 255.0;  /* 复根 */

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

    #pragma omp parallel for schedule(dynamic)
#endif
    for (int j = 0; j < height; ++j) {
        double y = ymax - j * dy;

        for (int i = 0; i < width; ++i) {
            double x = xmin + i * dx;

            Result res;
            if (cfg->method == METHOD_SECANT) {
                res = iterate_secant_xy(x, y, cfg->max_iter, cfg->tol);
            } else {
                res = iterate_newton_xy(x, y, cfg->max_iter, cfg->tol);
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
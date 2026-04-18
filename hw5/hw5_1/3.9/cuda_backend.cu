#include "cuda_backend.h"

#include <stdio.h>
#include <cuda_runtime.h>

typedef struct {
    int converged;
    int root_id;
    int iter;
} ResultGPU;

__device__ static inline unsigned char clamp_u8_gpu(double x) {
    if (x < 0.0) return 0;
    if (x > 255.0) return 255;
    return (unsigned char)(x + 0.5);
}

__device__ static inline int nearest_root_xy_gpu(double x, double y) {
    const double rx0 = 1.0,  ry0 = 0.0;
    const double rx1 = -0.5, ry1 =  0.8660254037844386;
    const double rx2 = -0.5, ry2 = -0.8660254037844386;

    double d0 = (x - rx0)*(x - rx0) + (y - ry0)*(y - ry0);
    double d1 = (x - rx1)*(x - rx1) + (y - ry1)*(y - ry1);
    double d2 = (x - rx2)*(x - rx2) + (y - ry2)*(y - ry2);

    if (d0 <= d1 && d0 <= d2) return 0;
    if (d1 <= d0 && d1 <= d2) return 1;
    return 2;
}

__device__ static inline ResultGPU iterate_newton_xy_gpu(double x0, double y0,
                                                         int max_iter, double tol) {
    ResultGPU res;
    res.converged = 0;
    res.root_id = -1;
    res.iter = max_iter;

    double x = x0, y = y0;
    double tol2 = tol * tol;

    for (int k = 0; k < max_iter; ++k) {
        double z2r = x*x - y*y;
        double z2i = 2.0*x*y;

        double z3r = z2r*x - z2i*y;
        double z3i = z2r*y + z2i*x;

        double fr = z3r - 1.0;
        double fi = z3i;

        double dfr = 3.0 * z2r;
        double dfi = 3.0 * z2i;

        double denom = dfr*dfr + dfi*dfi;
        if (denom < 1e-28) return res;

        double qr = (fr*dfr + fi*dfi) / denom;
        double qi = (fi*dfr - fr*dfi) / denom;

        double xn = x - qr;
        double yn = y - qi;

        double dx = xn - x;
        double dy = yn - y;

        if (dx*dx + dy*dy < tol2) {
            res.converged = 1;
            res.root_id = nearest_root_xy_gpu(xn, yn);
            res.iter = k + 1;
            return res;
        }

        x = xn;
        y = yn;
    }

    return res;
}

__device__ static inline ResultGPU iterate_secant_xy_gpu(double x0, double y0,
                                                         int max_iter, double tol) {
    ResultGPU res;
    res.converged = 0;
    res.root_id = -1;
    res.iter = max_iter;

    double x_prev = x0;
    double y_prev = y0;
    double x_cur  = x0 + 1e-3;
    double y_cur  = y0 + 1e-3;

    double tol2 = tol * tol;

    for (int k = 0; k < max_iter; ++k) {
        double z2r = x_prev*x_prev - y_prev*y_prev;
        double z2i = 2.0*x_prev*y_prev;
        double fpr = z2r*x_prev - z2i*y_prev - 1.0;
        double fpi = z2r*y_prev + z2i*x_prev;

        z2r = x_cur*x_cur - y_cur*y_cur;
        z2i = 2.0*x_cur*y_cur;
        double fcr = z2r*x_cur - z2i*y_cur - 1.0;
        double fci = z2r*y_cur + z2i*x_cur;

        double dr = fcr - fpr;
        double di = fci - fpi;
        double denom = dr*dr + di*di;
        if (denom < 1e-28) return res;

        double zr = x_cur - x_prev;
        double zi = y_cur - y_prev;
        double rr = (zr*dr + zi*di) / denom;
        double ri = (zi*dr - zr*di) / denom;

        double mr = fcr*rr - fci*ri;
        double mi = fcr*ri + fci*rr;

        double x_next = x_cur - mr;
        double y_next = y_cur - mi;

        double dx = x_next - x_cur;
        double dy = y_next - y_cur;

        if (dx*dx + dy*dy < tol2) {
            res.converged = 1;
            res.root_id = nearest_root_xy_gpu(x_next, y_next);
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

__device__ static inline void shade_two_color_gpu(unsigned char *r,
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

    *r = clamp_u8_gpu(base_r * brightness);
    *g = clamp_u8_gpu(base_g * brightness);
    *b = clamp_u8_gpu(base_b * brightness);
}

__global__ void render_kernel(unsigned char *img,
                              int width, int height,
                              double xmin, double xmax,
                              double ymin, double ymax,
                              int max_iter, double tol,
                              int method) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= width || j >= height) return;

    double dx = (xmax - xmin) / (double)(width - 1);
    double dy = (ymax - ymin) / (double)(height - 1);

    double x = xmin + i * dx;
    double y = ymax - j * dy;

    ResultGPU res;
    if (method == METHOD_SECANT) {
        res = iterate_secant_xy_gpu(x, y, max_iter, tol);
    } else {
        res = iterate_newton_xy_gpu(x, y, max_iter, tol);
    }

    unsigned char r, g, b;
    shade_two_color_gpu(&r, &g, &b, res.root_id, res.iter, max_iter, res.converged);

    size_t idx = (size_t)3 * ((size_t)j * (size_t)width + (size_t)i);
    img[idx + 0] = r;
    img[idx + 1] = g;
    img[idx + 2] = b;
}

int render_gpu(unsigned char *img, const FractalConfig *cfg) {
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count <= 0) {
        fprintf(stderr, "[gpu] No CUDA device found.\n");
        return 1;
    }

    size_t bytes = (size_t)3 * (size_t)cfg->width * (size_t)cfg->height;
    unsigned char *d_img = NULL;

    err = cudaMalloc((void **)&d_img, bytes);
    if (err != cudaSuccess) {
        fprintf(stderr, "[gpu] cudaMalloc failed: %s\n", cudaGetErrorString(err));
        return 1;
    }

    dim3 block(16, 16);
    dim3 grid((cfg->width + block.x - 1) / block.x,
              (cfg->height + block.y - 1) / block.y);

    render_kernel<<<grid, block>>>(
        d_img,
        cfg->width, cfg->height,
        cfg->xmin, cfg->xmax,
        cfg->ymin, cfg->ymax,
        cfg->max_iter, cfg->tol,
        (int)cfg->method
    );

    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "[gpu] kernel launch failed: %s\n", cudaGetErrorString(err));
        cudaFree(d_img);
        return 1;
    }

    err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        fprintf(stderr, "[gpu] kernel execution failed: %s\n", cudaGetErrorString(err));
        cudaFree(d_img);
        return 1;
    }

    err = cudaMemcpy(img, d_img, bytes, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "[gpu] cudaMemcpy failed: %s\n", cudaGetErrorString(err));
        cudaFree(d_img);
        return 1;
    }

    cudaFree(d_img);
    return 0;
}
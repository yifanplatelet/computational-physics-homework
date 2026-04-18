/*
 * newton_fractal.c
 *
 * 使用 Newton-Raphson 方法在复平面上绘制分形
 * 输出: PPM(P6) 图片，速度快，无需额外图形库
 *
 * 编译:
 *   gcc -O3 -std=c11 newton_fractal.c -lm -o newton_fractal
 *
 * 若支持 OpenMP 并想并行加速:
 *   gcc -O3 -std=c11 -fopenmp newton_fractal.c -lm -o newton_fractal
 *
 * 运行:
 *   ./newton_fractal
 *
 * 生成:
 *   newton_root.ppm    -> 按根着色
 *   newton_iter.ppm    -> 按迭代次数着色
 *
 * 你可以后续直接替换 equation 部分，快速改成其它方程。
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =========================================================
 * 基本数据结构
 * ========================================================= */

typedef struct {
    unsigned char r, g, b;
} RGB;

typedef struct {
    int width;
    int height;
    double xmin, xmax;
    double ymin, ymax;
    int max_iter;
    double tol;          // 收敛判据
    double deriv_eps;    // 导数接近0时停止
} RenderConfig;

typedef enum {
    COLOR_BY_ROOT = 0,
    COLOR_BY_ITER = 1
} ColorMode;

/* =========================================================
 * 方程接口：将方程和导数分离，方便替换
 * ========================================================= */

typedef double complex (*FuncComplex)(double complex z, void *params);

typedef struct {
    FuncComplex f;
    FuncComplex df;
    void *params;

    /* 若已知根，可填入用于分类着色 */
    double complex *roots;
    int root_count;
} Equation;

/* =========================================================
 * 示例 1：f(z) = z^n - 1
 * 这类方程非常适合做 Newton 分形，与你给的图高度一致
 * ========================================================= */

typedef struct {
    int n;
} ZnMinusOneParams;

static double complex f_zn_minus_one(double complex z, void *params) {
    ZnMinusOneParams *p = (ZnMinusOneParams *)params;
    return cpow(z, p->n) - 1.0;
}

static double complex df_zn_minus_one(double complex z, void *params) {
    ZnMinusOneParams *p = (ZnMinusOneParams *)params;
    if (p->n == 0) return 0.0;
    if (p->n == 1) return 1.0;
    return p->n * cpow(z, p->n - 1);
}

static double complex *build_roots_zn_minus_one(int n) {
    double complex *roots = (double complex *)malloc(sizeof(double complex) * n);
    if (!roots) return NULL;

    for (int k = 0; k < n; ++k) {
        double theta = 2.0 * M_PI * k / n;
        roots[k] = cos(theta) + I * sin(theta);
    }
    return roots;
}

/* =========================================================
 * 示例 2：任意复系数多项式
 *   p(z) = a[0] + a[1] z + a[2] z^2 + ... + a[n] z^n
 * 用 Horner 法评估
 * ========================================================= */

typedef struct {
    int degree;
    double complex *a;   // 长度 degree+1
    double complex *da;  // 导数系数，长度 degree
} PolyParams;

static double complex poly_eval(double complex z, const double complex *a, int degree) {
    double complex val = a[degree];
    for (int i = degree - 1; i >= 0; --i) {
        val = val * z + a[i];
    }
    return val;
}

static double complex f_poly(double complex z, void *params) {
    PolyParams *p = (PolyParams *)params;
    return poly_eval(z, p->a, p->degree);
}

static double complex df_poly(double complex z, void *params) {
    PolyParams *p = (PolyParams *)params;
    if (p->degree <= 0) return 0.0;
    return poly_eval(z, p->da, p->degree - 1);
}

static PolyParams *create_poly_params(const double complex *coeff, int degree) {
    PolyParams *p = (PolyParams *)malloc(sizeof(PolyParams));
    if (!p) return NULL;

    p->degree = degree;
    p->a = (double complex *)malloc(sizeof(double complex) * (degree + 1));
    p->da = (degree > 0) ? (double complex *)malloc(sizeof(double complex) * degree) : NULL;

    if (!p->a || (degree > 0 && !p->da)) {
        free(p->a);
        free(p->da);
        free(p);
        return NULL;
    }

    for (int i = 0; i <= degree; ++i) p->a[i] = coeff[i];
    for (int i = 1; i <= degree; ++i) p->da[i - 1] = i * coeff[i];

    return p;
}

static void free_poly_params(PolyParams *p) {
    if (!p) return;
    free(p->a);
    free(p->da);
    free(p);
}

/* =========================================================
 * 工具函数
 * ========================================================= */

static double clamp01(double x) {
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

static int nearest_root_index(double complex z, const double complex *roots, int n) {
    int best = -1;
    double bestd = 1e100;
    for (int i = 0; i < n; ++i) {
        double d = cabs(z - roots[i]);
        if (d < bestd) {
            bestd = d;
            best = i;
        }
    }
    return best;
}

/* HSV -> RGB，画出来更好看 */
static RGB hsv_to_rgb(double h, double s, double v) {
    RGB out = {0, 0, 0};

    h = fmod(h, 360.0);
    if (h < 0) h += 360.0;

    s = clamp01(s);
    v = clamp01(v);

    double c = v * s;
    double x = c * (1.0 - fabs(fmod(h / 60.0, 2.0) - 1.0));
    double m = v - c;

    double r1 = 0, g1 = 0, b1 = 0;

    if (h < 60)      { r1 = c; g1 = x; b1 = 0; }
    else if (h < 120){ r1 = x; g1 = c; b1 = 0; }
    else if (h < 180){ r1 = 0; g1 = c; b1 = x; }
    else if (h < 240){ r1 = 0; g1 = x; b1 = c; }
    else if (h < 300){ r1 = x; g1 = 0; b1 = c; }
    else             { r1 = c; g1 = 0; b1 = x; }

    out.r = (unsigned char)(255.0 * (r1 + m));
    out.g = (unsigned char)(255.0 * (g1 + m));
    out.b = (unsigned char)(255.0 * (b1 + m));
    return out;
}

static void write_ppm(const char *filename, const RGB *img, int w, int h) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "无法写入文件: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "P6\n%d %d\n255\n", w, h);
    fwrite(img, sizeof(RGB), (size_t)w * h, fp);
    fclose(fp);
}

/* =========================================================
 * Newton 迭代结果
 * ========================================================= */

typedef struct {
    int converged;
    int iter;
    int root_index;
    double complex z_final;
    double residual;
} NewtonResult;

static NewtonResult newton_iterate(double complex z0, const Equation *eq, const RenderConfig *cfg) {
    NewtonResult res;
    res.converged = 0;
    res.iter = 0;
    res.root_index = -1;
    res.z_final = z0;
    res.residual = 1e100;

    double complex z = z0;

    for (int k = 0; k < cfg->max_iter; ++k) {
        double complex fz = eq->f(z, eq->params);
        double complex dfz = eq->df(z, eq->params);
        double dfabs = cabs(dfz);

        if (dfabs < cfg->deriv_eps) {
            res.iter = k;
            res.z_final = z;
            res.residual = cabs(fz);
            return res;
        }

        double complex z_next = z - fz / dfz;

        if (cabs(z_next - z) < cfg->tol || cabs(fz) < cfg->tol) {
            z = z_next;
            res.converged = 1;
            res.iter = k + 1;
            res.z_final = z;
            res.residual = cabs(eq->f(z, eq->params));
            if (eq->roots && eq->root_count > 0) {
                res.root_index = nearest_root_index(z, eq->roots, eq->root_count);
            }
            return res;
        }

        z = z_next;
    }

    res.iter = cfg->max_iter;
    res.z_final = z;
    res.residual = cabs(eq->f(z, eq->params));
    if (eq->roots && eq->root_count > 0) {
        res.root_index = nearest_root_index(z, eq->roots, eq->root_count);
    }
    return res;
}

/* =========================================================
 * 着色方案
 * ========================================================= */

/*
 * 模式1：按收敛到哪个根着色
 * - 不同根给不同色相
 * - 收敛越快越亮
 * - 不收敛/导数太小 -> 白色或近黑色
 */
static RGB color_by_root(const NewtonResult *r, const Equation *eq, const RenderConfig *cfg) {
    if (!r->converged || r->root_index < 0 || eq->root_count <= 0) {
        RGB white = {255, 255, 255};
        return white;
    }

    double hue = 360.0 * r->root_index / eq->root_count;

    /* 快速收敛更亮，边界更暗，分形边界会更明显 */
    double t = (double)r->iter / cfg->max_iter;
    double v = 1.0 - 0.75 * pow(t, 0.65);   // 亮度
    double s = 0.80;

    return hsv_to_rgb(hue, s, v);
}

/*
 * 模式2：按迭代次数着色
 * - 颜色反映迭代次数
 * - 不同根仍可叠加轻微色相区分，增强层次
 */
static RGB color_by_iteration(const NewtonResult *r, const Equation *eq, const RenderConfig *cfg) {
    if (!r->converged) {
        RGB black = {0, 0, 0};
        return black;
    }

    double t = (double)r->iter / cfg->max_iter;

    /* 周期色带，能得到类似课件里的彩色层纹 */
    double hue = 720.0 * t;

    if (r->root_index >= 0 && eq->root_count > 0) {
        hue += 360.0 * r->root_index / eq->root_count * 0.25;
    }

    /* 中心和边界都尽量清晰 */
    double s = 0.85;
    double v = 0.95 - 0.25 * pow(t, 0.5);

    return hsv_to_rgb(hue, s, v);
}

/* =========================================================
 * 渲染
 * ========================================================= */

static void render_newton_fractal(
    RGB *img,
    const Equation *eq,
    const RenderConfig *cfg,
    ColorMode mode
) {
    int w = cfg->width;
    int h = cfg->height;

    double dx = (cfg->xmax - cfg->xmin) / (w - 1);
    double dy = (cfg->ymax - cfg->ymin) / (h - 1);

    /* 简单 2x2 超采样，可开关。为了速度这里默认关掉。 */
    const int antialias = 0;

    #pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < h; ++j) {
        for (int i = 0; i < w; ++i) {
            RGB color = {0, 0, 0};

            if (!antialias) {
                double x = cfg->xmin + i * dx;
                double y = cfg->ymax - j * dy;  // 图像y向下，复平面y向上
                double complex z0 = x + I * y;

                NewtonResult r = newton_iterate(z0, eq, cfg);

                if (mode == COLOR_BY_ROOT)
                    color = color_by_root(&r, eq, cfg);
                else
                    color = color_by_iteration(&r, eq, cfg);

            } else {
                /* 2x2 抗锯齿 */
                int count = 0;
                double rs = 0, gs = 0, bs = 0;
                for (int sy = 0; sy < 2; ++sy) {
                    for (int sx = 0; sx < 2; ++sx) {
                        double x = cfg->xmin + (i + (sx + 0.5) / 2.0) * dx;
                        double y = cfg->ymax - (j + (sy + 0.5) / 2.0) * dy;
                        double complex z0 = x + I * y;
                        NewtonResult r = newton_iterate(z0, eq, cfg);
                        RGB c = (mode == COLOR_BY_ROOT)
                                    ? color_by_root(&r, eq, cfg)
                                    : color_by_iteration(&r, eq, cfg);
                        rs += c.r;
                        gs += c.g;
                        bs += c.b;
                        count++;
                    }
                }
                color.r = (unsigned char)(rs / count);
                color.g = (unsigned char)(gs / count);
                color.b = (unsigned char)(bs / count);
            }

            img[j * w + i] = color;
        }
    }
}

/* =========================================================
 * 主程序
 * ========================================================= */

int main(void) {
    /* -----------------------------
     * 1) 配置渲染参数
     * ----------------------------- */
    RenderConfig cfg;
    cfg.width = 1600;
    cfg.height = 1600;

    /*
     * 对 z^n - 1，一般取 [-2,2] x [-2,2] 很经典
     * 也可以后面改成局部放大
     */
    cfg.xmin = -2.0;
    cfg.xmax =  2.0;
    cfg.ymin = -2.0;
    cfg.ymax =  2.0;

    cfg.max_iter = 40;
    cfg.tol = 1e-8;
    cfg.deriv_eps = 1e-14;

    //修改预想的方程：
    //1.若为z^n-1=0型：
    //ZnMinusOneParams p1;
    //p1.n = 6;
    //2.若为多项式型：
    //double complex coeff[]={-1.0,0.0,0.0,1.0};
    //PolyParams *pp=create_poly_params(coeff,3);

    
    Equation eq;
    //第一种情况
    //eq.f = f_zn_minus_one;
    //eq.df = df_zn_minus_one;
    //eq.params = &p1;
    //eq.roots = build_roots_zn_minus_one(pp.n);
    //eq.root_count = p1.n;

    //第二种情况
    double complex coeff[] = {-1.0, 0.0, 0.0, 1.0};
    PolyParams *pp = create_poly_params(coeff, 3);
    eq.f = f_poly;
    eq.df = df_poly;
    eq.params = pp;
    eq.roots = NULL;      // 若不知道根，可先不填
    eq.root_count = 0;

    /* -----------------------------
     * 3) 分配图像内存
     * ----------------------------- */
    RGB *img = (RGB *)malloc(sizeof(RGB) * cfg.width * cfg.height);
    if (!img) {
        fprintf(stderr, "图像内存分配失败\n");
        free(eq.roots);
        return EXIT_FAILURE;
    }

    /* -----------------------------
     * 4) 输出 1：按根着色
     * ----------------------------- */
    render_newton_fractal(img, &eq, &cfg, COLOR_BY_ROOT);
    write_ppm("newton_root.ppm", img, cfg.width, cfg.height);
    printf("已生成: newton_root.ppm\n");

    /* -----------------------------
     * 5) 输出 2：按迭代次数着色
     * ----------------------------- */
    render_newton_fractal(img, &eq, &cfg, COLOR_BY_ITER);
    write_ppm("newton_iter.ppm", img, cfg.width, cfg.height);
    printf("已生成: newton_iter.ppm\n");

    /* -----------------------------
     * 6) 清理
     * ----------------------------- */
    free(img);
    free(eq.roots);

    /*
     * 若你用了 create_poly_params(...)，记得:
     * free_poly_params(pp);
     */

    return 0;
}
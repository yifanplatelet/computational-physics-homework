#ifndef CPU_BACKEND_H
#define CPU_BACKEND_H

typedef enum {
    METHOD_NEWTON = 0,
    METHOD_SECANT = 1
} FractalMethod;

typedef enum {
    SYMMETRY_OFF = 0,
    SYMMETRY_AUTO = 1,
    SYMMETRY_ON = 2
} SymmetryMode;

typedef struct {
    int width;
    int height;
    int max_iter;
    double tol;
    double xmin, xmax, ymin, ymax;
    FractalMethod method;
    int threads;   /* 0 means use maximum available */
    SymmetryMode symmetry_mode; /* off/auto/on */
} FractalConfig;

void render_cpu(unsigned char *img, const FractalConfig *cfg);

#endif

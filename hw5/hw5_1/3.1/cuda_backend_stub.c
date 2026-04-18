#include "cuda_backend.h"
#include <stdio.h>

int render_gpu(unsigned char *img, const FractalConfig *cfg) {
    (void)img;
    (void)cfg;
    fprintf(stderr, "[gpu] CUDA backend not compiled in.\n");
    return 1;
}
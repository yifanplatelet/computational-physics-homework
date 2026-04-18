#ifndef CUDA_BACKEND_H
#define CUDA_BACKEND_H

#include "cpu_backend.h"

/* return 0 on success, non-zero on failure */
int render_gpu(unsigned char *img, const FractalConfig *cfg);

#endif
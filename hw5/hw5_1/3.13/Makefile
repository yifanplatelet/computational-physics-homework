CC      = gcc
NVCC    = nvcc

# 通用基础选项
BASE_CFLAGS  = -std=c11 -fopenmp
BASE_LDFLAGS = -lm -fopenmp

# Newton 更适合：
# - AVX/本机指令集
# - 更激进的展开
# - LTO
NEWTON_CFLAGS  = -Ofast -march=native -flto -funroll-loops $(BASE_CFLAGS)
NEWTON_LDFLAGS = $(BASE_LDFLAGS) -flto

# Secant 更适合：
# - 本机指令集
# - LTO
# - Ofast
# - 不强行再堆额外展开，避免代码体积过大带来的副作用
SECANT_CFLAGS  = -Ofast -march=native -flto $(BASE_CFLAGS)
SECANT_LDFLAGS = $(BASE_LDFLAGS) -flto

# 如果以后要编 CUDA 版，可继续保留
NVCCFLAGS = -O3 -Xcompiler -fopenmp

# 默认目标：给你一个通用版
all: fractal_cpu

# -------------------------
# 通用 CPU 版
# -------------------------
fractal_cpu: main.o cpu_backend.o cuda_backend_stub.o
	$(CC) $(SECANT_CFLAGS) -o $@ main.o cpu_backend.o cuda_backend_stub.o $(SECANT_LDFLAGS)

main.o: main.c cpu_backend.h cuda_backend.h
	$(CC) $(SECANT_CFLAGS) -c main.c -o main.o

cpu_backend.o: cpu_backend.c cpu_backend.h
	$(CC) $(SECANT_CFLAGS) -c cpu_backend.c -o cpu_backend.o

cuda_backend_stub.o: cuda_backend_stub.c cuda_backend.h cpu_backend.h
	$(CC) $(SECANT_CFLAGS) -c cuda_backend_stub.c -o cuda_backend_stub.o

# -------------------------
# Newton 专用编译
# -------------------------
fractal_newton: main_newton.o cpu_backend_newton.o cuda_backend_stub_newton.o
	$(CC) $(NEWTON_CFLAGS) -o $@ main_newton.o cpu_backend_newton.o cuda_backend_stub_newton.o $(NEWTON_LDFLAGS)

main_newton.o: main.c cpu_backend.h cuda_backend.h
	$(CC) $(NEWTON_CFLAGS) -c main.c -o main_newton.o

cpu_backend_newton.o: cpu_backend.c cpu_backend.h
	$(CC) $(NEWTON_CFLAGS) -c cpu_backend.c -o cpu_backend_newton.o

cuda_backend_stub_newton.o: cuda_backend_stub.c cuda_backend.h cpu_backend.h
	$(CC) $(NEWTON_CFLAGS) -c cuda_backend_stub.c -o cuda_backend_stub_newton.o

# -------------------------
# Secant 专用编译
# -------------------------
fractal_secant: main_secant.o cpu_backend_secant.o cuda_backend_stub_secant.o
	$(CC) $(SECANT_CFLAGS) -o $@ main_secant.o cpu_backend_secant.o cuda_backend_stub_secant.o $(SECANT_LDFLAGS)

main_secant.o: main.c cpu_backend.h cuda_backend.h
	$(CC) $(SECANT_CFLAGS) -c main.c -o main_secant.o

cpu_backend_secant.o: cpu_backend.c cpu_backend.h
	$(CC) $(SECANT_CFLAGS) -c cpu_backend.c -o cpu_backend_secant.o

cuda_backend_stub_secant.o: cuda_backend_stub.c cuda_backend.h cpu_backend.h
	$(CC) $(SECANT_CFLAGS) -c cuda_backend_stub.c -o cuda_backend_stub_secant.o

# -------------------------
# CUDA 版（保留接口）
# 如果以后你的老师机器装了 CUDA，再用
# -------------------------
cuda_backend_cuda.o: cuda_backend.cu cuda_backend.h cpu_backend.h
	$(NVCC) $(NVCCFLAGS) -c cuda_backend.cu -o cuda_backend_cuda.o

fractal_gpu: main.o cpu_backend.o cuda_backend_cuda.o
	$(NVCC) -o $@ main.o cpu_backend.o cuda_backend_cuda.o -Xcompiler -fopenmp -lm

# -------------------------
# 清理
# -------------------------
clean:
	rm -f *.o fractal_cpu fractal_newton fractal_secant fractal_gpu
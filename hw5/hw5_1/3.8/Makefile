CC      = gcc
CFLAGS  = -O3 -std=c11 -march=native -fopenmp
LDFLAGS = -lm -fopenmp

NVCC        = nvcc
NVCCFLAGS   = -O3 -Xcompiler -fopenmp

CPU_OBJS    = main.o cpu_backend.o cuda_backend_stub.o
CUDA_OBJS   = main.o cpu_backend.o cuda_backend_cuda.o

all: fractal_cpu

fractal_cpu: $(CPU_OBJS)
	$(CC) $(CFLAGS) -o $@ $(CPU_OBJS) $(LDFLAGS)

main.o: main.c cpu_backend.h cuda_backend.h
	$(CC) $(CFLAGS) -c main.c

cpu_backend.o: cpu_backend.c cpu_backend.h
	$(CC) $(CFLAGS) -c cpu_backend.c

cuda_backend_stub.o: cuda_backend_stub.c cuda_backend.h cpu_backend.h
	$(CC) $(CFLAGS) -c cuda_backend_stub.c

cuda_backend_cuda.o: cuda_backend.cu cuda_backend.h cpu_backend.h
	$(NVCC) $(NVCCFLAGS) -c cuda_backend.cu -o cuda_backend_cuda.o

fractal_gpu: $(CUDA_OBJS)
	$(NVCC) -o $@ $(CUDA_OBJS) -Xcompiler -fopenmp -lm

clean:
	rm -f *.o fractal_cpu fractal_gpu
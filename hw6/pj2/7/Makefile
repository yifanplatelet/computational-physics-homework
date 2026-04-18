CC=gcc
CFLAGS=-O3 -march=native -flto -std=c11 -fopenmp
LDLIBS=-lmpfr -lgmp

all: pi_bs_openmp_rewrite

pi_bs_openmp_rewrite: pi_bs_openmp_rewrite.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f pi_bs_openmp_rewrite

CC ?= gcc
CFLAGS ?= -O3 -march=native -std=c11 -fopenmp
LDFLAGS ?= -lmpfr -lgmp -lm

all: pi_bs_openmp_baseline pi_bs_openmp_factorized

pi_bs_openmp_baseline: src/pi_bs_openmp_baseline.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

pi_bs_openmp_factorized: src/pi_bs_openmp_factorized.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f pi_bs_openmp_baseline pi_bs_openmp_factorized

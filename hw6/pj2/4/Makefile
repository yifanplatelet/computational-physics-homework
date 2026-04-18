CC = gcc
CFLAGS = -O3 -march=native -std=c11
LIBS = -lmpfr -lgmp

all: pi_bs_single_opt pi_bs_openmp_opt

pi_bs_single_opt: bs_single/pi_bs_single_opt.c
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

pi_bs_openmp_opt: bs_openmp/pi_bs_openmp_opt.c
	$(CC) $(CFLAGS) -fopenmp $< -o $@ $(LIBS)

clean:
	rm -f pi_bs_single_opt pi_bs_openmp_opt

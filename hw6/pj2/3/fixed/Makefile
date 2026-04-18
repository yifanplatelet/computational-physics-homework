CC = gcc
CFLAGS = -O3 -march=native -std=c11
LIBS = -lmpfr -lgmp

all: pi_bs_single pi_bs_openmp

pi_bs_single: bs_single/pi_bs_single.c
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

pi_bs_openmp: bs_openmp/pi_bs_openmp.c
	$(CC) $(CFLAGS) -fopenmp $< -o $@ $(LIBS)

clean:
	rm -f pi_bs_single pi_bs_openmp

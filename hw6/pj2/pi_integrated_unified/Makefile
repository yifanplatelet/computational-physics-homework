CC=gcc
CFLAGS=-O3 -march=native -std=c11 -fopenmp
LDFLAGS=-lmpfr -lgmp

all: pi_bs_unified_exact

pi_bs_unified_exact: pi_bs_unified_exact.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f pi_bs_unified_exact

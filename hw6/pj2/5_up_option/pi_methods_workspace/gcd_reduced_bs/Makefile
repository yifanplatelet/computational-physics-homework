CC=gcc
CFLAGS=-O3 -march=native -std=c11
LIBS=-lmpfr -lgmp

all: pi_bs_gcd_reduced

pi_bs_gcd_reduced: pi_bs_gcd_reduced.c
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

clean:
	rm -f pi_bs_gcd_reduced

CC=gcc
CFLAGS=-O3 -march=native -std=c11
LIBS=-lmpfr -lgmp

all: pi_bs_skewed_backward

pi_bs_skewed_backward: pi_bs_skewed_backward.c
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

clean:
	rm -f pi_bs_skewed_backward

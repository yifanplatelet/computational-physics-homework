CC=gcc
CFLAGS=-O3 -march=native -std=c11
LIBS=-lmpfr -lgmp

all: pi_bs_top_truncated

pi_bs_top_truncated: pi_bs_top_truncated.c
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

clean:
	rm -f pi_bs_top_truncated

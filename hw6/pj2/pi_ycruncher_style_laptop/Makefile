CC=gcc
CFLAGS=-O3 -march=native -fopenmp -std=c11 -Wall -Wextra
LDFLAGS=-lmpfr -lgmp

all: pi_ycruncher_style

pi_ycruncher_style: pi_ycruncher_style.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f pi_ycruncher_style

CC = gcc
CFLAGS = -O3 -march=native -std=c11 -fopenmp
LDLIBS = -lmpfr -lgmp

all: pi_bs_baseline_from_main pi_bs_yc_inspired

pi_bs_baseline_from_main: pi_bs_baseline_from_main.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

pi_bs_yc_inspired: pi_bs_yc_inspired.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f pi_bs_baseline_from_main pi_bs_yc_inspired

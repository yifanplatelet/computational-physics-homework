CC ?= gcc
CFLAGS ?= -Ofast -march=native -flto -fomit-frame-pointer -funroll-loops -fopenmp -std=c11 -DNDEBUG
LDLIBS ?= -lgmp

all: pi_bs_extreme_mpf

pi_bs_extreme_mpf: pi_bs_extreme_mpf.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

clean:
	rm -f pi_bs_extreme_mpf

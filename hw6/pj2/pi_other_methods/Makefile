CC = gcc
CFLAGS = -O3 -march=native -std=c11 -Wall -Wextra
OMPFLAGS = -fopenmp
LIBS = -lmpfr -lgmp

all: \
	pi_chudnovsky_openmp_baseline \
	pi_ramanujan_openmp_bs \
	pi_gauss_legendre \
	pi_borwein_quartic \
	pi_machin_mpfr

pi_chudnovsky_openmp_baseline:
	$(CC) $(CFLAGS) $(OMPFLAGS) chudnovsky_openmp_baseline/pi_chudnovsky_openmp_baseline.c -o $@ $(LIBS)

pi_ramanujan_openmp_bs:
	$(CC) $(CFLAGS) $(OMPFLAGS) ramanujan_openmp_bs/pi_ramanujan_openmp_bs.c -o $@ $(LIBS)

pi_gauss_legendre:
	$(CC) $(CFLAGS) gauss_legendre/pi_gauss_legendre.c -o $@ $(LIBS)

pi_borwein_quartic:
	$(CC) $(CFLAGS) borwein_quartic/pi_borwein_quartic.c -o $@ $(LIBS)

pi_machin_mpfr:
	$(CC) $(CFLAGS) machin_mpfr/pi_machin_mpfr.c -o $@ $(LIBS)

clean:
	rm -f pi_chudnovsky_openmp_baseline pi_ramanujan_openmp_bs pi_gauss_legendre pi_borwein_quartic pi_machin_mpfr

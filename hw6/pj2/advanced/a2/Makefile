CC ?= gcc
CFLAGS ?= -O3 -march=native -fopenmp -std=c11
LDLIBS ?= -lgmp -lmpfr

all: pi_bs_baseline_from_user pi_bs_optimized_yc_learned pi_bs_fixedpoint_yc_adaptive pi_bs_gmp_kernel_tuned

pi_bs_baseline_from_user: pi_bs_baseline_from_user.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

pi_bs_optimized_yc_learned: pi_bs_optimized_yc_learned.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

pi_bs_fixedpoint_yc_adaptive: pi_bs_fixedpoint_yc_adaptive.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

pi_bs_gmp_kernel_tuned: pi_bs_gmp_kernel_tuned.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

clean:
	rm -f pi_bs_baseline_from_user pi_bs_optimized_yc_learned pi_bs_fixedpoint_yc_adaptive pi_bs_gmp_kernel_tuned

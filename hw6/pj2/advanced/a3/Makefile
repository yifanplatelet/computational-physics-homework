CC ?= gcc
CFLAGS ?= -Ofast -march=native -flto -fomit-frame-pointer -funroll-loops -fopenmp -std=c11 -DNDEBUG
LDLIBS ?= -lmpfr -lgmp

all: pi_bs_qt_tuned pi_bs_baseline_plus pi_bs_extreme_trimmed

pi_bs_qt_tuned: pi_bs_qt_tuned.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

pi_bs_baseline_plus: pi_bs_baseline_plus.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

pi_bs_extreme_trimmed: pi_bs_extreme_trimmed.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

clean:
	rm -f pi_bs_qt_tuned pi_bs_baseline_plus pi_bs_extreme_trimmed

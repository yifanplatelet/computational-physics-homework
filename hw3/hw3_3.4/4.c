#include <stdio.h>

int main() {
    float feps = 1.0f;
    double deps = 1.0;
    long double leps = 1.0L;

    while ((1.0f + feps / 2.0f) != 1.0f) {
        feps /= 2.0f;
    }

    while ((1.0 + deps / 2.0) != 1.0) {
        deps /= 2.0;
    }

    while ((1.0L + leps / 2.0L) != 1.0L) {
        leps /= 2.0L;
    }

    printf("Machine epsilon for float      = %.20e\n", feps);
    printf("Machine epsilon for double     = %.20e\n", deps);
    printf("Machine epsilon for long double= %.20Le\n", leps);

    return 0;
}
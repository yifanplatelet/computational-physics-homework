#include <stdio.h>
#include <math.h>

double f(double x) {
    return x * exp(x);
}

int main() {
    double x = 2.0;
    double exact = 3.0 * exp(2.0);

    printf("Exact f'(2) = %.15f\n\n", exact);

    printf("%-12s %-15s %-15s %-15s %-15s %-15s %-15s\n",
           "h",
           "forward",
           "backward",
           "central",
           "forward3",
           "backward3",
           "central5");

    for (int k = 1; k <= 12; k++) {
        double h = pow(10.0, -k);

        double forward =
            (f(x + h) - f(x)) / h;

        double backward =
            (f(x) - f(x - h)) / h;

        double central =
            (f(x + h) - f(x - h)) / (2.0 * h);

        double forward3 =
            (-3.0 * f(x) + 4.0 * f(x + h) - f(x + 2.0 * h)) / (2.0 * h);

        double backward3 =
            (3.0 * f(x) - 4.0 * f(x - h) + f(x - 2.0 * h)) / (2.0 * h);

        double central5 =
            (f(x - 2.0 * h)
             - 8.0 * f(x - h)
             + 8.0 * f(x + h)
             - f(x + 2.0 * h)) / (12.0 * h);

        printf("%-12.1e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n",
               h,
               fabs(forward - exact),
               fabs(backward - exact),
               fabs(central - exact),
               fabs(forward3 - exact),
               fabs(backward3 - exact),
               fabs(central5 - exact));
    }

    return 0;
}
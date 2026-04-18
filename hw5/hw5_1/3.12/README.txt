Changes in this package:
- Secant now uses RGB root coloring (root 0 red, root 1 green, root 2 blue).
- Secant uses a purely real perturbation z1 = z0 + eps, with eps = 4*dx.
- Degenerate/non-convergent secant points fall back to a short Newton probe for root classification only,
  which removes the large black center disk while keeping secant as the main iteration.
- Main defaults were reduced to 4000x4000, max_iter=60 so the program is runnable by default.

Recommended secant command:
./fractal_cpu --backend=cpu --method=secant --width=4000 --height=4000 --max_iter=100 --tol=1e-6 --output=secant_rgb.ppm

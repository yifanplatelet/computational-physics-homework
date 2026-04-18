Recommended commands
====================

1) Build CPU version
   make clean
   make fractal_cpu

2) Newton (keep current optimized two-color output)
   ./fractal_cpu --backend=cpu --method=newton --width=4000 --height=4000 --max_iter=50 --tol=1e-6 --output=newton.ppm

3) Secant with teacher-like RGB coloring
   ./fractal_cpu --backend=cpu --method=secant --width=4000 --height=4000 --max_iter=100 --tol=1e-6 --output=secant_rgb.ppm

Notes
-----
- The secant renderer now uses a purely real second initial guess:
    z1 = z0 + eps, eps = 4 * dx
  which preserves the real-axis symmetry and gives a cleaner RGB basin image.
- Secant is intentionally colored in three RGB basins blended toward white near
  slow-converging boundaries, matching the teacher-style figure more closely.

# CPU DIC speedups

The `speedup` branch improves the CPU implementation of `ncorr_alg_rgdic` without reducing DIC resolution or changing convergence settings.

## Changes

- Parallelized interpolation-table and reference-gradient precomputation with OpenMP.
- Removed redundant zero-initialization of buffers that are immediately overwritten. This avoids roughly 18–19 GiB of unnecessary memory writes for an 8K image.
- Replaced dynamically allocated seed and deformation vectors with fixed-size arrays.
- Rewrote the biquintic interpolation hot loop in Horner form to reduce arithmetic and memory traffic.
- Replaced repeated `pow(x,2)` operations with direct multiplication and reused the Hessian scale factor.

All interpolation coefficients and solver calculations remain double precision. Tests produced identical valid-point masks, with displacement differences below approximately `3e-14` pixels.

The included `ncorr_alg_rgdic.mexw64` was compiled with OpenMP support. When recompiling Ncorr, OpenMP must remain enabled to benefit from parallel precomputation.

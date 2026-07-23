# Neural Spline Flow (NSF) conditional density estimator

An autoregressive flow whose per-dimension transform is a monotonic
rational-quadratic spline (Durkan et al., 2019) instead of MAF's affine
shift-and-scale. Splines are far more expressive per layer, which helps
on sharply non-Gaussian posteriors (SLCP, two moons). We reuse the MADE
masking machinery from `R/flows.R`; each MADE outputs `3K - 1` spline
parameters per dimension (K bin widths, K bin heights, K - 1 interior
derivatives). The spline acts on `[-B, B]` and is the identity outside
(linear tails), so the standard-normal base distribution is unaffected
in the tails. Note: NSF implementations elsewhere often use coupling
layers; ours is autoregressive — same density family, different
conditioning structure.

## Details

Select with `npe(..., density_estimator = "nsf")`.

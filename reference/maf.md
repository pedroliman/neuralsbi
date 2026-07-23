# Masked Autoregressive Flow (MAF) conditional density estimator

A normalizing flow maps parameters \\\theta\\ to a standard-normal base
variable through a stack of invertible transforms, giving exact
densities by the change of variables. The MAF (Papamakarios et al.,
2017) uses masked autoregressive networks (MADE, Germain et al., 2015):
each transform is \$\$u_d = (\theta_d - \mu_d(\theta\_{\<d}, x))
\exp(-\alpha_d(\theta\_{\<d}, x)),\$\$ where the masks guarantee that
\\\mu_d, \alpha_d\\ depend only on earlier dimensions of \\\theta\\ (and
freely on the conditioning data `x`). Density evaluation is a single
forward pass; sampling inverts the transform one dimension at a time.
Between transforms the parameter order is reversed so every dimension
gets conditioned on every other across the stack.

## Details

This is `sbi`'s default flow family, and the default estimator in
`neuralsbi` too. It handles non-Gaussian posteriors that the MDN
struggles with. It is selected by default, or explicitly with
`npe(..., density_estimator = "maf")`.

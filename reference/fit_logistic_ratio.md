# Fit the closed-form logistic ratio estimator on standardized (theta, x)

The atomic objective at two atoms has a closed form. For one simulation
\\(\theta_i, x_i)\\ and one contrast \\\theta_j\\, the two-way softmax
reduces to \\\log \sigma\\\left(f(\theta_i, x_i) - f(\theta_j,
x_i)\right)\\, and with \\f = w'\varphi\\ linear in the features that is
a logistic regression on the *difference* \\\varphi(\theta_i, x_i) -
\varphi(\theta_j, x_i)\\ with every label positive. So the whole
estimator is one ridge-penalized IRLS on an `n * (num_atoms - 1)` by
`ncol(Phi)` design, no optimizer and no torch.

## Usage

``` r
fit_logistic_ratio(theta, x, num_atoms = 10L, ridge = 1e-06, verbose = FALSE)
```

## Details

Working in differences is what makes the fit *exact* for a
linear-Gaussian simulator rather than merely close. The differences
cancel every term that depends on `x` alone, including the evidence
\\\log p(x)\\, which is the one part of the log ratio a quadratic basis
cannot represent. What is left is the parameter dependence \\\log p(x
\mid \theta)\\, which for a linear-Gaussian model lies exactly in the
span of
[`nre_features()`](https://neuralsbi.pedrodelima.com/reference/nre_features.md).
The price is the one the atomic objective always pays: the level of the
fitted ratio at a given `x` is arbitrary.

Contrasts come from cyclic shifts of the parameter rows rather than a
random draw. The rows are independent prior draws in random order
already, so a shift is as good a scramble, and it makes the fit
deterministic: the same simulations give the same estimator whether or
not a seed was set.

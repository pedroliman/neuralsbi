# Ridge-penalized logistic regression by iteratively reweighted least squares

[`stats::glm.fit()`](https://rdrr.io/r/stats/glm.html) would do this,
but it warns and wanders off when the two classes are separable, which a
sharp ratio makes easy to hit. The ridge term keeps the normal equations
solvable and the coefficients finite, the same role it plays in
`fit_linear_gaussian()`, and it is measured against each column's own
scale for the same reason: under `standardize = FALSE` the quadratic
features carry the fourth power of the data's units, so an absolute 1e-6
is either nothing at all or the only thing left. On a simulator whose
output has sd 5e-4 the absolute version shrank the fit to noise; the
relative one leaves it alone.

## Usage

``` r
irls_logistic(X, y, ridge = 1e-06, max_iter = 100L, tol = 1e-08)
```

# Standardization (z-scoring) helpers

Neural density estimators train far more reliably when inputs and
targets are standardized to roughly zero mean and unit variance.
`neuralsbi` learns these transforms from the training simulations,
applies them internally, and inverts them when returning posterior draws
/ densities.

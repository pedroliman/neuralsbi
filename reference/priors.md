# Priors for neural simulation-based inference

A prior in `neuralsbi` is a lightweight object (class `nsbi_prior`) that
knows how to (a) draw samples and (b) evaluate its log-density. Bounded
priors also carry `lower`/`upper` support limits, which are used to
reject out-of-support posterior samples ("leakage" correction).

# Two-moons calibration study (verification-roadmap Level 2, milestone M2).
#
# Two moons has no closed-form posterior, so we cannot score it against an exact
# answer. Instead we check that the fitted posterior is *calibrated*: SBC ranks
# should be uniform, expected coverage should track the diagonal, and TARP --
# a joint test -- should agree. TARP matters here because the two-moons
# posterior is a thin crescent: a fit can have calibrated per-parameter
# marginals while getting the joint (the correlation along the crescent) wrong,
# which only the joint test sees.
#
# Run with libtorch available:
#     R CMD INSTALL --no-docs . && Rscript inst/benchmarks/two_moons_calibration.R
# Figures are written to docs/figures/.

library(neuralsbi)
stopifnot(requireNamespace("torch", quietly = TRUE), torch::torch_is_installed())

set.seed(1)
out_dir <- "docs/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

task <- task_two_moons()

# A neural spline flow captures the crescent most cleanly (see the
# density-estimators vignette); the MDN is a cheaper alternative.
fit <- npe(task$prior, task$simulator, n_simulations = 2500,
           density_estimator = "nsf", max_epochs = 150, seed = 1)

# --- Simulation-based calibration + expected coverage ------------------------
res <- sbc(fit, task$simulator, n_sbc = 150, n_posterior_samples = 300,
           seed = 2)
cat("SBC per-parameter uniformity p-values:",
    paste(sprintf("%.3f", res$uniformity_pvalue), collapse = "  "), "\n")

png(file.path(out_dir, "two_moons_sbc.png"), width = 900, height = 400, res = 110)
op <- par(mfrow = c(1, 2))
plot_sbc(res, param = 1)
plot_sbc(res, param = 2)
par(op)
dev.off()

png(file.path(out_dir, "two_moons_coverage.png"), width = 500, height = 500, res = 110)
plot_coverage(res)
dev.off()

# --- TARP joint coverage -----------------------------------------------------
tr <- tarp(fit, task$simulator, n_tarp = 150, n_posterior_samples = 300,
           seed = 3)
png(file.path(out_dir, "two_moons_tarp.png"), width = 500, height = 500, res = 110)
plot_tarp(tr)
dev.off()

cat("Wrote two_moons_sbc.png, two_moons_coverage.png, two_moons_tarp.png to",
    out_dir, "\n")

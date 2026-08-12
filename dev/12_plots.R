# 12_plots.R -----------------------------------------------------------------
#
# Every plotting function the package exports, what each one is for, and how to
# get at the numbers behind it. All five return their data invisibly and draw
# as a side effect, so they work in a script and in a pipeline.
#
#   pairplot()                   posterior draws: marginals and 2-D regions
#   plot_sbc()                   SBC rank histogram
#   plot_coverage()              nominal vs empirical credible-interval coverage
#   plot_tarp()                  TARP expected-coverage curve
#   plot_posterior_predictive()  predictive draws against the observation
#
# Task source
#   Lueckmann, J.-M., Boelts, J., Greenberg, D., Goncalves, P. and Macke, J.
#   "Benchmarking Simulation-Based Inference", AISTATS 2021. The sbibm suite.
#   https://github.com/sbi-benchmark/sbibm
#   task_sir() is sbibm's sir task: an SIR model with contact rate beta and
#   recovery rate gamma under log-normal priors, beta ~ LogNormal(log 0.4, 0.5)
#   and gamma ~ LogNormal(log 0.125, 0.2), solved by Euler steps for a
#   population of one million over 160 days. The data are binomial subsamples
#   of the infected fraction at 10 evenly spaced times. Two parameters and ten
#   outcomes, which is a convenient size for every plot here.
#
# Runtime: about 5 minutes on a laptop CPU.

library(neuralsbi)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("This script needs ggplot2. install.packages('ggplot2')")
}
has_pairplot <- requireNamespace("GGally", quietly = TRUE) &&
  requireNamespace("ggdensity", quietly = TRUE)
has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()

outdir <- file.path(tempdir(), "neuralsbi-plots")
dir.create(outdir, showWarnings = FALSE)
save_plot <- function(p, name) {
  ggplot2::ggsave(file.path(outdir, name), p, width = 6, height = 4.5,
                  dpi = 110)
  invisible(p)
}

# ---------------------------------------------------------------------------
# Fit something to plot
# ---------------------------------------------------------------------------

task <- task_sir()
print(task)

# task_sir()'s prior is a prior_custom() without parameter names, so name the
# columns here and they will label every axis from now on.
prior <- prior_custom(
  sample_fn = task$prior$sample,
  log_prob_fn = task$prior$log_prob,
  dim = 2L, lower = c(0, 0),
  param_names = c("beta", "gamma")
)
simulator <- function(theta) {
  stats::setNames(task$simulator(theta), paste0("t", 1:10))
}

fit <- npe(prior, simulator, n_simulations = 4000,
           density_estimator = if (has_torch) "maf" else "linear_gaussian",
           max_epochs = 150L, patience = 10L, seed = 1)

theta_true <- c(beta = 0.4, gamma = 0.125)
set.seed(9)
x_obs <- simulator(theta_true)
post <- posterior(fit, x_obs = x_obs)
draws <- sample(post, 4000)

# ---------------------------------------------------------------------------
# print() and summary() on posterior draws
# ---------------------------------------------------------------------------
#
# Draws come back as an nsbi_samples matrix. It prints as a short report, not
# as 4000 rows of numbers, and it carries the rejection-sampling acceptance
# rate that produced it.

print(draws)
print(summary(draws))
print(summary(draws, probs = c(0.05, 0.5, 0.95)))

# ---------------------------------------------------------------------------
# pairplot()
# ---------------------------------------------------------------------------
#
# Marginal densities on the diagonal, highest-density regions in the lower
# triangle, shaded by probability level (99, 95, 80, 50%). Optionally a marker
# for a reference value.

if (has_pairplot) {
  p <- pairplot(draws, truth = theta_true)
  save_plot(p, "pairplot.png")

  # Labels default to the column names, which came from the prior. Anything
  # that parses as R syntax renders as its plotmath symbol, so "beta" draws as
  # the Greek letter and "R[0]" draws as a subscript.
  pairplot(draws, truth = theta_true, labels = c("beta", "gamma"))

  # limits fixes the axes across every panel, which matters when comparing two
  # fits side by side: without it each panel picks its own range.
  pairplot(draws, truth = theta_true,
           limits = list(c(0.2, 0.8), c(0.05, 0.25)),
           col = "darkorange", alpha = 0.3)
} else {
  cat("GGally and ggdensity are not installed: skipping pairplot().\n")
  cat("  install.packages(c('GGally', 'ggdensity'))\n")
}

# ---------------------------------------------------------------------------
# plot_posterior_predictive()
# ---------------------------------------------------------------------------
#
# One histogram per data dimension, with the observation as a vertical line.
# The invisible return is the fraction of predictive draws below the
# observation in each dimension, which is the tail probability the plot is
# showing. Values near 0 or 1 mark the dimensions the model does not reproduce,
# and are much easier to scan than ten panels.

pred <- posterior_predictive(post, simulator, n = 500)
q <- plot_posterior_predictive(pred, x_obs)
print(round(q, 2))
save_plot(ggplot2::last_plot(), "posterior_predictive.png")

# Expect several of these to come back at 0.00 or 1.00. The posterior here is
# very tight (ten binomial observations of 1000 draws each are informative), so
# the predictive band is narrow, and a fit that is off by a fraction of a
# posterior sd in gamma misses the observation at the ends of the time grid.
# The SBC below says this fit is calibrated ON AVERAGE, which is not the same
# as being right at this particular observation. Both statements are true and
# both are worth reporting.

# labels default to the simulator's output names.
plot_posterior_predictive(pred, x_obs,
                          labels = paste0("day", round(seq(1, 160, length.out = 10))),
                          bins = 20L)

# ---------------------------------------------------------------------------
# plot_sbc()
# ---------------------------------------------------------------------------

sbc_res <- sbc(fit, simulator, n_sbc = 200L, n_posterior_samples = 400L,
               seed = 1)
print(sbc_res)

# param takes an index or a parameter name.
plot_sbc(sbc_res, param = "beta")
save_plot(ggplot2::last_plot(), "sbc_beta.png")
plot_sbc(sbc_res, param = 2L, bins = 10L)

# The dashed line is the expected count per bin; the dotted lines are a 99%
# band from the binomial noise of a finite n_sbc. Bars outside the band are the
# ones to look at.

# ---------------------------------------------------------------------------
# plot_coverage()
# ---------------------------------------------------------------------------
#
# The same ranks, read as intervals, one line per parameter. The invisible
# return is the expected_coverage() data frame.

cov <- plot_coverage(sbc_res)
print(round(head(cov), 3))
save_plot(ggplot2::last_plot(), "coverage.png")

# A coarser grid of levels, when the default 19 points is more than you need.
plot_coverage(sbc_res, levels = c(0.5, 0.8, 0.9, 0.95))

# ---------------------------------------------------------------------------
# plot_tarp()
# ---------------------------------------------------------------------------

tarp_res <- tarp(fit, simulator, n_tarp = 200L, n_posterior_samples = 400L,
                 seed = 1)
line <- plot_tarp(tarp_res)
print(round(head(line), 3))
save_plot(ggplot2::last_plot(), "tarp.png")

# ---------------------------------------------------------------------------
# Building on the returned objects
# ---------------------------------------------------------------------------
#
# Every function returns a ggplot (or a GGally ggmatrix), so the usual
# additions work.

p <- plot_coverage(sbc_res)          # draws, and returns the data frame
gg <- ggplot2::last_plot() +
  ggplot2::labs(title = "SIR: expected coverage",
                subtitle = sprintf("%d SBC trials, %d draws each",
                                   sbc_res$n_sbc, sbc_res$n_posterior_samples))
print(gg)
save_plot(gg, "coverage_titled.png")

# And the data is always there if you would rather draw it yourself.
df <- as.data.frame(draws)
gg2 <- ggplot2::ggplot(df, ggplot2::aes(beta, gamma)) +
  ggplot2::geom_bin2d(bins = 40) +
  ggplot2::geom_point(x = theta_true[["beta"]], y = theta_true[["gamma"]],
                      colour = "firebrick", size = 3, shape = 4, stroke = 1.5) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::theme_minimal()
print(gg2)
save_plot(gg2, "hex.png")

cat("\nplots written to", outdir, "\n")
print(list.files(outdir))

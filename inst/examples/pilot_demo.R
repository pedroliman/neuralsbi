# Pilot demo: end-to-end NPE on two tasks, saving figures to docs/figures/.
# 1) Linear Gaussian  -> posterior matches the analytic Gaussian (validation).
# 2) Two Moons        -> neural MDN recovers a bimodal, non-Gaussian posterior.
suppressMessages({ library(neuralsbi); library(torch) })
set.seed(1); torch::torch_manual_seed(1)
outdir <- "docs/figures"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Linear Gaussian ---------------------------------------------------
d <- 2; sigma <- 0.5
prior_lg <- prior_normal(mean = c(0, 0), sd = 1)
sim_lg <- function(theta) theta + matrix(rnorm(length(theta), sd = sigma),
                                         nrow = nrow(theta))
fit_lg <- npe(prior_lg, sim_lg, n_simulations = 6000,
              density_estimator = "mdn", n_components = 2L,
              hidden = c(50L, 50L), max_epochs = 300L, seed = 1)
x_obs <- c(1.0, -0.5)
post_lg <- posterior(fit_lg, x_obs = x_obs)
draws_lg <- sample(post_lg, 8000)

prec <- diag(d) + diag(d) / sigma^2
Sigma <- solve(prec); mu <- as.numeric(Sigma %*% (x_obs / sigma^2))
z <- matrix(rnorm(8000 * d), ncol = d)
analytic <- sweep(z %*% chol(Sigma), 2, mu, `+`)
cat(sprintf("Linear-Gaussian C2ST vs analytic: %.3f\n",
            c2st(draws_lg, analytic, seed = 1)$accuracy))

png(file.path(outdir, "linear_gaussian_posterior.png"), width = 900,
    height = 450, res = 110)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (j in 1:2) {
  plot(density(draws_lg[, j]), col = "steelblue", lwd = 2,
       main = bquote(theta[.(j)]), xlab = "")
  lines(density(analytic[, j]), col = "firebrick", lwd = 2, lty = 2)
  abline(v = mu[j], col = "grey40", lty = 3)
  legend("topright", c("neuralsbi", "analytic"), col = c("steelblue", "firebrick"),
         lty = c(1, 2), bty = "n", cex = 0.8)
}
par(op); dev.off()

## ---- 2. Two Moons ---------------------------------------------------------
prior_tm <- prior_uniform(low = c(-1, -1), high = c(1, 1))
two_moons_sim <- function(theta) {
  n <- nrow(theta)
  a <- runif(n, -pi / 2, pi / 2); r <- rnorm(n, 0.1, 0.01)
  cbind(r * cos(a) + 0.25 - abs(theta[, 1] + theta[, 2]) / sqrt(2),
        r * sin(a) + (-theta[, 1] + theta[, 2]) / sqrt(2))
}
fit_tm <- npe(prior_tm, two_moons_sim, n_simulations = 15000,
              density_estimator = "mdn", n_components = 8L,
              hidden = c(64L, 64L), max_epochs = 300L, seed = 1)
post_tm <- posterior(fit_tm, x_obs = c(0, 0))
draws_tm <- sample(post_tm, 5000)
proj <- (draws_tm[, 1] + draws_tm[, 2]) / sqrt(2)
cat(sprintf("Two-Moons bimodal: %s (frac |proj|>0.15 = %.2f)\n",
            mean(abs(proj) > 0.15) > 0.6, mean(abs(proj) > 0.15)))

png(file.path(outdir, "two_moons_posterior.png"), width = 500, height = 500,
    res = 110)
plot(draws_tm[, 1], draws_tm[, 2], pch = 16, cex = 0.3,
     col = adjustcolor("steelblue", 0.4), xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = expression(theta[1]), ylab = expression(theta[2]),
     main = "Two Moons posterior (x_obs = 0)")
dev.off()
cat("Saved figures to", outdir, "\n")

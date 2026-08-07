# Two Moons: a standard SBI benchmark with a bimodal, crescent-shaped posterior.
# Demonstrates that the neural MDN captures non-Gaussian, multimodal posteriors
# that the closed-form linear_gaussian baseline cannot.
suppressMessages({
  library(neuralsbi)
  library(torch)
})

set.seed(1)
torch::torch_manual_seed(1)

task <- task_two_moons()

fit <- npe(task$prior, task$simulator, n_simulations = 20000,
           density_estimator = "mdn", n_components = 8L,
           hidden = c(64L, 64L), max_epochs = 400L, seed = 1, verbose = TRUE)

x_obs <- c(0, 0)
post <- posterior(fit, x_obs = x_obs)
draws <- sample(post, 5000)

cat("\nPosterior draws summary (expect two crescent modes):\n")
cat("  theta1 range:", round(range(draws[, 1]), 2), "\n")
cat("  theta2 range:", round(range(draws[, 2]), 2), "\n")

# The bimodal direction is theta1 + theta2: the simulator's abs(theta1 + theta2)
# term makes +s and -s equally consistent with the observation, so the
# posterior has two symmetric modes at theta1 + theta2 ~ +/- c.
proj <- (draws[, 1] + draws[, 2]) / sqrt(2)
cat("  bimodal (two clusters in theta1 + theta2):",
    mean(abs(proj) > 0.15) > 0.6, "\n")

if (interactive()) {
  pairplot(draws, limits = list(c(-1, 1), c(-1, 1)))
}

library(neuralsbi)
library(deSolve)   # install.packages("deSolve")
library(future)
library(dplyr)


# Set up ------------------------------------------------------------------

# An SIR ODE:
ode_model <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dS <- -beta * S * I / N
    dI <-  beta * S * I / N - gamma * I
    dCumInfections <- beta * S * I / N
    list(c(dS, dI, dCumInfections))
  })
}

# Simulator function
# Here we are already only making the simulator return only 10 observations
my_simulator <- function(theta, times = 0:60, N = 100000, I0 = 5) {
  # Let's use just a deterministic ODE, so there's no randomness in transmission:
  # WE could make it fully stochastic if we want. 
  # I just want to make this as fast as it can be for this vignette build
  ode_result <- ode(y = c(S = N - I0, I = I0, CumInfections = 0),
                    times = times,
                    func = ode_model,
                    parms = c(beta = theta[["beta"]],
                              gamma = theta[["gamma"]],
                              N = N))
  
  cases <- round(diff(ode_result[, "CumInfections"]))
  obs_cases <- rbinom(n = length(cases), size= cases, prob = theta[["rho"]])
  
  # return output vector
  return(obs_cases)
}

# test simulator
my_simulator(theta=c(beta = 0.5, gamma = 1/7, rho = 0.5))

# Train neural posterior estimator (NPE) ----------------------------------

# Set up our priors
prior <- prior_uniform(low = c(beta = 0.1, gamma = 1/10, rho = 0.1),
                       high = c(beta = 0.7, gamma = 1/4, rho = 0.9))

# Set up parallel simulations:
plan(multisession, workers = 8)

# Observation times: the simulator returns I at these points.
times <- 0:60
sim_args <- list(times = times, N = 100000, I0 = 5)

# Here we fit the neural posterior estimator with neuralsbi's npe function:

fit   <- npe(prior, my_simulator, n_simulations = 5000, sim_args = sim_args, seed = 1)

true_r0 = 2
true_gamma = 1/7
true_params <- c(beta = true_r0 * true_gamma, gamma = true_gamma, rho = 0.7)

y_obs   <- my_simulator(true_params, times = times, N = 100000, I0 = 5)

post <- posterior(fit, x_obs = y_obs)

post_draws <- sample(post, 1000)

pairplot(post_draws, true_params)

# Condition NPE on many observations: one fit, reused --------------------

# The point of amortized NPE: train once, then condition on any number of
# observations for free. Draw 6 parameter vectors from the prior, treat each
# as a "true" data-generating process, and recover its posterior predictive.
library(ggplot2)

# Set gamma fixed and let beta vary
set.seed(42)
R0s <- c(1.5, 2, 2.5, 3)
gamma <- rep(1/7,4)
beta <- R0s * gamma
rho <- rep(0.2, 4)
theta_grid <- matrix(data = c(beta, gamma, rho), nrow = 4, dimnames = list(1:4, c("beta", "gamma", "rho")))


# For each drawn theta: simulate an observation, condition the (single) fit on
# it, then build a posterior-predictive band. Everything reuses `fit`.
pp_summary <- do.call(rbind, lapply(seq_len(nrow(theta_grid)), function(i) {
  
  theta_i <- theta_grid[i, ]
  R0_i <- R0s[i]
  y_obs   <- my_simulator(theta_i, times = times, N = 100000, I0 = 5)

  post <- posterior(fit, x_obs = y_obs)
  pp   <- posterior_predictive(post, my_simulator, 1000)

  data.frame(
    obs_id = sprintf("Case %d: R0=%.2f",
                     i, R0_i),
    time   = times[-1],
    mean   = colMeans(pp),
    lower  = apply(pp, 2, quantile, probs = 0.005),
    upper  = apply(pp, 2, quantile, probs = 0.995),
    obs    = y_obs
  )
}))

# Ribbon plot faceted 2x3, one legend, one line/point style shared across panels.
pp_plot <- ggplot(pp_summary, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "Posterior predictive 90%"),
              alpha = 0.3) +
  geom_line(aes(y = mean, colour = "Posterior predictive mean")) +
  geom_point(aes(y = obs, colour = "Observed (Ground Truth)")) +
  facet_wrap(~ obs_id, nrow = 2, ncol = 2, scales = "free_y", axes = "all") +
  scale_colour_manual(values = c("Observed (Ground Truth)" = "black",
                                 "Posterior predictive mean" = "steelblue")) +
  scale_fill_manual(values = c("Posterior predictive 99%" = "steelblue")) +
  labs(x = "Time", y = "Observed cases", colour = NULL, fill = NULL) +
  theme_classic() + 
  theme(legend.position = "bottom", strip.background = element_blank())

print(pp_plot)

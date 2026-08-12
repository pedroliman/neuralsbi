# 14_utilities.R -------------------------------------------------------------
#
# The parts that are not a method: the simulator contract, sim_args, parallel
# simulation, progress reporting, seeds and reproducibility, saving and
# reloading all three kinds of fit, and the tidy accessors on posterior draws.
#
# Task source
#   Lueckmann, J.-M., Boelts, J., Greenberg, D., Goncalves, P. and Macke, J.
#   "Benchmarking Simulation-Based Inference", AISTATS 2021.
#   https://github.com/sbi-benchmark/sbibm
#   task_gaussian_linear() is used where the point is the plumbing rather than
#   the science, because it is fast and has an analytic posterior to check
#   against.
#
#   The parallel section uses the Sick-Sicker cohort model from
#   Alarid-Escudero, F. et al. "An Introductory Tutorial on Cohort
#   State-Transition Models in R for Cost-Effectiveness Analysis", Medical
#   Decision Making 43(1), 2023,
#   https://github.com/DARTH-git/cohort-modeling-tutorial-intro/blob/main/analysis/cSTM_time_indep.R
#   run on monthly cycles, because a simulator has to be slow enough for
#   parallelism to be worth anything and a 900-cycle Markov trace is.
#
# Runtime: under a minute on a laptop CPU.

library(neuralsbi)

has_torch <- requireNamespace("torch", quietly = TRUE) &&
  torch::torch_is_installed()

# ---------------------------------------------------------------------------
# 1. The simulator contract
# ---------------------------------------------------------------------------
#
# A simulator is a function called once per parameter set, returning ONE
# simulated observation. How it receives its parameters is decided once per
# run from formals(), never by probing:
#
#   formals exactly match the prior's parameter names -> one scalar per formal
#   anything else                                     -> the named parameter
#                                                        vector as the first
#                                                        argument

prior <- prior_uniform(low = c(mu = -2, nu = -2), high = c(mu = 2, nu = 2))

named_form  <- function(mu, nu) c(a = mu + rnorm(1, sd = .2),
                                  b = nu + rnorm(1, sd = .2))
vector_form <- function(theta)  c(a = theta[["mu"]] + rnorm(1, sd = .2),
                                  b = theta[["nu"]] + rnorm(1, sd = .2))

str(simulate_for_sbi(named_form,  prior, n = 4, seed = 1)$x)
str(simulate_for_sbi(vector_form, prior, n = 4, seed = 1)$x)

# A partial name match is the trap. Some formals match and some do not, so both
# signatures are plausible, the vector form wins, and the simulator trains on
# nonsense without erroring. One typo in a prior name is enough. It warns.
half_match <- function(mu, nyu) c(a = mu, b = nyu)
invisible(tryCatch(simulate_for_sbi(half_match, prior, n = 2),
                   warning = function(w) cat("warning:", conditionMessage(w),
                                             "\n")))

# Accepted return shapes: a numeric vector, a scalar, a one-row matrix, a
# one-row data frame. Names on the output become the outcome names used in
# plots.
shapes <- list(
  vector    = function(mu, nu) c(a = mu, b = nu),
  scalar    = function(mu, nu) mu + nu,
  matrix1   = function(mu, nu) matrix(c(mu, nu), nrow = 1,
                                      dimnames = list(NULL, c("a", "b"))),
  dataframe = function(mu, nu) data.frame(a = mu, b = nu)
)
for (nm in names(shapes)) {
  cat(sprintf("%-10s -> %d column(s), names %s\n", nm,
              ncol(simulate_for_sbi(shapes[[nm]], prior, n = 2)$x),
              paste(colnames(simulate_for_sbi(shapes[[nm]], prior, n = 2)$x),
                    collapse = ",")))
}

# ---------------------------------------------------------------------------
# 2. Simulations that fail
# ---------------------------------------------------------------------------
#
# A simulation whose output is not finite is dropped together with its
# parameters, with a warning and a count.

flaky <- function(mu, nu) {
  if (mu > 1.5) return(c(a = NA_real_, b = NA_real_))
  c(a = mu + rnorm(1, sd = .2), b = nu + rnorm(1, sd = .2))
}
sims <- suppressWarnings(simulate_for_sbi(flaky, prior, n = 400, seed = 1))
cat("\ndropped:", sims$n_dropped, "of 400\n")

# Worth stopping over rather than shrugging at. Dropping conditions on the
# simulator having succeeded, and when failure depends on the parameters the
# survivors are no longer a sample from the prior: the fit then targets the
# posterior GIVEN SUCCESS. Here every draw with mu > 1.5 is gone, so the fit
# would never propose one. The honest fix is to handle the failure inside the
# simulator.

# ---------------------------------------------------------------------------
# 3. sim_args: everything the simulator needs that is not a parameter
# ---------------------------------------------------------------------------
#
# Observed covariates, a time grid, a population size, solver settings. Passed
# to every call. Prefer this to closing over a large object: under a parallel
# plan the simulator and everything its environment captures is shipped to a
# worker once per batch.

with_grid <- function(mu, nu, times, noise) {
  stats::setNames(mu * exp(-abs(nu) * times) + rnorm(length(times), sd = noise),
                  paste0("t", times))
}
sims <- simulate_for_sbi(with_grid, prior, n = 5,
                         sim_args = list(times = c(0.5, 1, 2, 4),
                                         noise = 0.05),
                         seed = 1)
print(round(sims$x, 3))

# Every function that calls a simulator takes sim_args: npe(), nle(), nre(),
# npe_sequential(), sbc(), tarp(), posterior_predictive().

# ---------------------------------------------------------------------------
# 4. Parallel simulation
# ---------------------------------------------------------------------------
#
# Declare a future plan and every simulator call in the package spreads across
# workers. There is no argument to set and no variant of the function to call.

# The Sick-Sicker model again, on monthly rather than annual cycles (the DARTH
# tutorial's own alternative, cycle_length = 1/12). 900 cycles per call is slow
# enough for the comparison to mean something.
cycle_length <- 1 / 12
n_cycles <- 900
rate_to_prob <- function(r, t = 1) 1 - exp(-r * t)

sick_sicker <- function(p_S1S2, hr_S1, hr_S2) {
  p_HS1 <- rate_to_prob(0.15, cycle_length)
  p_S1H <- rate_to_prob(0.5, cycle_length)
  p_HD  <- rate_to_prob(0.002, cycle_length)
  p_S1D <- rate_to_prob(0.002 * hr_S1, cycle_length)
  p_S2D <- rate_to_prob(0.002 * hr_S2, cycle_length)
  P <- rbind(c((1 - p_HD) * (1 - p_HS1), (1 - p_HD) * p_HS1, 0, p_HD),
             c((1 - p_S1D) * p_S1H, (1 - p_S1D) * (1 - p_S1H - p_S1S2),
               (1 - p_S1D) * p_S1S2, p_S1D),
             c(0, 0, 1 - p_S2D, p_S2D),
             c(0, 0, 0, 1))
  m <- c(1, 0, 0, 0)
  keep <- c(120, 300, 480, 660)          # ages 35, 50, 65, 80
  out <- numeric(0)
  for (t in seq_len(n_cycles)) {
    m <- m %*% P
    if (t %in% keep) out <- c(out, 1 - m[4], m[2] + m[3])
  }
  stats::setNames(out, paste0(rep(c("surv", "prev"), 4),
                              rep(c(35, 50, 65, 80), each = 2)))
}

ss_prior <- prior_uniform(low  = c(p_S1S2 = 0.001, hr_S1 = 1.0, hr_S2 =  5),
                          high = c(p_S1S2 = 0.050, hr_S1 = 4.5, hr_S2 = 15))

if (requireNamespace("future", quietly = TRUE)) {
  n_sim <- 8000
  t_seq <- system.time(
    simulate_for_sbi(sick_sicker, ss_prior, n = n_sim, seed = 1)
  )[["elapsed"]]

  old_plan <- future::plan()
  workers <- max(2L, min(4L, parallel::detectCores()))
  future::plan(future::multisession, workers = workers)
  t_par <- system.time(
    simulate_for_sbi(sick_sicker, ss_prior, n = n_sim, seed = 1)
  )[["elapsed"]]
  future::plan(old_plan)

  cat(sprintf("\n%d simulations: %.1f s sequential, %.1f s on %d workers\n",
              n_sim, t_seq, t_par, workers))
} else {
  cat("\nfuture is not installed; simulation stays sequential.\n")
}

# Expect roughly a 2x speedup on four workers, not 4x. A multisession plan
# starts fresh R processes and ships the simulator to each one, which costs a
# couple of seconds before any simulation runs. Below about ten seconds of
# total simulation time you will lose money on the trade; the crossover is
# worth measuring once for your own simulator rather than assuming.

# Random numbers do not depend on the plan. Each simulation gets its own
# L'Ecuyer-CMRG stream derived from the session RNG state at the moment
# simulation starts, so set.seed(42) (or seed = 42) gives the same simulations
# whatever the worker count.
a <- simulate_for_sbi(sick_sicker, ss_prior, n = 50, seed = 42)$x
b <- simulate_for_sbi(sick_sicker, ss_prior, n = 50, seed = 42)$x
cat("same seed, same simulations:", identical(a, b), "\n")

# Training always runs in the calling process: torch modules are external
# pointers and cannot be shipped to a worker. libtorch threads internally,
# through torch::torch_set_num_threads().

# The one-per-session hint about going parallel is silenced with
options(neuralsbi.parallel_hint = FALSE)

# ---------------------------------------------------------------------------
# 5. Progress reporting
# ---------------------------------------------------------------------------
#
# options(neuralsbi.progress = ):
#   "auto"     (default) report interactively, stay quiet in scripts and knitr
#   TRUE       always report
#   FALSE      never
#   "builtin"  always report with the built-in bar, ignoring progressr
#
# With progressr installed, the package emits standard progressr updates, so
# progressr::handlers() controls the look:
#
#   library(progressr)
#   handlers(global = TRUE)
#   handlers("cli")
#   fit <- npe(prior, simulator, n_simulations = 10000)

old_progress <- getOption("neuralsbi.progress")
options(neuralsbi.progress = FALSE)

# While a neural estimator trains, one step is one epoch. The total is not
# known in advance because training stops early, so the bar targets
# best epoch + patience and the target moves out whenever the loss improves.

# ---------------------------------------------------------------------------
# 6. Seeds and reproducibility
# ---------------------------------------------------------------------------
#
# Every fitting function takes seed. It covers the simulations AND the network
# initialization, so two calls with the same seed give the same fit.

task <- task_gaussian_linear(dim = 3L)
f1 <- npe(task$prior, task$simulator, n_simulations = 1000,
          density_estimator = "linear_gaussian", seed = 7)
f2 <- npe(task$prior, task$simulator, n_simulations = 1000,
          density_estimator = "linear_gaussian", seed = 7)
x_obs <- task$simulator(sample_prior(task$prior, 1)[1, ])

# Sampling is separately random, so pin the session seed before each draw or
# the comparison tests the sampler rather than the fit.
set.seed(1); d1 <- sample(posterior(f1, x_obs = x_obs), 100)
set.seed(1); d2 <- sample(posterior(f2, x_obs = x_obs), 100)
cat("\nsame seed, same fit:",
    isTRUE(all.equal(d1, d2, check.attributes = FALSE)), "\n")

set.seed(1); d3 <- sample(posterior(f1, x_obs = x_obs), 100)
d4 <- sample(posterior(f1, x_obs = x_obs), 100)
cat("no session seed, different draws from the same fit:",
    !isTRUE(all.equal(d3, d4, check.attributes = FALSE)), "\n")

# ---------------------------------------------------------------------------
# 7. Saving and reloading
# ---------------------------------------------------------------------------
#
# saveRDS() on a torch-backed fit writes the external pointer, not the network.
# The file reloads without complaint, the object prints normally, and the first
# call that touches the network fails with "external pointer is not valid".
# save_npe() writes the weights via torch_save() on the state_dict alongside
# everything else the fit carries as ordinary R objects.
#
# save_nle()/load_nle() and save_nre()/load_nre() are aliases. Every pair
# handles every kind of fit; the extra names exist so calling code reads the
# way the fit was made.

fit_npe <- npe(task$prior, task$simulator, n_simulations = 1000,
               density_estimator = if (has_torch) "maf" else "linear_gaussian",
               seed = 1)
fit_nle <- nle(task$prior, task$simulator, n_simulations = 1000,
               density_estimator = "linear_gaussian", seed = 1)
fit_nre <- nre(task$prior, task$simulator, n_simulations = 1000,
               classifier = "logistic", seed = 1)

dir <- file.path(tempdir(), "nsbi-fits")
dir.create(dir, showWarnings = FALSE)
save_npe(fit_npe, file.path(dir, "npe.rds"))
save_nle(fit_nle, file.path(dir, "nle.rds"))
save_nre(fit_nre, file.path(dir, "nre.rds"))
print(file.size(list.files(dir, full.names = TRUE)))

back <- load_npe(file.path(dir, "npe.rds"))
print(back)
cat("reloaded posterior works:",
    nrow(sample(posterior(back, x_obs = x_obs), 50)) == 50, "\n")

# What goes into the file: a format marker, the package version that wrote it,
# a timestamp, the fit with its network dropped, and the weights as raw bytes.
# Weights, not code, so a fit saved by one version loads into a later one as
# long as the architecture has not changed. A fit trained on "cuda" or "mps"
# always reloads onto CPU; move it back with fit$de$net$to(device = "cuda").
bundle <- readRDS(file.path(dir, "npe.rds"))
str(bundle, max.level = 1)
unlink(dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# 8. Tidy accessors
# ---------------------------------------------------------------------------

post <- posterior(fit_npe, x_obs = x_obs)
draws <- sample(post, 2000)

print(draws)                       # short report, not 2000 rows
print(summary(draws))              # one row per parameter
print(head(as.data.frame(draws)))  # straight into dplyr or ggplot2

# summary() on the posterior samples for you.
print(summary(post, n = 500))

# summary() on a fit prints it and returns the training metadata invisibly.
str(summary(fit_npe))
str(summary(fit_nre))              # records `classifier`, not `density_estimator`

# sample_posterior() is sample() under a name that does not mask base::sample.
cat("\n", nrow(sample_posterior(post, n = 100)), "draws\n")

# ---------------------------------------------------------------------------
# 9. The device argument
# ---------------------------------------------------------------------------
#
#   "cpu"          the default, on purpose, matching Python sbi
#   "cuda", "mps"  error if that device is not actually available, so a typo or
#                  a missing driver is not mistaken for a slow CPU run
#   "gpu", "auto"  resolve CUDA -> MPS -> CPU, falling back quietly
#
# For the small networks typical of SBI, a GPU is often SLOWER than CPU:
# per-kernel launch overhead dominates until the network and batch are large.
# Try CPU first and switch only if profiling says otherwise.

bad <- try(npe(task$prior, task$simulator, n_simulations = 100,
               density_estimator = "linear_gaussian", device = "gpu-ish"),
           silent = TRUE)
cat(conditionMessage(attr(bad, "condition")), "\n")

# ---------------------------------------------------------------------------
# 10. Prior helpers, once more
# ---------------------------------------------------------------------------

print(round(sample_prior(ss_prior, 3), 3))
print(within_support(ss_prior, rbind(c(0.2, 3, 10), c(0.2, 3, 99))))

options(neuralsbi.progress = old_progress)

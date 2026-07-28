#!/usr/bin/env Rscript
# Fast structural checks for the benchmark harness. Minutes, not hours.
#
# It does not measure anything: it checks that each piece is wired to the right
# thing. The interesting checks are the ones that can actually fail silently in
# a real run:
#
#   * every simulator produces data of the right shape;
#   * every simulator is *calibrated against sbibm's own observations* -- we
#     simulate at the true parameters that generated each shipped observation
#     and check the observation looks like a draw from our simulator;
#   * the Bernoulli GLM summary statistics reproduce sbibm's, exactly, from its
#     shipped raw spike train;
#   * the C2ST scores two halves of a reference posterior at chance and a prior
#     sample near 1;
#   * the paper's budget labels decode to the right numbers;
#   * a whole cell runs end to end, for both NPE and NLE.
#
#   Rscript dev/benchmarks/smoke_test.R
#   Rscript dev/benchmarks/smoke_test.R --tasks two_moons,sir

.here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) normalizePath(dirname(sub("^--file=", "", a[1])))
  else normalizePath("dev/benchmarks")
})
for (f in c("utils.R", "pt_io.R", "sbibm_data.R", "tasks.R", "c2st.R",
            "runner.R")) {
  source(file.path(.here, "R", f))
}

opts <- parse_args(list(
  tasks = sbibm_task_names(),
  n_replicates = 200,
  end_to_end = TRUE
))

load_neuralsbi()
set.seed(20240101)

PASS <- 0L; FAIL <- 0L
check <- function(label, ok, detail = "") {
  if (isTRUE(ok)) {
    PASS <<- PASS + 1L
    cat(sprintf("  ok    %-58s %s\n", label, detail))
  } else {
    FAIL <<- FAIL + 1L
    cat(sprintf("  FAIL  %-58s %s\n", label, detail))
  }
}

# --- 1. published results decode correctly -----------------------------------

cat("\npaper results\n")
res <- paper_results()
check("main_paper.csv has 2400 rows", nrow(res) == 2400)
check("three budgets decode to 1e3/1e4/1e5",
      identical(sort(unique(res$num_simulations)), c(1e3, 1e4, 1e5)))
# Sanity on the decoding: C2ST must improve with budget for gaussian_linear NPE.
gl <- vapply(c(1e3, 1e4, 1e5), function(n) {
  stats::median(paper_c2st("gaussian_linear", "NPE", n, res)$C2ST)
}, numeric(1))
check("gaussian_linear NPE C2ST falls with budget",
      all(diff(gl) < 0), paste(sprintf("%.3f", gl), collapse = " > "))

# --- 2. frozen constants read out of sbibm's .pt files ------------------------

cat("\nfrozen constants\n")
X <- glm_design_matrix()
check("GLM design matrix is 100 x 10", identical(dim(X), c(100L, 10L)))
check("GLM design matrix column 1 is the intercept", all(X[, 1] == 1))
check("GLM design matrix columns 2:10 are lagged copies of one stimulus",
      all(abs(X[2:100, 3] - X[1:99, 2]) < 1e-6) && X[1, 3] == 0)

perm <- slcp_permutation()
check("SLCP distractor permutation is a permutation of 1:100",
      setequal(perm, 1:100))
gmm <- slcp_noise_params()
check("SLCP distractor mixture has 20 components in 92 dimensions",
      identical(dim(gmm$loc), c(20L, 92L)) && length(gmm$scale_tril) == 20 &&
        identical(dim(gmm$scale_tril[[1]]), c(92L, 92L)))
check("SLCP distractor scale factors are lower triangular",
      all(abs(gmm$scale_tril[[1]][upper.tri(gmm$scale_tril[[1]])]) < 1e-12))

# --- 3. the GLM summary statistics are exactly sbibm's ------------------------

cat("\nbernoulli_glm summary statistics\n")
glm_task <- sbibm_task("bernoulli_glm")
glm_raw <- sbibm_task("bernoulli_glm_raw")
worst <- 0
for (i in 1:10) {
  y <- sbibm_observation(glm_raw, i)
  want <- as.numeric(sbibm_observation(glm_task, i))
  got <- as.numeric(y %*% X)
  worst <- max(worst, max(abs(got - want)))
}
check("summary(raw spike train) reproduces sbibm's observation.csv",
      worst < 1e-4, sprintf("max abs error %.2e over 10 observations", worst))

# --- 4. every simulator: shapes, and agreement with sbibm's observations -----

cat("\nsimulators\n")
tasks <- intersect(sbibm_task_names(), opts$tasks)
for (nm in tasks) {
  tk <- sbibm_task(nm)
  th <- sample_prior(tk$prior, 32)
  x <- tk$simulate(th)
  check(sprintf("%s: prior draws are %d-dimensional", nm, tk$dim_theta),
        ncol(th) == tk$dim_theta)
  check(sprintf("%s: simulator returns %d columns", nm, tk$dim_x),
        is.matrix(x) && ncol(x) == tk$dim_x && nrow(x) == 32)

  obs1 <- sbibm_observation(tk, 1)
  check(sprintf("%s: sbibm observation has matching width", nm),
        ncol(obs1) == tk$dim_x)

  # Simulate replicates at the parameters that generated each observation and
  # ask where the observation falls in the simulated marginals. If our
  # simulator matched sbibm's, those ranks are uniform; if it did not, they pile
  # up at the ends. Mid-ranks, because the counts in sir and the spike trains in
  # bernoulli_glm_raw are discrete and heavily tied at zero.
  ranks <- c()
  for (i in 1:10) {
    theta_true <- sbibm_true_parameters(tk, i)
    obs <- as.numeric(sbibm_observation(tk, i))
    rep_theta <- matrix(rep(as.numeric(theta_true), each = opts$n_replicates),
                        nrow = opts$n_replicates)
    sims <- tk$simulate(rep_theta)
    sims <- sims[is.finite(rowSums(sims)), , drop = FALSE]
    if (!nrow(sims)) next
    below <- colMeans(sweep(sims, 2, obs, `<`))
    tied <- colMeans(sweep(sims, 2, obs, `==`))
    ranks <- c(ranks, below + 0.5 * tied)
  }
  # The tail band cannot be finer than the replicate grid, so it and the
  # threshold both scale with `n_replicates`: under uniform ranks a fraction
  # 2 * band lands outside, and we allow three times that.
  band <- max(0.005, 1 / opts$n_replicates)
  outside <- mean(ranks < band | ranks > 1 - band)
  check(sprintf("%s: sbibm observations look like our simulator's draws", nm),
        outside < max(0.10, 6 * band) && abs(mean(ranks) - 0.5) < 0.15,
        sprintf("%.1f%% of coordinates in the outer %.1f%%, mean rank %.2f",
                100 * outside, 100 * band, mean(ranks)))
}

# --- 5. reference posteriors and the C2ST ------------------------------------

cat("\nc2st\n")
tm <- sbibm_task("two_moons")
ref <- sbibm_reference_posterior(tm, 1)
check("two_moons reference posterior is 10000 x 2",
      identical(dim(ref), c(10000L, 2L)))

half <- 1500
a <- ref[1:half, , drop = FALSE]
b <- ref[(half + 1):(2 * half), , drop = FALSE]
t0 <- Sys.time()
null_c2st <- c2st_sbibm(a, b, seed = 1)
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
check("C2ST of two halves of the reference posterior is at chance",
      abs(null_c2st - 0.5) < 0.05,
      sprintf("%.3f (%.0fs for 2x%d samples)", null_c2st, dt, half))

prior_draws <- sample_prior(tm$prior, half)
alt_c2st <- c2st_sbibm(a, prior_draws, seed = 1)
check("C2ST of reference posterior vs. prior is near 1",
      alt_c2st > 0.9, sprintf("%.3f", alt_c2st))

# --- 6. end to end -----------------------------------------------------------

if (isTRUE(opts$end_to_end)) {
  cat("\nend to end (linear_gaussian estimator, so no torch needed)\n")
  cell <- run_cell("gaussian_linear", "npe", num_simulations = 500,
                   num_observation = 1, seed = 1,
                   num_posterior_samples = 1000, estimator = "linear_gaussian",
                   verbose = FALSE)
  check("NPE cell returns 1000 posterior draws",
        identical(dim(cell$samples), c(1000L, 10L)))
  cell <- score_cell(cell)
  # linear_gaussian is exact for this task, so even 500 simulations should land
  # close to the analytic reference.
  check("NPE on gaussian_linear scores near chance against the reference",
        cell$C2ST < 0.7, sprintf("C2ST %.3f", cell$C2ST))
  row <- cell_row(cell)
  check("metrics row has the expected columns",
        all(c("task", "algorithm", "num_simulations", "num_observation",
              "C2ST") %in% names(row)))

  cell <- run_cell("gaussian_linear", "nle", num_simulations = 500,
                   num_observation = 1, seed = 1, num_posterior_samples = 500,
                   estimator = "linear_gaussian", verbose = FALSE,
                   mcmc = list(n_chains = 10L, warmup = 50L, thin = 2L,
                               init_strategy = "resample"))
  check("NLE cell returns 500 posterior draws",
        identical(dim(cell$samples), c(500L, 10L)))
  cell <- score_cell(cell)
  check("NLE on gaussian_linear scores near chance against the reference",
        cell$C2ST < 0.75, sprintf("C2ST %.3f", cell$C2ST))

  # The real configuration, if torch is here. Deliberately undertrained: this
  # asks whether nsf and maf are wired up, not how well they do.
  if (requireNamespace("torch", quietly = TRUE) && torch::torch_is_installed()) {
    cat("\nend to end (the paper's estimators, undertrained on purpose)\n")
    for (alg in c("npe", "nle")) {
      cell <- run_cell("two_moons", alg, num_simulations = 200,
                       num_observation = 1, seed = 1,
                       num_posterior_samples = 200, max_epochs = 15L,
                       verbose = FALSE,
                       mcmc = list(n_chains = 10L, warmup = 20L, thin = 1L,
                                   init_strategy = "resample"))
      check(sprintf("%s runs with the paper's estimator (%s)", toupper(alg),
                    if (alg == "npe") "nsf" else "maf"),
            identical(dim(cell$samples), c(200L, 2L)),
            sprintf("%.0fs", cell$time_total))
    }
  } else {
    cat("\ntorch not available: skipping the nsf/maf wiring check\n")
  }
}

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)

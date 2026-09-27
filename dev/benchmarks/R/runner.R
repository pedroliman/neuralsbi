# Running one cell of the benchmark grid with neuralsbi, under the same
# settings sbibm used.
#
# The hyperparameters below come from `config/algorithm/{npe,nle}.yaml` in the
# sbi-benchmark results repository, which is what produced `main_paper.csv`:
#
#   NPE: neural_net = nsf,  hidden_features = 50, training_batch_size = 10000,
#        z_score_x = z_score_theta = true, automatic_transforms_enabled = false,
#        num_rounds = 1
#   NLE: neural_net = maf,  hidden_features = 50, training_batch_size = 10000,
#        slice sampling with 100 chains, thin 10, 100 warmup steps
#
# `sbi`'s flow defaults fill in the rest: 5 transforms, 10 spline bins, tail
# bound 3, Adam at 5e-4, 10% validation split, early stopping after 20 epochs
# without improvement. neuralsbi shares those defaults, so only the estimator,
# the hidden width and the batch size need setting here.

NUM_POSTERIOR_SAMPLES <- 10000L

npe_settings <- function(num_simulations) {
  list(density_estimator = "nsf",
       hidden = c(50L, 50L),
       n_transforms = 5L,
       n_bins = 10L,
       tail_bound = 3,
       batch_size = as.integer(min(10000, num_simulations)))
}

nle_settings <- function(num_simulations) {
  list(density_estimator = "maf",
       hidden = c(50L, 50L),
       n_transforms = 5L,
       batch_size = as.integer(min(10000, num_simulations)))
}

# sbibm's `mcmc_parameters` for (S)NLE, mapped onto neuralsbi's slice sampler.
NLE_MCMC <- list(n_chains = 100L, warmup = 100L, thin = 10L,
                 init_strategy = "resample")

#' Run one (task, algorithm, budget, observation) cell.
#'
#' @return A list with the posterior `samples`, timings, and the metadata the
#'   report needs.
#' @param estimator Overrides the paper's estimator choice. Only the smoke test
#'   uses this, to exercise the plumbing with the torch-free `linear_gaussian`
#'   baseline; a real benchmark run must leave it alone.
run_cell <- function(task_name, algorithm, num_simulations, num_observation,
                     seed = 1L, num_posterior_samples = NUM_POSTERIOR_SAMPLES,
                     max_epochs = 2000L, verbose = TRUE, estimator = NULL,
                     mcmc = NLE_MCMC) {
  algorithm <- match.arg(tolower(algorithm), c("npe", "nle"))
  task <- sbibm_task(task_name)
  num_simulations <- as.integer(num_simulations)

  if (verbose) {
    say(sprintf("%s / %s / %s sims / observation %d",
                task_name, toupper(algorithm), budget_label(num_simulations),
                num_observation))
  }

  set.seed(seed)
  t0 <- Sys.time()
  theta <- sample_prior(task$prior, num_simulations)
  x <- simulate_chunked(task, theta, verbose = verbose && num_simulations > 20000)
  t_sim <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (verbose) say(sprintf("  simulated %d draws in %.1fs", num_simulations, t_sim))

  obs <- sbibm_observation(task, num_observation)

  t1 <- Sys.time()
  if (algorithm == "npe") {
    cfg <- npe_settings(num_simulations)
    if (!is.null(estimator)) cfg <- list(density_estimator = estimator)
    fit <- do.call(npe, c(list(prior = task$prior, theta = theta, x = x,
                               max_epochs = max_epochs, seed = seed,
                               verbose = verbose), cfg))
    t_train <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    t2 <- Sys.time()
    post <- posterior(fit, x_obs = obs)
    samples <- sample(post, n = num_posterior_samples)
    acceptance <- attr(samples, "acceptance_rate")
  } else {
    cfg <- nle_settings(num_simulations)
    if (!is.null(estimator)) cfg <- list(density_estimator = estimator)
    fit <- do.call(nle, c(list(prior = task$prior, theta = theta, x = x,
                               max_epochs = max_epochs, seed = seed,
                               verbose = verbose), cfg))
    t_train <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    t2 <- Sys.time()
    post <- do.call(posterior, c(list(fit = fit, x_obs = obs, sampler = "slice",
                                      seed = seed), mcmc))
    samples <- sample(post, n = num_posterior_samples, verbose = verbose)
    acceptance <- NA_real_
  }
  t_sample <- as.numeric(difftime(Sys.time(), t2, units = "secs"))

  list(samples = unname(as.matrix(samples)),
       task = task_name, algorithm = toupper(algorithm),
       num_simulations = num_simulations, num_observation = num_observation,
       seed = seed,
       n_dropped = fit$n_dropped %||% 0L,
       acceptance_rate = acceptance,
       time_simulate = t_sim, time_train = t_train, time_sample = t_sample,
       time_total = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

#' Score a cell's posterior samples against sbibm's reference posterior.
score_cell <- function(cell, c2st_seed = 1L, max_epochs = 1000L) {
  task <- sbibm_task(cell$task)
  reference <- sbibm_reference_posterior(task, cell$num_observation)
  value <- c2st_sbibm(X = reference, Y = cell$samples, z_score = TRUE,
                      seed = c2st_seed, max_epochs = max_epochs)
  cell$C2ST <- as.numeric(value)
  cell$C2ST_folds <- paste(sprintf("%.4f", attr(value, "folds")), collapse = ";")
  cell
}

#' Flatten a scored cell into a one-row data frame for `metrics.csv`.
cell_row <- function(cell) {
  data.frame(
    task = cell$task,
    algorithm = cell$algorithm,
    num_simulations = cell$num_simulations,
    num_observation = cell$num_observation,
    seed = cell$seed,
    C2ST = cell$C2ST %||% NA_real_,
    C2ST_folds = cell$C2ST_folds %||% NA_character_,
    n_dropped = cell$n_dropped,
    acceptance_rate = cell$acceptance_rate,
    time_simulate = cell$time_simulate,
    time_train = cell$time_train,
    time_sample = cell$time_sample,
    time_total = cell$time_total,
    neuralsbi_version = as.character(utils::packageVersion("neuralsbi")),
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    stringsAsFactors = FALSE
  )
}

results_dir <- function() {
  d <- Sys.getenv("NEURALSBI_BENCH_OUT", "")
  if (!nzchar(d)) d <- file.path(BENCH_DIR, "results")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

metrics_path <- function() file.path(results_dir(), "metrics.csv")

#' Append a row to `metrics.csv`, replacing any earlier row for the same cell.
append_metrics <- function(row) {
  path <- metrics_path()
  if (file.exists(path)) {
    old <- utils::read.csv(path, stringsAsFactors = FALSE)
    key <- function(d) paste(d$task, d$algorithm, d$num_simulations,
                             d$num_observation, sep = "|")
    old <- old[!key(old) %in% key(row), , drop = FALSE]
    row <- rbind(old, row)
  }
  utils::write.csv(row, path, row.names = FALSE)
  invisible(path)
}

#' Store posterior draws so a run can be re-scored without re-fitting.
save_samples <- function(cell) {
  d <- file.path(results_dir(), "samples", cell$task)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(d, sprintf("%s_%s_obs%02d.csv.gz", tolower(cell$algorithm),
                               format(cell$num_simulations, scientific = FALSE,
                                      trim = TRUE),
                               cell$num_observation))
  con <- gzfile(path, "w")
  utils::write.csv(cell$samples, con, row.names = FALSE)
  close(con)
  invisible(path)
}

# Access to the reference material shipped with sbibm and with the
# `sbi-benchmark/results` repository.
#
# We never re-derive the ground truth. Observations, true parameters and
# reference posterior samples come straight from the sbibm checkout, and the
# published C2ST numbers come straight from `main_paper.csv`. That is the whole
# point: neuralsbi is scored against exactly the same targets the paper used.

#' Root of the sbibm checkout (the directory containing `sbibm/tasks`).
sbibm_root <- function() {
  p <- Sys.getenv("SBIBM_PATH", "")
  if (!nzchar(p)) p <- file.path(BENCH_DIR, "external", "sbibm")
  if (!dir.exists(file.path(p, "sbibm", "tasks"))) {
    stop("sbibm not found at '", p, "'.\n",
         "Run `Rscript dev/benchmarks/00_setup.R` or set SBIBM_PATH.",
         call. = FALSE)
  }
  normalizePath(p)
}

#' Root of the sbi-benchmark results checkout.
results_root <- function() {
  p <- Sys.getenv("SBIBM_RESULTS_PATH", "")
  if (!nzchar(p)) p <- file.path(BENCH_DIR, "external", "results")
  if (!file.exists(file.path(p, "benchmarking_sbi", "results", "main_paper.csv"))) {
    stop("sbi-benchmark results not found at '", p, "'.\n",
         "Run `Rscript dev/benchmarks/00_setup.R` or set SBIBM_RESULTS_PATH.",
         call. = FALSE)
  }
  normalizePath(p)
}

task_files_dir <- function(task) {
  file.path(sbibm_root(), "sbibm", "tasks", task$files_dir, "files")
}

obs_dir <- function(task, num_observation) {
  file.path(task_files_dir(task), paste0("num_observation_", num_observation))
}

#' The observation sbibm conditions on, as a 1 x dim_x matrix.
sbibm_observation <- function(task, num_observation) {
  path <- file.path(obs_dir(task, num_observation), task$observation_file)
  as.matrix(utils::read.csv(path))
}

#' The parameter that generated the observation, as a 1 x dim_theta matrix.
sbibm_true_parameters <- function(task, num_observation) {
  path <- file.path(obs_dir(task, num_observation), "true_parameters.csv")
  as.matrix(utils::read.csv(path))
}

#' The 10k reference posterior samples for an observation.
sbibm_reference_posterior <- function(task, num_observation) {
  path <- file.path(obs_dir(task, num_observation),
                    "reference_posterior_samples.csv.bz2")
  as.matrix(utils::read.csv(bzfile(path)))
}

# --- published results -------------------------------------------------------

# main_paper.csv writes budgets as "10" followed by a unicode superscript digit
# (U+00B3, U+2074, U+2075). Decoding them by their final byte sidesteps every
# locale and encoding question a literal comparison would raise.
decode_budget <- function(labels) {
  last_byte <- vapply(labels, function(s) {
    r <- charToRaw(s)
    as.integer(r[length(r)])
  }, integer(1), USE.NAMES = FALSE)
  out <- rep(NA_real_, length(labels))
  out[last_byte == 0xb3] <- 1e3
  out[last_byte == 0xb4] <- 1e4
  out[last_byte == 0xb5] <- 1e5
  out
}

#' The paper's per-run results, tidied: budgets as numbers, algorithms upper case.
paper_results <- function() {
  path <- file.path(results_root(), "benchmarking_sbi", "results",
                    "main_paper.csv")
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  d$num_simulations <- decode_budget(d$num_simulations)
  if (anyNA(d$num_simulations)) {
    stop("Unrecognised simulation-budget labels in main_paper.csv", call. = FALSE)
  }
  d
}

#' The paper's C2ST for one cell of the grid, over all 10 observations.
#'
#' Returns a data frame with `num_observation` and `C2ST`.
paper_c2st <- function(task_name, algorithm, num_simulations,
                       results = paper_results()) {
  sel <- results$task == task_name &
    results$algorithm == toupper(algorithm) &
    results$num_simulations == num_simulations
  out <- results[sel, c("num_observation", "C2ST")]
  out[order(out$num_observation), ]
}

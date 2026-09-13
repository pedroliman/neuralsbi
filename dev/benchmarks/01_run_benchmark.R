#!/usr/bin/env Rscript
# Run part of the benchmark grid with neuralsbi and score it against sbibm's
# reference posteriors.
#
# Every argument accepts a comma-separated list, and the defaults are the full
# grid the paper reports for NPE and NLE: 10 tasks x 2 algorithms x 3 budgets x
# 10 observations = 600 fits. Start smaller.
#
#   Rscript 01_run_benchmark.R --tasks two_moons --algorithms npe \
#       --budgets 1000 --observations 1
#
#   Rscript 01_run_benchmark.R --tasks gaussian_linear,two_moons,slcp \
#       --algorithms npe,nle --budgets 1000,10000 --observations 1,2,3
#
# Results are appended to results/metrics.csv, one row per cell, replacing any
# earlier row for the same cell. Posterior draws go to results/samples/ so a run
# can be re-scored without re-fitting. Interrupt and resume freely: with
# --skip-existing the script leaves finished cells alone.

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
  algorithms = c("npe", "nle"),
  budgets = c(1e3, 1e4, 1e5),
  observations = 1:10,
  seed = 1,
  posterior_samples = 10000,
  max_epochs = 2000,
  c2st_max_epochs = 1000,
  skip_existing = FALSE,
  save_samples = TRUE,
  verbose = TRUE,
  estimator = ""
))

load_neuralsbi()

estimator <- if (nzchar(opts$estimator)) opts$estimator else NULL
if (is.null(estimator)) {
  if (!requireNamespace("torch", quietly = TRUE) || !torch::torch_is_installed()) {
    stop("This benchmark needs torch: the paper's estimators are nsf (NPE) and ",
         "maf (NLE), both torch-backed. Run torch::install_torch().\n",
         "To exercise the pipeline without torch, pass ",
         "--estimator linear_gaussian; the resulting numbers are not a ",
         "reproduction of the paper.", call. = FALSE)
  }
} else {
  say("NOTE: --estimator ", estimator, " overrides the paper's estimator. ",
      "These numbers do not reproduce the paper.")
}

done <- if (opts$skip_existing && file.exists(metrics_path())) {
  d <- utils::read.csv(metrics_path(), stringsAsFactors = FALSE)
  paste(d$task, toupper(d$algorithm), d$num_simulations, d$num_observation,
        sep = "|")
} else character(0)

grid <- expand.grid(num_observation = as.integer(opts$observations),
                    num_simulations = as.integer(opts$budgets),
                    algorithm = tolower(opts$algorithms),
                    task = opts$tasks,
                    stringsAsFactors = FALSE)
grid <- grid[, c("task", "algorithm", "num_simulations", "num_observation")]

say(sprintf("%d cells to run", nrow(grid)))
failures <- list()

for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  key <- paste(g$task, toupper(g$algorithm), g$num_simulations,
               g$num_observation, sep = "|")
  if (key %in% done) {
    say(sprintf("[%d/%d] skipping (already done) %s", i, nrow(grid), key))
    next
  }
  say(sprintf("[%d/%d] %s", i, nrow(grid), key))
  ok <- tryCatch({
    cell <- run_cell(g$task, g$algorithm, g$num_simulations, g$num_observation,
                     seed = as.integer(opts$seed),
                     num_posterior_samples = as.integer(opts$posterior_samples),
                     max_epochs = as.integer(opts$max_epochs),
                     verbose = isTRUE(opts$verbose), estimator = estimator)
    if (isTRUE(opts$save_samples)) save_samples(cell)
    cell <- score_cell(cell, c2st_seed = as.integer(opts$seed),
                       max_epochs = as.integer(opts$c2st_max_epochs))
    append_metrics(cell_row(cell))
    say(sprintf("  C2ST = %.3f  (%.0fs)", cell$C2ST, cell$time_total))
    TRUE
  }, error = function(e) {
    say("  FAILED: ", conditionMessage(e))
    failures[[length(failures) + 1L]] <<- list(key = key,
                                               message = conditionMessage(e))
    FALSE
  })
}

if (length(failures)) {
  say(sprintf("%d cell(s) failed:", length(failures)))
  for (f in failures) say("  ", f$key, ": ", f$message)
}
say("metrics written to ", metrics_path())
say("now run: Rscript 03_report.R")

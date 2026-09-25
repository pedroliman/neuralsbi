#!/usr/bin/env Rscript
# Drive the whole grid, one task at a time, in a separate R process each.
#
# Two reasons for the subprocesses: a task that blows up (an ODE that will not
# integrate, an out-of-memory flow) does not take the rest of the run with it,
# and torch releases its memory between tasks. Each subprocess is
# `01_run_benchmark.R` with --skip-existing, so the whole thing is resumable:
# rerun it after an interruption and it picks up where it stopped.
#
# Tasks run cheapest first, so problems surface in minutes rather than hours.
#
#   Rscript 02_run_all.R                         # everything (days of compute)
#   Rscript 02_run_all.R --budgets 1000 --observations 1,2,3
#   Rscript 02_run_all.R --algorithms npe

.here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) normalizePath(dirname(sub("^--file=", "", a[1])))
  else normalizePath("dev/benchmarks")
})
for (f in c("utils.R", "pt_io.R", "sbibm_data.R", "tasks.R")) {
  source(file.path(.here, "R", f))
}

# Roughly ascending cost: low-dimensional and analytic first, ODEs and
# 100-dimensional data last.
TASK_ORDER <- c("gaussian_linear", "gaussian_linear_uniform", "gaussian_mixture",
                "two_moons", "slcp", "bernoulli_glm", "sir", "lotka_volterra",
                "bernoulli_glm_raw", "slcp_distractors")

opts <- parse_args(list(
  tasks = TASK_ORDER,
  algorithms = c("npe", "nle"),
  budgets = c(1e3, 1e4, 1e5),
  observations = 1:10,
  seed = 1,
  max_epochs = 2000,
  dry_run = FALSE
))

tasks <- TASK_ORDER[TASK_ORDER %in% opts$tasks]
unknown <- setdiff(opts$tasks, TASK_ORDER)
if (length(unknown)) stop("Unknown task(s): ", paste(unknown, collapse = ", "))

log_dir <- file.path(.here, "results", "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

script <- file.path(.here, "01_run_benchmark.R")
failed <- character(0)

for (tk in tasks) {
  args <- c(script,
            "--tasks", tk,
            "--algorithms", paste(opts$algorithms, collapse = ","),
            "--budgets", paste(format(opts$budgets, scientific = FALSE,
                                      trim = TRUE), collapse = ","),
            "--observations", paste(opts$observations, collapse = ","),
            "--seed", opts$seed,
            "--max-epochs", opts$max_epochs,
            "--skip-existing")
  log_file <- file.path(log_dir, paste0(tk, ".log"))
  say("=== ", tk, "  (log: ", log_file, ")")
  if (isTRUE(opts$dry_run)) {
    say("  Rscript ", paste(args, collapse = " "))
    next
  }
  status <- system2("Rscript", args, stdout = log_file, stderr = log_file)
  if (status != 0) {
    say("  task exited with status ", status, "; see the log")
    failed <- c(failed, tk)
  }
}

if (length(failed)) {
  say("tasks that exited non-zero: ", paste(failed, collapse = ", "))
}
say("done. Now run: Rscript 03_report.R")

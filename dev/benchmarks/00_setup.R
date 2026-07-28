#!/usr/bin/env Rscript
# Fetch the two upstream repositories the benchmark scores against.
#
#   sbibm                  - task definitions, observations, reference posteriors
#   sbi-benchmark/results  - the published numbers from Lueckmann et al. (2021)
#
# Neither is used as Python. We read CSVs and a handful of pickled tensors.
# Both land in dev/benchmarks/external/, which is git-ignored; set SBIBM_PATH
# and SBIBM_RESULTS_PATH if you already have checkouts elsewhere.

.here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) normalizePath(dirname(sub("^--file=", "", a[1])))
  else normalizePath("dev/benchmarks")
})
for (f in c("utils.R", "pt_io.R", "sbibm_data.R")) source(file.path(.here, "R", f))

opts <- parse_args(list(force = FALSE))

REPOS <- list(
  sbibm = "https://github.com/sbi-benchmark/sbibm.git",
  results = "https://github.com/sbi-benchmark/results.git"
)

ext <- file.path(BENCH_DIR, "external")
dir.create(ext, showWarnings = FALSE, recursive = TRUE)

for (nm in names(REPOS)) {
  dest <- file.path(ext, nm)
  if (dir.exists(dest) && !opts$force) {
    say("already present: ", dest)
    next
  }
  if (dir.exists(dest)) unlink(dest, recursive = TRUE)
  say("cloning ", REPOS[[nm]])
  status <- system2("git", c("clone", "--depth", "1", REPOS[[nm]], shQuote(dest)))
  if (status != 0) {
    stop("git clone failed for ", nm, ". Clone it by hand and point ",
         toupper(nm), "_PATH at it.", call. = FALSE)
  }
}

say("checking that the files the benchmark needs are readable")
res <- paper_results()
say(sprintf("  main_paper.csv: %d rows, %d tasks, budgets %s",
            nrow(res), length(unique(res$task)),
            paste(sort(unique(res$num_simulations)), collapse = ", ")))

ref <- utils::read.csv(bzfile(file.path(sbibm_root(), "sbibm", "tasks",
                                        "two_moons", "files",
                                        "num_observation_1",
                                        "reference_posterior_samples.csv.bz2")))
say(sprintf("  two_moons reference posterior: %d x %d", nrow(ref), ncol(ref)))

say("checking neuralsbi and torch")
load_neuralsbi()
say("  neuralsbi ", as.character(utils::packageVersion("neuralsbi")))
if (!requireNamespace("torch", quietly = TRUE) || !torch::torch_is_installed()) {
  say("  torch/libtorch NOT available. The benchmark needs it: the estimators ",
      "the paper used (nsf for NPE, maf for NLE) are torch-backed.")
  say("  install.packages('torch'); torch::install_torch()")
} else {
  say("  torch ", as.character(utils::packageVersion("torch")), " ready")
}

say("setup complete")

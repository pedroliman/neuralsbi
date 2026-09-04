#!/usr/bin/env Rscript
# Compare what neuralsbi produced against the numbers in Lueckmann et al. (2021)
# and say, per benchmark, whether we reproduce the paper.
#
# The comparison is per (task, algorithm, simulation budget). For each cell we
# take the median C2ST over the observations we actually ran, and the paper's
# median over *those same observations*, so a partial run is still a fair
# comparison. C2ST is an accuracy: 0.5 means the approximate posterior is
# indistinguishable from the reference, 1.0 means perfectly separable, so lower
# is better and "we beat the paper" is a pass.
#
#   Rscript 03_report.R
#   Rscript 03_report.R --tolerance 0.03 --detail
#
# Writes results/comparison.csv and results/report.md.

.here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) normalizePath(dirname(sub("^--file=", "", a[1])))
  else normalizePath("dev/benchmarks")
})
for (f in c("utils.R", "pt_io.R", "sbibm_data.R", "runner.R")) {
  source(file.path(.here, "R", f))
}

opts <- parse_args(list(
  tolerance = 0.05,
  min_observations = 1,
  detail = FALSE,
  out = ""
))

if (!file.exists(metrics_path())) {
  stop("No results at ", metrics_path(), ". Run 01_run_benchmark.R first.",
       call. = FALSE)
}

ours <- utils::read.csv(metrics_path(), stringsAsFactors = FALSE)
ours <- ours[!is.na(ours$C2ST), , drop = FALSE]
paper <- paper_results()

cells <- unique(ours[, c("task", "algorithm", "num_simulations")])
cells <- cells[order(cells$task, cells$algorithm, cells$num_simulations), ]

rows <- lapply(seq_len(nrow(cells)), function(i) {
  ce <- cells[i, ]
  mine <- ours[ours$task == ce$task & ours$algorithm == ce$algorithm &
                 ours$num_simulations == ce$num_simulations, ]
  ref <- paper[paper$task == ce$task & paper$algorithm == ce$algorithm &
                 paper$num_simulations == ce$num_simulations, ]
  # Match on observation so a partial run compares like with like.
  ref <- ref[ref$num_observation %in% mine$num_observation, ]
  mine <- mine[mine$num_observation %in% ref$num_observation, ]
  if (!nrow(mine)) return(NULL)
  data.frame(
    task = ce$task,
    algorithm = ce$algorithm,
    num_simulations = ce$num_simulations,
    n_observations = nrow(mine),
    c2st_neuralsbi = stats::median(mine$C2ST),
    c2st_paper = stats::median(ref$C2ST),
    delta = stats::median(mine$C2ST) - stats::median(ref$C2ST),
    c2st_neuralsbi_min = min(mine$C2ST),
    c2st_neuralsbi_max = max(mine$C2ST),
    minutes = sum(mine$time_total) / 60,
    stringsAsFactors = FALSE
  )
})
cmp <- do.call(rbind, rows)
cmp <- cmp[cmp$n_observations >= opts$min_observations, , drop = FALSE]

tol <- opts$tolerance
cmp$verdict <- ifelse(cmp$delta < -tol, "BETTER",
               ifelse(cmp$delta <= tol, "MATCH", "WORSE"))

# --- console output ----------------------------------------------------------

fmt <- function(d) {
  data.frame(
    task = d$task,
    alg = d$algorithm,
    sims = vapply(d$num_simulations, budget_label, character(1)),
    obs = d$n_observations,
    ours = sprintf("%.3f", d$c2st_neuralsbi),
    range = sprintf("%.2f-%.2f", d$c2st_neuralsbi_min, d$c2st_neuralsbi_max),
    paper = sprintf("%.3f", d$c2st_paper),
    delta = sprintf("%+.3f", d$delta),
    verdict = d$verdict,
    stringsAsFactors = FALSE
  )
}

cat("\nC2ST against sbibm reference posteriors (median over observations).\n")
cat("Lower is better. A cell passes when it is no more than ", tol,
    " above the paper.\n\n", sep = "")
print(fmt(cmp), row.names = FALSE)

n_ok <- sum(cmp$verdict %in% c("MATCH", "BETTER"))
cat(sprintf("\n%d / %d cells reproduce the paper (tolerance %.3f).\n",
            n_ok, nrow(cmp), tol))
if (n_ok < nrow(cmp)) {
  bad <- cmp[cmp$verdict == "WORSE", ]
  cat("\nCells where neuralsbi is worse than the paper:\n")
  print(fmt(bad), row.names = FALSE)
}

missing <- setdiff(
  paste(paper$task, paper$algorithm, paper$num_simulations, sep = "|"),
  paste(cmp$task, cmp$algorithm, cmp$num_simulations, sep = "|")
)
missing <- unique(missing[grepl("\\|(NPE|NLE)\\|", missing)])
if (length(missing)) {
  cat(sprintf("\n%d NPE/NLE cell(s) of the paper's grid have no results yet.\n",
              length(missing)))
}

if (isTRUE(opts$detail)) {
  cat("\nPer-observation detail:\n")
  det <- merge(
    ours[, c("task", "algorithm", "num_simulations", "num_observation", "C2ST")],
    paper[, c("task", "algorithm", "num_simulations", "num_observation", "C2ST")],
    by = c("task", "algorithm", "num_simulations", "num_observation"),
    suffixes = c("_neuralsbi", "_paper")
  )
  det <- det[order(det$task, det$algorithm, det$num_simulations,
                   det$num_observation), ]
  det$delta <- det$C2ST_neuralsbi - det$C2ST_paper
  print(det, row.names = FALSE, digits = 3)
}

# --- files -------------------------------------------------------------------

out_csv <- file.path(results_dir(), "comparison.csv")
utils::write.csv(cmp, out_csv, row.names = FALSE)

md_path <- if (nzchar(opts$out)) opts$out else file.path(results_dir(), "report.md")
md <- c(
  "# neuralsbi vs. Benchmarking Simulation-Based Inference",
  "",
  sprintf("Generated %s with neuralsbi %s.", format(Sys.time(), "%Y-%m-%d %H:%M"),
          as.character(utils::packageVersion("neuralsbi"))),
  "",
  paste0("Metric: C2ST between the approximate posterior and sbibm's reference ",
         "posterior, z-scored, 10k samples on each side, median over the ",
         "observations run. Reference column is the published median over the ",
         "same observations. A cell passes when it sits no more than ",
         sprintf("%.3f", tol), " above the paper."),
  "",
  sprintf("%d of %d cells reproduce the paper.", n_ok, nrow(cmp)),
  "",
  "| task | algorithm | simulations | obs | neuralsbi | range | paper | delta | verdict |",
  "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
  apply(fmt(cmp), 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
)
writeLines(md, md_path)

cat("\nwrote ", out_csv, "\n", sep = "")
cat("wrote ", md_path, "\n", sep = "")

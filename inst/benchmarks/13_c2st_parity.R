#!/usr/bin/env Rscript
# Generate the sample-set pairs the C2ST parity check runs on, and score them
# with neuralsbi's c2st().
#
# c2st() claims to run the same test as sbibm/metrics/c2st.py. This script and
# its Python twin (14_c2st_parity.py) check that claim on pairs whose answer is
# known by construction: two draws from the same distribution should score near
# 0.5, and a shift, a change of scale or a higher-dimensional change of scale
# should score above it by the same amount on both sides. The two fits use
# different random number streams, so agreement means agreement to Monte-Carlo
# noise, not to the digit.
#
# Usage:
#   Rscript 13_c2st_parity.R
#   python  14_c2st_parity.py
suppressMessages(library(neuralsbi))

out_dir <- file.path("results", "c2st_parity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)
cases <- list(
  same_d2  = list(x = matrix(rnorm(2000 * 2), ncol = 2),
                  y = matrix(rnorm(2000 * 2), ncol = 2)),
  shift_d2 = list(x = matrix(rnorm(2000 * 2), ncol = 2),
                  y = matrix(rnorm(2000 * 2, mean = 0.5), ncol = 2)),
  scale_d2 = list(x = matrix(rnorm(2000 * 2), ncol = 2),
                  y = matrix(rnorm(2000 * 2, sd = 2), ncol = 2)),
  shift_d1 = list(x = matrix(rnorm(2000), ncol = 1),
                  y = matrix(rnorm(2000, mean = 0.3), ncol = 1)),
  scale_d5 = list(x = matrix(rnorm(1500 * 5), ncol = 5),
                  y = matrix(rnorm(1500 * 5, sd = 1.3), ncol = 5))
)

rows <- list()
for (nm in names(cases)) {
  x <- cases[[nm]]$x
  y <- cases[[nm]]$y
  write.csv(x, file.path(out_dir, paste0(nm, "_x.csv")), row.names = FALSE)
  write.csv(y, file.path(out_dir, paste0(nm, "_y.csv")), row.names = FALSE)
  res <- c2st(x, y, seed = 1)
  rows[[length(rows) + 1L]] <- data.frame(case = nm, accuracy = res$accuracy,
                                          auc = res$auc)
}
tab <- do.call(rbind, rows)
print(tab, digits = 4, row.names = FALSE)
write.csv(tab, file.path(out_dir, "neuralsbi.csv"), row.names = FALSE)
cat("wrote", file.path(out_dir, "neuralsbi.csv"), "\n")
cat("now run: python 14_c2st_parity.py\n")

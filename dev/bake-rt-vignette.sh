#!/usr/bin/env bash
# Rebuild vignette("sir-time-varying-beta") from its .Rmd.orig source.
#
#     bash dev/bake-rt-vignette.sh
#
# Installs whatever Suggests are missing, installs the package, then bakes the
# one vignette (not all six). Needs internet: the article pulls NYT case counts
# and Census population estimates at knit time.
#
# Cost: three amortized NPE fits. On 4 CPU cores that is roughly 25-40 minutes
# each, so budget about two hours. Set NSBI_SMOKE=1 for a ~10-minute pass that
# only checks the article runs end to end -- its numbers are not publishable.
set -euo pipefail
cd "$(dirname "$0")/.."

echo "==> checking R packages"
Rscript -e '
need <- setdiff(c("torch", "ggplot2", "GGally", "knitr", "rmarkdown"),
                rownames(installed.packages()))
if (length(need)) {
  message("installing: ", paste(need, collapse = ", "))
  install.packages(need, repos = "https://cloud.r-project.org")
}
if (!torch::torch_is_installed()) {
  message("installing libtorch (one-off, ~100 MB)")
  torch::install_torch()
}
cat("torch ready:", torch::torch_is_installed(), "\n")'

echo "==> installing neuralsbi"
R CMD INSTALL --no-docs . >/dev/null

if [ "${NSBI_SMOKE:-0}" = "1" ]; then
  echo "==> SMOKE MODE: cutting the training budget, results are throwaway"
  cp vignettes/sir-time-varying-beta.Rmd.orig /tmp/rt-vignette-backup.Rmd.orig
  trap 'mv /tmp/rt-vignette-backup.Rmd.orig vignettes/sir-time-varying-beta.Rmd.orig' EXIT
  sed -i.bak 's/^n_sim <- 15000$/n_sim <- 1500/; s/max_epochs = 300/max_epochs = 30/; s/n_sbc = 150/n_sbc = 30/' \
    vignettes/sir-time-varying-beta.Rmd.orig
  rm -f vignettes/sir-time-varying-beta.Rmd.orig.bak
fi

echo "==> baking (this is the slow part)"
time Rscript vignettes/precompute.R sir-time-varying

echo
echo "==> done. Regenerated:"
echo "    vignettes/sir-time-varying-beta.Rmd"
echo "    vignettes/figures/sir-rt-*.png"
echo
echo "Preview with:"
echo "    Rscript -e 'rmarkdown::render(\"vignettes/sir-time-varying-beta.Rmd\")'"

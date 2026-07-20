#!/bin/bash
# Install the R toolchain this package needs for development: R itself,
# testthat for the test suite, and knitr/rmarkdown/pkgdown for vignettes and
# the website. The container state is cached after the hook completes, so the
# expensive first run pays off across sessions.
set -euo pipefail

# Only needed on Claude Code on the web; local machines manage their own R.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

export DEBIAN_FRONTEND=noninteractive

# System layer: R, pandoc, and the headers the R packages compile against
# (libuv for fs/testthat; fontconfig/freetype/png/tiff/jpeg/harfbuzz/fribidi
# for the pkgdown graphics stack; xml2/ssl/curl for its web stack).
if ! command -v R > /dev/null 2>&1; then
  apt-get update -qq
  apt-get install -y -qq \
    r-base-core pandoc \
    libuv1-dev libfontconfig1-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev libharfbuzz-dev libfribidi-dev \
    libxml2-dev libssl-dev libcurl4-openssl-dev
fi

# R package layer. Install only what is missing so warm starts are fast.
Rscript --no-init-file -e '
  pkgs <- c("testthat", "knitr", "rmarkdown", "pkgdown")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    options(Ncpus = max(1L, parallel::detectCores() - 1L))
    install.packages(missing, repos = "https://cloud.r-project.org")
    still <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
    if (length(still)) stop("failed to install: ", paste(still, collapse = ", "))
  }
'

echo "R toolchain ready: $(R --version | head -1)"

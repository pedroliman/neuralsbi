# Precompute (cache) the package vignettes.
#
# The neural estimators shown in these vignettes need `torch`/libtorch and take
# too long to train on CI (which usually has neither libtorch nor a GPU). So we
# evaluate the expensive `*.Rmd.orig` sources *here*, once, and commit the
# resulting `*.Rmd` -- with output and figures already baked in -- together with
# the figures under `vignettes/figures/`.
#
# On CI, in `R CMD check`, and in the pkgdown build, the committed `.Rmd` files
# are plain Markdown: every chunk has already been turned into a static code
# block, so nothing is re-evaluated, no torch is required, and the build stays
# fast.
#
# Re-run this whenever you edit a `*.Rmd.orig`:
#
#     R CMD INSTALL --no-docs .        # make the current sources importable
#     Rscript vignettes/precompute.R   # bake the vignettes
#
# It needs the package installed (so `library(neuralsbi)` resolves) and torch
# available (`torch::install_torch()`). Install `future` too: each vignette
# declares `plan(multisession)`, which is the difference between minutes and
# hours on the two SIR articles. Commit the regenerated `*.Rmd` and
# `vignettes/figures/` alongside your source change.

if (!requireNamespace("knitr", quietly = TRUE)) {
  stop("knitr is required to precompute the vignettes.")
}
if (!requireNamespace("svglite", quietly = TRUE)) {
  stop("svglite is required to bake the vignettes' figures as SVG.")
}
# Figures are SVG by default, except for the handful of panels that draw one
# mark per posterior draw. Those chunks set dev = "ragg_png", because a scatter
# of several thousand points is many times larger as SVG than as PNG.
if (!requireNamespace("ragg", quietly = TRUE)) {
  stop("ragg is required: the dense scatter panels bake to PNG with it.")
}
if (!requireNamespace("torch", quietly = TRUE) || !torch::torch_is_installed()) {
  stop("torch (libtorch) is required to bake the neural vignettes. ",
       "Install with install.packages('torch'); torch::install_torch().")
}

# Every vignette declares its own `future` plan in its setup chunk, so the
# simulations run across cores. `future` is a Suggests; without it the bake
# still works, one core at a time, and the SIR articles take hours instead of
# minutes. Warn rather than stop, so a quick bake of a cheap article is not
# blocked on installing it.
if (!requireNamespace("future", quietly = TRUE)) {
  warning("future is not installed, so the vignettes will bake sequentially. ",
          "install.packages('future') to use all your cores.", call. = FALSE)
}

# Show the simulation and training bars while baking. They are off under knitr
# by default, which is right for the articles but unhelpful when you are
# watching a run that takes hours. The bar writes to stderr, so it reaches the
# console without landing in the baked .Rmd. Each bake runs in its own R
# process, so the option is set again there rather than inherited.

# knit each source with the working directory inside vignettes/, so the baked
# figure paths ("figures/<name>-1.svg") are relative to the vignette itself and
# resolve the same way under R CMD build and pkgdown.
vign_dir <- if (basename(getwd()) == "vignettes") getwd() else file.path(getwd(), "vignettes")
old_wd <- setwd(vign_dir)
on.exit(setwd(old_wd), add = TRUE)

origs <- sort(list.files(".", pattern = "\\.Rmd\\.orig$"))
if (length(origs) == 0) stop("No *.Rmd.orig sources found in vignettes/.")

# Optional: bake only the vignettes named on the command line, matched as
# substrings of the file name. Baking one article beats baking all six when
# you are iterating on a single source.
#     Rscript vignettes/precompute.R sir-time-varying
selected <- commandArgs(trailingOnly = TRUE)
if (length(selected)) {
  keep <- Reduce(`|`, lapply(selected, function(p) grepl(p, origs, fixed = TRUE)))
  if (!any(keep)) {
    stop("No *.Rmd.orig matched ", paste(selected, collapse = ", "), ". Available: ",
         paste(origs, collapse = ", "), call. = FALSE)
  }
  origs <- origs[keep]
  message("Baking only: ", paste(origs, collapse = ", "))
}

rscript <- file.path(R.home("bin"), "Rscript")

for (orig in origs) {
  out <- sub("\\.orig$", "", orig)
  message("Baking ", orig, " -> ", out)

  # One R process per article, rather than one knit() call after another in
  # this session. A fresh `envir` isolates the objects a chunk creates but not
  # the packages it attaches, and an article whose `library(neuralsbi)` runs
  # second prints no startup message, so its baked text would lose the note
  # that `sample()` masks the base function. Baking every article has to give
  # the same file as baking that one article alone.
  code <- sprintf(
    'options(neuralsbi.progress = "builtin"); knitr::knit(%s, output = %s, envir = new.env(parent = globalenv()), quiet = FALSE)',
    encodeString(orig, quote = '"'), encodeString(out, quote = '"'))
  status <- system2(rscript, c("-e", shQuote(code)))
  if (status != 0) {
    stop(sprintf("Baking %s failed (Rscript exited %d).", orig, status),
         call. = FALSE)
  }

  # The chunks set error = FALSE, so a genuine failure aborts the knit above.
  # As a backstop, refuse to leave a baked vignette that carries an error
  # trace (e.g. a torch chunk that silently degraded) -- better to fail here
  # than to commit a broken article.
  baked <- readLines(out, warn = FALSE)
  bad <- grep("^#> Error|libtorch is not installed|not found$", baked)
  if (length(bad)) {
    stop(sprintf("%s still contains error output at line(s) %s. Fix the ",
                 out, paste(head(bad, 5), collapse = ", ")),
         "environment (torch installed? package up to date?) and re-run.",
         call. = FALSE)
  }
}

message("Done. Review the regenerated *.Rmd and figures/, then commit them.")

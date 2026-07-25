## Fast standalone comparison of the three epidemic models in
## vignette("sir-time-varying-beta"), without baking the whole article.
##
##     Rscript dev/compare-behavioral-models.R [n_simulations] [max_epochs] [model]
##
## Defaults to 15000 / 300 / all, which is what the vignette uses. For a quick
## look, try:
##
##     Rscript dev/compare-behavioral-models.R 2000 40 EFB
##
## It reuses the vignette source as the single definition of the models: the
## chunks are extracted and evaluated here, so this script and the article can
## never disagree about what the simulators are.

args    <- commandArgs(trailingOnly = TRUE)
n_sim   <- if (length(args) >= 1) as.integer(args[1]) else 15000L
max_ep  <- if (length(args) >= 2) as.integer(args[2]) else 300L
which_m <- if (length(args) >= 3) args[3] else "all"

src <- "vignettes/sir-time-varying-beta.Rmd.orig"
if (!file.exists(src)) {
  stop("run this from the package root: Rscript dev/compare-behavioral-models.R",
       call. = FALSE)
}

## Evaluate the article's setup chunks, stopping before it starts training.
lines <- readLines(src, warn = FALSE)
starts <- grep("^```\\{r", lines)
ends   <- grep("^```$", lines)
stop_at <- grep("^```\\{r fit-npe", lines)
if (!length(stop_at)) stop("could not find the fit-npe chunk", call. = FALSE)

env <- new.env(parent = globalenv())
for (s in starts[starts < stop_at]) {
  e <- ends[ends > s][1]
  code <- lines[(s + 1):(e - 1)]
  eval(parse(text = paste(code, collapse = "\n")), envir = env)
}
message("model definitions loaded from ", basename(src))

library(neuralsbi)
sims <- get("sims", envir = env)
if (which_m != "all") sims <- sims[which_m]
smooth7      <- get("smooth7", envir = env)
case_daily_s <- get("case_daily_s", envir = env)
pop          <- get("pop", envir = env)
days         <- get("days", envir = env)
gamma_fixed  <- get("gamma_fixed", envir = env)
embed        <- get("embed", envir = env)

obs_for <- function(s) c(log(pop[[s]]), log1p(case_daily_s[s, ]))
highlight <- c("Washington", "New York", "California", "Louisiana",
               "Michigan", "West Virginia")

for (nm in names(sims)) {
  m <- sims[[nm]]
  cat(sprintf("\n########## %s ##########\n", nm))
  t0 <- Sys.time()
  fit <- npe(m$prior, m$sim, n_simulations = n_sim, density_estimator = "maf",
             embedding_net = embed, max_epochs = max_ep, n_restarts = 1,
             seed = 1, verbose = FALSE)
  cat(sprintf("fit time: %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  print(fit)
  tot <- 0
  for (s in highlight) {
    draws <- sample(posterior(fit, x_obs = obs_for(s)), 400)
    r <- m$sim(draws, N_fixed = pop[[s]], full = TRUE)
    pred <- smooth7(r$reported)
    obs <- case_daily_s[s, ]; j <- which.max(obs)
    med <- apply(pred, 2, stats::median)
    rmse <- sqrt(mean((log1p(med) - log1p(obs))^2)); tot <- tot + rmse
    cat(sprintf("  %-14s rmse %.3f | peak %8.0f vs %8.0f (%.2fx) | 90%% %7.0f-%8.0f | R0 %.2f | attack %.3f\n",
                s, rmse, med[j], obs[j], med[j] / obs[j],
                stats::quantile(pred[, j], .05), stats::quantile(pred[, j], .95),
                stats::median(r$beta_eff[, 1]) / gamma_fixed,
                stats::median(1 - r$S_frac[, days])))
  }
  cat(sprintf("  == %s mean rmse %.3f over %d states ==\n",
              nm, tot / length(highlight), length(highlight)))
}

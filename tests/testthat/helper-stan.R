# Shared across test files: skip the Stan round-trip tests unless a working
# CmdStan is present. Same contract as skip_if_no_torch() -- the suite runs
# everywhere, and the tests that need a toolchain say so instead of failing.
skip_if_no_cmdstan <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  ok <- tryCatch(!is.null(cmdstanr::cmdstan_version()), error = function(e) FALSE)
  testthat::skip_if_not(isTRUE(ok), "CmdStan not installed")
}

# posterior(sampler = "stan") prefers cmdstanr and falls back to rstan, so on a
# machine with cmdstanr the fallback is never reached by the dispatch. The
# rstan test therefore calls stan_run_rstan() directly.
#
# Installing rstan is not enough to run it. rstan compiles the model at call
# time and needs the headers it links against plus a C++ compiler, none of
# which come with a binary install -- CI hit this as "Boost not found; call
# install.packages('BH')" from inside stan_model(), which is a missing
# toolchain reported as a test failure. Check for the toolchain up front so it
# skips like every other optional dependency.
skip_if_no_rstan <- function() {
  testthat::skip_if_not_installed("rstan")
  for (pkg in c("BH", "RcppEigen", "StanHeaders")) {
    testthat::skip_if_not_installed(pkg)
  }
  testthat::skip_if_not(
    nzchar(Sys.which("g++")) || nzchar(Sys.which("clang++")),
    "no C++ compiler for rstan to build against"
  )
}

# Compile the generated `functions` block on its own and evaluate it at fixed
# (theta, x). This is the only way to check the transpiled code actually means
# what log_lik() means: a whole model would fold in the prior and the
# constraint Jacobians.
#
# sig_figs = 18 matters. CmdStan writes 6 significant digits by default, which
# is coarse enough to look like a code-generation bug. The sampler itself is
# vestigial -- one iteration over a dummy parameter -- so its diagnostics are
# switched off rather than warned about.
stan_eval_log_lik <- function(fit, theta, x_obs) {
  code <- paste0(
    stan_code(fit, model = FALSE),
    sprintf('
data {
  int<lower=1> M;
  int<lower=1> N;
  array[M] vector[%d] th;
  matrix[N, %d] x;
  int<lower=1> nsbi_nw;
  vector[nsbi_nw] nsbi_w;
}
parameters { real dummy; }
model { dummy ~ std_normal(); }
generated quantities {
  vector[M] lp_sum;
  vector[M] lp_one;
  for (m in 1:M) {
    lp_sum[m] = nsbi_log_lik_sum_lpdf(x | th[m], nsbi_w);
    lp_one[m] = nsbi_log_lik_lpdf(x[1]\' | th[m], nsbi_w);
  }
}
', fit$dim_theta, fit$dim_x))

  file <- tempfile(fileext = ".stan")
  on.exit(unlink(file), add = TRUE)
  writeLines(code, file)

  model <- cmdstanr::cmdstan_model(file, quiet = TRUE)
  data <- stan_data(fit, x_obs)
  data$M <- nrow(theta)
  data$th <- theta
  run <- model$sample(data = data, chains = 1, iter_warmup = 1,
                      iter_sampling = 1, refresh = 0, show_messages = FALSE,
                      seed = 1, sig_figs = 18, diagnostics = NULL)
  list(sum = as.numeric(run$draws("lp_sum")),
       one = as.numeric(run$draws("lp_one")))
}

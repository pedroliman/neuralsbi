#' Visualize posterior samples
#'
#' A pair plot built on [GGally::ggpairs()]: 1-D marginal densities on the
#' diagonal and 2-D scatter in the lower triangle, with optional markers for a
#' reference (e.g. true) parameter value. Analogous to Python `sbi`'s
#' `pairplot`.
#'
#' @param samples A matrix of posterior draws (rows = draws), or an
#'   `nsbi_samples` object.
#' @param truth Optional reference parameter vector to overlay.
#' @param labels Optional parameter labels.
#' @param limits Optional list (one `c(lo, hi)` per parameter, in column
#'   order) or matrix of per-parameter axis limits.
#' @param col Point and density fill colour.
#' @param alpha Point and density fill transparency.
#' @param ... Passed to the lower-triangle [ggplot2::geom_point()] layer.
#' @return A `ggmatrix` object (also drawn as a side effect), invisibly.
#' @export
pairplot <- function(samples, truth = NULL, labels = NULL, limits = NULL,
                     col = "steelblue", alpha = 0.4, ...) {
  require_ggplot2(ggally = TRUE)
  X <- as_theta_matrix(samples)
  d <- ncol(X)
  if (is.null(labels)) {
    labels <- colnames(X) %||% paste0("theta[", seq_len(d), "]")
  }
  colnames(X) <- labels
  df <- as.data.frame(X, check.names = FALSE)

  lims <- NULL
  if (!is.null(limits)) {
    lims <- stats::setNames(
      lapply(seq_len(d), function(j) if (is.list(limits)) limits[[j]] else limits[j, ]),
      labels
    )
  }
  truth_df <- NULL
  if (!is.null(truth)) {
    truth_df <- as.data.frame(as.list(stats::setNames(as.numeric(truth), labels)),
                              check.names = FALSE)
  }

  lower_fn <- function(data, mapping, ...) {
    xn <- rlang::as_name(mapping$x)
    yn <- rlang::as_name(mapping$y)
    p <- ggplot2::ggplot(data, mapping) +
      ggplot2::geom_point(colour = col, alpha = alpha, size = 0.6, ...)
    if (!is.null(truth_df)) {
      p <- p +
        ggplot2::geom_vline(xintercept = truth_df[[xn]], colour = "firebrick", linewidth = 0.6) +
        ggplot2::geom_hline(yintercept = truth_df[[yn]], colour = "firebrick", linewidth = 0.6) +
        ggplot2::geom_point(data = truth_df, mapping = ggplot2::aes(x = .data[[xn]], y = .data[[yn]]),
                            colour = "firebrick", shape = 4, size = 2.5, stroke = 1,
                            inherit.aes = FALSE)
    }
    if (!is.null(lims)) p <- p + ggplot2::coord_cartesian(xlim = lims[[xn]], ylim = lims[[yn]])
    p
  }
  diag_fn <- function(data, mapping, ...) {
    xn <- rlang::as_name(mapping$x)
    p <- ggplot2::ggplot(data, mapping) +
      ggplot2::geom_density(colour = col, fill = col, alpha = alpha)
    if (!is.null(truth_df)) {
      p <- p + ggplot2::geom_vline(xintercept = truth_df[[xn]], colour = "firebrick", linewidth = 0.6)
    }
    if (!is.null(lims)) p <- p + ggplot2::coord_cartesian(xlim = lims[[xn]])
    p
  }

  p <- GGally::ggpairs(df,
    lower = list(continuous = lower_fn),
    diag  = list(continuous = diag_fn),
    upper = list(continuous = "blank"),
    progress = FALSE
  ) + ggplot2::theme_minimal()

  print(p)
  invisible(p)
}

#' Plot an SBC rank histogram
#'
#' Uniform bars indicate calibration; a U shape means the posterior is too
#' narrow (overconfident); an inverted-U means it is too wide.
#'
#' @param sbc_result An `nsbi_sbc` object from [sbc()].
#' @param param Which parameter index to plot (default 1).
#' @param bins Number of histogram bins.
#' @return A `ggplot` object (also drawn as a side effect), invisibly.
#' @export
plot_sbc <- function(sbc_result, param = 1L, bins = 20L) {
  stopifnot(inherits(sbc_result, "nsbi_sbc"))
  require_ggplot2()
  r <- sbc_result$ranks[, param]
  L <- sbc_result$n_posterior_samples
  breaks <- seq(0, L, length.out = bins + 1L)
  expected <- sbc_result$n_sbc / bins
  ci <- stats::qbinom(c(0.005, 0.995), sbc_result$n_sbc, 1 / bins)

  p <- ggplot2::ggplot(data.frame(rank = r), ggplot2::aes(x = .data$rank)) +
    ggplot2::geom_histogram(breaks = breaks, fill = "grey80", colour = "white") +
    ggplot2::geom_hline(yintercept = expected, colour = "firebrick", linewidth = 0.7, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = ci, colour = "firebrick", linewidth = 0.5, linetype = "dotted") +
    ggplot2::labs(title = sprintf("SBC ranks: parameter %d", param),
                 x = "rank of true value", y = "count") +
    ggplot2::theme_minimal()

  print(p)
  invisible(r)
}

#' Plot nominal vs. empirical credible-interval coverage
#'
#' Well-calibrated posteriors lie on the diagonal. Curves above the diagonal
#' mean the posterior is too wide (conservative); below means overconfident.
#' A shaded band shows the Monte-Carlo uncertainty from the finite number of
#' SBC trials.
#'
#' @param sbc_result An `nsbi_sbc` object from [sbc()].
#' @param levels Nominal credibility levels to evaluate.
#' @return A `ggplot` object (also drawn as a side effect), invisibly.
#' @export
plot_coverage <- function(sbc_result, levels = seq(0.05, 0.95, by = 0.05)) {
  stopifnot(inherits(sbc_result, "nsbi_sbc"))
  require_ggplot2()
  cov <- expected_coverage(sbc_result, levels = levels)
  d <- ncol(cov) - 1L
  n <- sbc_result$n_sbc
  band <- data.frame(
    nominal = cov$nominal,
    lo = stats::qbinom(0.005, n, cov$nominal) / n,
    hi = stats::qbinom(0.995, n, cov$nominal) / n
  )
  long <- do.call(rbind, lapply(seq_len(d), function(j) {
    data.frame(nominal = cov$nominal, empirical = cov[[j + 1L]],
              parameter = colnames(cov)[j + 1L])
  }))

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = band, ggplot2::aes(x = .data$nominal, ymin = .data$lo, ymax = .data$hi),
                         fill = "grey60", alpha = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = "dashed") +
    ggplot2::geom_line(data = long,
                       ggplot2::aes(x = .data$nominal, y = .data$empirical, colour = .data$parameter),
                       linewidth = 0.8) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(title = "Expected coverage", x = "nominal credibility level",
                 y = "empirical coverage", colour = NULL) +
    ggplot2::theme_minimal()

  print(p)
  invisible(cov)
}

#' Plot TARP expected coverage
#'
#' Draws the expected coverage probability (ECP) curve from [tarp()] against
#' the nominal credibility level. A calibrated posterior lies on the diagonal;
#' a curve above the diagonal means the posterior is too wide (conservative),
#' below means overconfident. The shaded band shows the Monte-Carlo uncertainty
#' from the finite number of TARP trials.
#'
#' @param tarp_result An `nsbi_tarp` object from [tarp()].
#' @return A `ggplot` object (also drawn as a side effect), invisibly.
#' @export
plot_tarp <- function(tarp_result) {
  stopifnot(inherits(tarp_result, "nsbi_tarp"))
  require_ggplot2()
  lev <- tarp_result$levels
  n <- tarp_result$n_tarp
  band <- data.frame(
    nominal = lev,
    lo = stats::qbinom(0.005, n, lev) / n,
    hi = stats::qbinom(0.995, n, lev) / n
  )
  line <- data.frame(nominal = lev, ecp = tarp_result$ecp)

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = band, ggplot2::aes(x = .data$nominal, ymin = .data$lo, ymax = .data$hi),
                         fill = "grey60", alpha = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = "dashed") +
    ggplot2::geom_line(data = line, ggplot2::aes(x = .data$nominal, y = .data$ecp),
                       colour = "steelblue", linewidth = 0.8) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(title = "TARP coverage", x = "nominal credibility level",
                 y = "expected coverage probability") +
    ggplot2::theme_minimal()

  print(p)
  invisible(line)
}

#' Plot posterior predictive checks
#'
#' Compares data simulated from posterior parameter draws (see
#' [posterior_predictive()]) with the observed data, one marginal histogram per
#' data dimension with the observation marked. If the observation falls in the
#' tails of the predictive distribution, the model (or the fit) does not
#' reproduce the data it is conditioned on.
#'
#' @param pred A matrix of predictive draws from [posterior_predictive()].
#' @param x_obs The observed data vector the posterior was conditioned on.
#' @param labels Optional labels for the data dimensions.
#' @param bins Number of histogram bins.
#' @return A `ggplot` object (also drawn as a side effect), invisibly.
#' @export
plot_posterior_predictive <- function(pred, x_obs, labels = NULL, bins = 30L) {
  require_ggplot2()
  pred <- as_theta_matrix(pred)
  d <- ncol(pred)
  x_obs <- as.numeric(x_obs)
  stopifnot(length(x_obs) == d)
  if (is.null(labels)) {
    labels <- colnames(pred) %||% paste0("x[", seq_len(d), "]")
  }
  q <- stats::setNames(vapply(seq_len(d), function(j) mean(pred[, j] < x_obs[j]), numeric(1)), labels)

  long <- data.frame(
    value = as.vector(pred),
    dim = factor(rep(labels, each = nrow(pred)), levels = labels)
  )
  obs_df <- data.frame(dim = factor(labels, levels = labels), obs = x_obs)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = bins, fill = "grey80", colour = "white") +
    ggplot2::geom_vline(data = obs_df, ggplot2::aes(xintercept = .data$obs),
                        colour = "firebrick", linewidth = 0.8) +
    ggplot2::facet_wrap(~dim, scales = "free") +
    ggplot2::labs(x = "predictive draws", y = "count") +
    ggplot2::theme_minimal()

  print(p)
  invisible(q)
}

#' @export
print.nsbi_samples <- function(x, ...) {
  acc <- attr(x, "acceptance_rate")
  cat(sprintf("<nsbi_samples> %d draws x %d parameters\n", nrow(x), ncol(x)))
  if (!is.null(acc)) cat(sprintf("  support acceptance rate: %.3f\n", acc))
  cls <- setdiff(class(x), "nsbi_samples")
  print(utils::head(`class<-`(unclass(x), cls)))
  invisible(x)
}

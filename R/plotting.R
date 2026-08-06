#' Visualize posterior samples
#'
#' A pair plot built on [GGally::ggpairs()]: 1-D marginal densities on the
#' diagonal and 2-D highest-density regions (via [ggdensity::geom_hdr()]) in
#' the lower triangle, with optional markers for a reference (e.g. true)
#' parameter value. Analogous to Python `sbi`'s `pairplot`.
#'
#' @param samples A matrix of posterior draws (rows = draws), or an
#'   `nsbi_samples` object.
#' @param truth Optional reference parameter vector to overlay.
#' @param labels Optional parameter labels. Defaults to `colnames(samples)` (set
#'   automatically for [sample()] draws from a fit with named parameters --
#'   see [prior_uniform()]/[prior_normal()]) or `theta[1]`, `theta[2]`, ....
#'   Labels that parse as R syntax (`"beta[1]"`, `"rho"`) render as their
#'   plotmath symbol.
#' @param limits Optional list (one `c(lo, hi)` per parameter, in column
#'   order) or matrix of per-parameter axis limits.
#' @param col Density-region and marginal-density fill colour.
#' @param alpha Marginal-density fill transparency. The lower-triangle
#'   highest-density regions shade themselves by probability level (99/95/80/50%)
#'   instead, via `ggdensity`'s own `alpha` mapping.
#' @param ... Passed to the lower-triangle [ggdensity::geom_hdr()] layer.
#' @return A `ggmatrix` object (also drawn as a side effect), invisibly.
#' @export
pairplot <- function(samples, truth = NULL, labels = NULL, limits = NULL,
                     col = "steelblue", alpha = 0.4, ...) {
  require_ggplot2(ggally = TRUE, ggdensity = TRUE)
  X <- as_theta_matrix(samples)
  d <- ncol(X)
  if (is.null(labels)) {
    labels <- colnames(X) %||% paste0("theta[", seq_len(d), "]")
  }
  labels <- math_safe_text(labels)
  colnames(X) <- labels
  df <- as.data.frame(X, check.names = FALSE)

  lims <- NULL
  if (!is.null(limits)) {
    # One c(lo, hi) per parameter, in column order. A list or matrix of the
    # wrong length used to index past its end and report "subscript out of
    # bounds", which says nothing about the argument or about how many
    # parameters there are to give limits for.
    n_lim <- if (is.list(limits)) length(limits) else
      if (is.matrix(limits)) nrow(limits) else NA_integer_
    if (is.na(n_lim)) {
      stop(sprintf(paste0("`limits` must be a list of %s (one c(lo, hi) per ",
                          "parameter) or a matrix with %s, not %s."),
                   n_things(d, "element"), n_things(d, "row"),
                   describe_value(limits)),
           call. = FALSE)
    }
    if (n_lim != d) {
      stop(sprintf(paste0("`limits` has %s but `samples` has %s. Give one ",
                          "c(lo, hi) per parameter, in column order."),
                   n_things(n_lim, if (is.list(limits)) "element" else "row"),
                   n_things(d, "parameter")),
           call. = FALSE)
    }
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
      ggdensity::geom_hdr(fill = col, ...)
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
    progress = FALSE,
    labeller = ggplot2::label_parsed
  ) + ggplot2::theme_minimal()

  print(p)
  invisible(p)
}

#' 99% Monte-Carlo binomial band for a calibration curve
#'
#' At each nominal level, how far empirical coverage over `n` trials can
#' wander from the diagonal by chance alone, under a `Binomial(n, nominal)`
#' model of the count that lands inside its interval. [plot_coverage()] and
#' [plot_tarp()] shade this as the ribbon behind their curve; [plot_sbc()]
#' uses the same interval, scaled back up to a count, for its histogram's
#' dashed reference lines. One formula in one place means all three figures'
#' bands agree on what "real" departure from calibration looks like.
#'
#' @param nominal Nominal level(s) the band is evaluated at.
#' @param n Number of trials the empirical coverage was computed over.
#' @param level Two-sided coverage of the band, e.g. `0.99` for the
#'   `c(0.005, 0.995)` binomial quantiles.
#' @return A data frame with columns `nominal`, `lo` and `hi`, all as
#'   fractions of `n`.
#' @keywords internal
binom_band <- function(nominal, n, level = 0.99) {
  tail <- (1 - level) / 2
  data.frame(
    nominal = nominal,
    lo = stats::qbinom(tail, n, nominal) / n,
    hi = stats::qbinom(1 - tail, n, nominal) / n
  )
}

#' Shared skeleton for a nominal-vs-empirical calibration plot
#'
#' The Monte-Carlo band, the diagonal reference line, the equal-aspect
#' `[0, 1] x [0, 1]` frame and the minimal theme are the same figure in
#' [plot_coverage()] and [plot_tarp()]; only the curve drawn on top (one per
#' parameter, with a legend, for [plot_coverage()]; a single curve for
#' [plot_tarp()]) and the title differ. This builds the shared part and
#' leaves the caller to add its own [ggplot2::geom_line()] and
#' [ggplot2::labs()] on top.
#'
#' @param df Base data for the plot, inherited by whatever curve layer the
#'   caller adds afterward.
#' @param band A data frame with `nominal`, `lo` and `hi` columns (see
#'   [binom_band()]), shaded as the Monte-Carlo uncertainty ribbon.
#' @param xlab,ylab Axis labels.
#' @return A `ggplot` object carrying the shared layers, ready for the
#'   caller's curve, scale and title.
#' @keywords internal
calibration_plot <- function(df, band, xlab, ylab) {
  ggplot2::ggplot(df) +
    ggplot2::geom_ribbon(data = band, ggplot2::aes(x = .data$nominal, ymin = .data$lo, ymax = .data$hi),
                         fill = "grey60", alpha = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = "dashed") +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_minimal()
}

#' Plot an SBC rank histogram
#'
#' Uniform bars indicate calibration; a U shape means the posterior is too
#' narrow (overconfident); an inverted-U means it is too wide.
#'
#' @param sbc_result An `nsbi_sbc` object from [sbc()].
#' @param param Which parameter to plot: an index (default 1), or a parameter
#'   name when [sbc()] was run against a fit with named parameters (see
#'   [prior_uniform()]/[prior_normal()]). The title uses that parameter's name
#'   when there is one, rendered as a plotmath symbol if the name parses as R
#'   syntax.
#' @param bins Number of histogram bins.
#' @return A `ggplot` object (also drawn as a side effect), invisibly.
#' @export
plot_sbc <- function(sbc_result, param = 1L, bins = 20L) {
  stopifnot(inherits(sbc_result, "nsbi_sbc"))
  param <- check_index(param, "param", colnames(sbc_result$ranks),
                       ncol(sbc_result$ranks))
  require_ggplot2()
  r <- sbc_result$ranks[, param]
  L <- sbc_result$n_posterior_samples
  breaks <- seq(0, L, length.out = bins + 1L)
  expected <- sbc_result$n_sbc / bins
  band <- binom_band(1 / bins, sbc_result$n_sbc)
  ci <- c(band$lo, band$hi) * sbc_result$n_sbc

  param_name <- colnames(sbc_result$ranks)[param]
  title <- if (is.null(param_name)) {
    sprintf("SBC ranks: parameter %d", param)
  } else {
    bquote("SBC ranks:" ~ .(math_expr(param_name)))
  }

  p <- ggplot2::ggplot(data.frame(rank = r), ggplot2::aes(x = .data$rank)) +
    ggplot2::geom_histogram(breaks = breaks, fill = "grey80", colour = "white") +
    ggplot2::geom_hline(yintercept = expected, colour = "firebrick", linewidth = 0.7, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = ci, colour = "firebrick", linewidth = 0.5, linetype = "dotted") +
    ggplot2::labs(title = title, x = "rank of true value", y = "count") +
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
#' @param levels Nominal credibility levels to evaluate, each strictly between
#'   0 and 1. Passed to [expected_coverage()], which checks them.
#' @return A `ggplot` object (also drawn as a side effect), invisibly.
#' @export
plot_coverage <- function(sbc_result, levels = seq(0.05, 0.95, by = 0.05)) {
  stopifnot(inherits(sbc_result, "nsbi_sbc"))
  require_ggplot2()
  cov <- expected_coverage(sbc_result, levels = levels)
  d <- ncol(cov) - 1L
  band <- binom_band(cov$nominal, sbc_result$n_sbc)
  long <- do.call(rbind, lapply(seq_len(d), function(j) {
    data.frame(nominal = cov$nominal, empirical = cov[[j + 1L]],
              parameter = colnames(cov)[j + 1L])
  }))

  p <- calibration_plot(long, band, "nominal credibility level", "empirical coverage") +
    ggplot2::geom_line(ggplot2::aes(x = .data$nominal, y = .data$empirical, colour = .data$parameter),
                       linewidth = 0.8) +
    ggplot2::scale_colour_discrete(labels = math_labels) +
    ggplot2::labs(title = "Expected coverage", colour = NULL)

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
  band <- binom_band(lev, tarp_result$n_tarp)
  line <- data.frame(nominal = lev, ecp = tarp_result$ecp)

  p <- calibration_plot(line, band, "nominal credibility level", "expected coverage probability") +
    ggplot2::geom_line(ggplot2::aes(x = .data$nominal, y = .data$ecp),
                       colour = "steelblue", linewidth = 0.8) +
    ggplot2::labs(title = "TARP coverage")

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
#' @param labels Optional labels for the data dimensions. Defaults to
#'   `colnames(pred)` (set automatically when the simulator names its output
#'   columns) or `x[1]`, `x[2]`, .... Labels that parse as R syntax render as
#'   their plotmath symbol.
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
  labels <- math_safe_text(labels)

  long <- data.frame(
    value = as.vector(pred),
    dim = factor(rep(labels, each = nrow(pred)), levels = labels)
  )
  obs_df <- data.frame(dim = factor(labels, levels = labels), obs = x_obs)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = bins, fill = "grey80", colour = "white") +
    ggplot2::geom_vline(data = obs_df, ggplot2::aes(xintercept = .data$obs),
                        colour = "firebrick", linewidth = 0.8) +
    ggplot2::facet_wrap(~dim, scales = "free", labeller = ggplot2::label_parsed) +
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

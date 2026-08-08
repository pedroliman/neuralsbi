#' Save and reload a fitted model
#'
#' A fit whose estimator is `"maf"`, `"mdn"` or `"nsf"` holds a `torch` module,
#' and a torch module is an external pointer. `saveRDS()` writes the pointer,
#' not the network: the file reloads without complaint and the object prints
#' normally, but the first call that touches the network fails with
#' `external pointer is not valid`. `save_npe()` and `load_npe()` are the
#' round trip that works.
#'
#' `save_npe()` writes one `.rds` file holding the network's weights (via
#' [torch::torch_save()] on its `state_dict`) alongside everything else the fit
#' carries as ordinary R objects: the prior, the standardization centers and
#' scales, parameter and outcome names, the simulation count, and the training
#' history. `load_npe()` rebuilds the network from the recorded architecture
#' and restores the weights, returning an `nsbi_npe` that behaves exactly like
#' the one you trained.
#'
#' A `"linear_gaussian"` fit holds no torch objects and round-trips through
#' `saveRDS()` unharmed; `save_npe()` accepts it anyway, so saving code does
#' not have to know which estimator was used.
#'
#' Weights are saved, not code. A fit saved by one version of `neuralsbi` loads
#' into a later one as long as the estimator's architecture has not changed;
#' `load_npe()` reports the version that wrote the file when the rebuild fails.
#'
#' `save_nle()` and `load_nle()` are aliases; both pairs handle either kind of
#' fit, and the two names exist only so calling code reads the way the fit was
#' made.
#'
#' `save_npe()` always writes the network's weights from the CPU, regardless
#' of what device the fit trained on, so the file itself does not pin you to
#' the machine (or the GPU) that produced it. `load_npe()` defaults to
#' rebuilding the network on the CPU for the same reason -- it is the one
#' device every machine has -- and takes a `device` argument for when you do
#' want it back on a GPU.
#'
#' @param fit An `nsbi_npe` object from [npe()] or [npe_sequential()], or an
#'   `nsbi_nle` object from [nle()].
#' @param path File to write to (or read from). The convention is `.rds`.
#'   `load_npe()` says so when there is no such file.
#' @param device Where to rebuild the network: `"cpu"` (the default), `"cuda"`
#'   or `"mps"`. See the `device` argument of [npe()]/[nle()]; the same
#'   fallback-with-a-warning applies when the requested backend is not
#'   available here.
#' @return `save_npe()` returns `path` invisibly. `load_npe()` returns the fit.
#'
#' @examples
#' prior <- prior_uniform(c(mu = -2, nu = -2), c(mu = 2, nu = 2))
#' simulator <- function(mu, nu) c(a = mu + rnorm(1, sd = 0.1),
#'                                 b = nu + rnorm(1, sd = 0.1))
#' fit <- npe(prior, simulator, n_simulations = 500,
#'            density_estimator = "linear_gaussian")
#'
#' path <- tempfile(fileext = ".rds")
#' save_npe(fit, path)
#' fit2 <- load_npe(path)
#' sample(posterior(fit2, x_obs = c(0.8, 0.6)), 100)
#' unlink(path)
#' @name save_npe
NULL

#' @rdname save_npe
#' @export
save_npe <- function(fit, path) {
  stopifnot(inherits(fit, "nsbi_npe") || inherits(fit, "nsbi_nle"))
  check_path(path)
  net <- fit$de$net
  weights <- NULL
  if (!is.null(net)) {
    if (!torch_net_alive(net)) {
      stop("This fit's network is a dangling external pointer, so there are ",
           "no weights to save. It came from readRDS(); refit, or reload the ",
           "original with load_npe().", call. = FALSE)
    }
    require_torch()
    tmp <- tempfile(fileext = ".pt")
    on.exit(unlink(tmp), add = TRUE)
    # Saved to CPU regardless of the device the fit trained on, so the file
    # itself is portable: reloading it on a machine with no GPU (or a
    # different one) never has to allocate a CUDA/MPS tensor it cannot honor.
    state <- lapply(net$state_dict(), function(t) t$cpu())
    torch::torch_save(state, tmp)
    weights <- readBin(tmp, "raw", n = file.size(tmp))
  }
  bundle <- list(
    nsbi_save_format = 1L,
    package_version = as.character(utils::packageVersion("neuralsbi")),
    saved_at = Sys.time(),
    fit = de_drop_net(fit),
    weights = weights
  )
  saveRDS(bundle, path)
  invisible(path)
}

#' @rdname save_npe
#' @export
load_npe <- function(path, device = "cpu") {
  # readRDS() on a path that is not a file says "cannot open the connection",
  # and the file it could not open is named only in the warning that comes with
  # it. save_npe() has checked its path since it was written; this is the same
  # check, plus the file has to be there.
  check_path(path, must_exist = TRUE)
  check_device(device)
  bundle <- readRDS(path)
  if (!is.list(bundle) || !identical(bundle$nsbi_save_format, 1L)) {
    stop(sprintf("'%s' was not written by save_npe().", path), call. = FALSE)
  }
  fit <- bundle$fit
  if (is.null(bundle$weights)) return(fit)

  require_torch()
  device <- resolve_device(device)
  tmp <- tempfile(fileext = ".pt")
  on.exit(unlink(tmp), add = TRUE)
  writeBin(bundle$weights, tmp)
  # Rebuilt (and its saved state loaded) on the CPU, matching how it was
  # saved, then moved to the requested device as a separate step -- state
  # dicts load into a network on the same device they describe.
  net <- de_rebuild_net(fit$de)
  state <- torch::torch_load(tmp)
  tryCatch(
    net$load_state_dict(state),
    error = function(e) {
      stop(sprintf(
        paste0("Could not restore the network saved by neuralsbi %s: %s\n",
               "The estimator's architecture has changed since the fit was ",
               "saved; retrain with the current version."),
        bundle$package_version %||% "(unknown)", conditionMessage(e)),
        call. = FALSE)
    }
  )
  net$to(device = device)
  net$eval()
  fit$de$net <- net
  fit$de$device <- device
  fit
}

#' @rdname save_npe
#' @export
save_nle <- function(fit, path) save_npe(fit, path)

#' @rdname save_npe
#' @export
load_nle <- function(path, device = "cpu") load_npe(path, device)

#' The fit without its torch module, for R-level serialization
#' @keywords internal
de_drop_net <- function(fit) {
  fit$de$net <- NULL
  fit
}

#' Rebuild an estimator's network from the architecture recorded on the fit
#'
#' The one place that knows how to turn a stored estimator back into a torch
#' module. Every field it reads is set by the matching `fit_*()`, so adding an
#' estimator means adding a branch here.
#' @keywords internal
de_rebuild_net <- function(de) {
  kind <- class(de)[1L]
  switch(
    kind,
    nsbi_de_mdn = mdn_module(de$dim_x, de$dim_theta, de$n_components,
                             de$hidden, de$embedding)(),
    nsbi_de_maf = maf_module(de$dim_x, de$dim_theta, de$n_transforms,
                             de$hidden, de$embedding)(),
    nsbi_de_nsf = nsf_module(de$dim_x, de$dim_theta, de$n_transforms,
                             de$hidden, de$n_bins, de$tail_bound,
                             de$embedding)(),
    stop(sprintf("Cannot rebuild a network for estimator class '%s'.", kind),
         call. = FALSE)
  )
}

#' Is this torch module still backed by a live pointer?
#'
#' `readRDS()` on a torch-backed fit returns a module whose external pointer is
#' nil. Nothing about the R object says so, so the only way to find out is to
#' touch a tensor.
#' @keywords internal
torch_net_alive <- function(net) {
  if (is.null(net)) return(TRUE)
  tryCatch({
    params <- net$parameters
    if (length(params) > 0L) invisible(dim(params[[1L]]))
    TRUE
  }, error = function(e) FALSE)
}

#' Fail at the door, not three calls later, on a fit that lost its network
#' @keywords internal
check_fit_alive <- function(fit) {
  if (torch_net_alive(fit$de$net)) return(invisible(TRUE))
  stop("This fit's neural network is no longer usable: its torch module is a ",
       "dangling external pointer.\n",
       "A torch-backed fit does not survive saveRDS()/readRDS(). Save it with ",
       "save_npe(fit, path) and reload it with load_npe(path).",
       call. = FALSE)
}

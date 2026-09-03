#' Neural Ratio Estimation (NRE)
#'
#' `nre()` learns the third factorization of the joint. [npe()] learns the
#' posterior and [nle()] learns the likelihood; `nre()` learns neither density
#' but their ratio, \eqn{r(\theta, x) = p(x \mid \theta) / p(x)}, by training a
#' binary classifier to tell \eqn{(\theta, x)} pairs drawn from the joint apart
#' from pairs whose parameter came from a different simulation. The posterior
#' follows from Bayes' rule,
#' \eqn{p(\theta \mid x) \propto r(\theta, x)\,p(\theta)}, and is sampled with
#' MCMC by [posterior()].
#'
#' @section Why a classifier instead of a density:
#'
#' A density estimator has to spend capacity describing the shape of a
#' distribution -- normalization, tails, the lot -- even where that shape does
#' not affect the posterior. A classifier only has to say which of two
#' parameter values explains the data better, and the optimal classifier's
#' logit *is* the log ratio. That makes NRE the natural choice when the data
#' are high-dimensional or awkward to model directly (discrete counts, mixed
#' types, anything a flow handles badly) but easy to discriminate.
#'
#' Like [nle()] and unlike [npe()], the ratio is learned for a single
#' observation, so the log-likelihood of \eqn{n} independent trials from the
#' same parameter is a sum:
#'
#' \deqn{\log \frac{p(x_1, \ldots, x_n \mid \theta)}{\prod_i p(x_i)} = \sum_{i=1}^{n} \log r(\theta, x_i).}
#'
#' Train once, condition on as many trials as you like. The price is the same
#' one [nle()] pays: posterior draws cost an MCMC run rather than a forward
#' pass.
#'
#' @section The training objective:
#'
#' The default is the *atomic* loss of Durkan et al. (2020), which is what
#' `sbi`'s `NRE` (an alias for `NRE_B`) trains. For each simulation
#' \eqn{(\theta_i, x_i)} in a minibatch, `num_atoms - 1` contrasting parameters
#' are taken from the other simulations in that batch, and the classifier is
#' scored on a `num_atoms`-way softmax over which of them produced \eqn{x_i}:
#'
#' \deqn{\mathcal{L} = -\frac{1}{b}\sum_{i} \left[ f(\theta_i, x_i) - \log \sum_{k} \exp f(\theta_{ik}, x_i) \right].}
#'
#' More atoms mean a harder discrimination problem and a sharper ratio, at a
#' linear cost in forward passes per epoch. `num_atoms = 10` is `sbi`'s
#' default and this one. The softmax is invariant to adding any function of
#' \eqn{x} to \eqn{f}, so the learned ratio is calibrated up to a constant at
#' each fixed observation -- exactly what a posterior needs, and the reason
#' [log_prob()] on an NRE posterior is unnormalized.
#'
#' @section Standardization has no Jacobian here:
#'
#' The estimators in [npe()] and [nle()] train on z-scored data and need a
#' change-of-variables term to report densities in the original units. A ratio
#' needs none: both \eqn{p(x \mid \theta)} and \eqn{p(x)} pick up the same
#' Jacobian factor and it cancels, so `log_ratio()` in standardized space is
#' already the ratio in the units the simulator returned.
#'
#' @inheritParams npe
#' @param classifier One of `"resnet"` (a residual MLP, needs `torch`; the
#'   default, matching Python `sbi`), `"mlp"` (a plain MLP, needs `torch`),
#'   `"linear"` (a single linear layer on the raw inputs, needs `torch`), or
#'   `"logistic"` (a closed-form logistic regression on quadratic features, no
#'   `torch`), or a function `function(theta, x)` returning a fitted ratio
#'   estimator -- one whose class has a `de_log_ratio()` method, which is the
#'   only thing the rest of the pipeline asks of it.
#' @param num_atoms Number of parameter values the classifier compares per
#'   simulation: one true and `num_atoms - 1` contrasts. Clamped to the
#'   minibatch size, as in `sbi`.
#' @param hidden Width of the classifier's hidden layers. One number, not a
#'   per-layer vector as in [npe()]: `sbi`'s classifiers use a single width
#'   throughout, and `n_blocks` is what sets the depth.
#' @param n_blocks Depth of the classifier: residual blocks for `"resnet"`,
#'   hidden layers for `"mlp"`. Ignored by `"linear"` and `"logistic"`.
#'   `sbi` fixes its MLP at two hidden layers and only lets `num_blocks` reach
#'   the residual net; here the one argument sets both.
#' @param embedding_net Optional summary network built with [embedding_mlp()].
#'   The classifier then sees \eqn{(\theta, f_\psi(x))}, with the embedding
#'   trained jointly. Ignored (with a warning) by `"logistic"`.
#'
#' @return An object of class `nsbi_nre`. Evaluate the learned ratio with
#'   [log_ratio()] or turn it into a posterior with [posterior()].
#'
#' @seealso [log_ratio()] to evaluate the ratio, [posterior()] to sample it,
#'   [nle()] for the likelihood factorization and [npe()] for the posterior
#'   one.
#'
#' @references Hermans, J., Begy, V. and Louppe, G. (2020). Likelihood-free
#'   MCMC with Amortized Approximate Ratio Estimators. *ICML*.
#'   \doi{10.48550/arXiv.1903.04057}
#'
#'   Durkan, C., Murray, I. and Papamakarios, G. (2020). On Contrastive
#'   Learning for Likelihood-free Inference. *ICML*.
#'   \doi{10.48550/arXiv.2002.03712}
#'
#' @examples
#' # One noisy measurement per simulator call; the observation is 50 of them.
#' prior <- prior_uniform(c(mu = -3), c(mu = 3))
#' simulator <- function(mu) c(y = rnorm(1, mu, 0.5))
#'
#' fit <- nre(prior, simulator, n_simulations = 2000, classifier = "logistic")
#'
#' x_obs <- matrix(rnorm(50, mean = 1, sd = 0.5), ncol = 1)
#' log_ratio(fit, theta = c(1), x = x_obs)
#' @export
nre <- function(prior, simulator = NULL, n_simulations = 1000,
                sim_args = list(), theta = NULL, x = NULL,
                classifier = c("resnet", "mlp", "linear", "logistic"),
                num_atoms = 10L, hidden = 50L, n_blocks = 2L,
                embedding_net = NULL,
                max_epochs = 2000L, batch_size = 200L, lr = 5e-4,
                validation_fraction = 0.1, patience = 20L,
                n_restarts = 1L, clip_grad_norm = 5,
                standardize = TRUE, device = "cpu",
                seed = NULL, verbose = FALSE) {
  # See npe(): everything here is checked before the simulator runs, so a typo
  # does not cost the simulation budget first.
  if (is.function(classifier)) {
    check_function(classifier, "classifier", what = "theta and x matrices")
  } else {
    classifier <- match.arg(classifier)
  }
  device <- check_device_arg(device)
  check_prior(prior)
  if (!is.null(embedding_net) && !inherits(embedding_net, "nsbi_embedding")) {
    stop("`embedding_net` must be built with embedding_mlp().", call. = FALSE)
  }
  if (!is.null(embedding_net) && identical(classifier, "logistic")) {
    warning("`embedding_net` is ignored by the logistic classifier.",
            call. = FALSE)
  }
  hidden <- check_count(hidden, "hidden")
  n_blocks <- check_count(n_blocks, "n_blocks")
  num_atoms <- check_count(
    num_atoms, "num_atoms", min = 2L,
    why = "since one atom is the true parameter and the rest are contrasts")
  # min_val_rows only binds for the neural classifiers fit_nre_net() builds:
  # the closed-form "logistic" fit never splits off a validation set, and a
  # caller-supplied classifier function may not either, so both keep the
  # default floor of 1L (see fit_nre_net(), GitHub #188). `n` is filled in
  # whenever it is already known exactly -- from `theta` when it was passed
  # directly, or from `n_simulations` when it is a valid count -- so the
  # common trigger, nre(prior, simulator, n_simulations = 15), is caught here
  # rather than after the simulation budget is spent.
  min_val_rows <- if (is.function(classifier) || identical(classifier, "logistic")) {
    1L
  } else {
    2L
  }
  n_hint <- if (!is.null(theta) && !is.null(x)) {
    tryCatch(nrow(as_theta_matrix(theta, prior$dim)), error = function(e) NULL)
  } else if (is.numeric(n_simulations) && length(n_simulations) == 1L &&
             is.finite(n_simulations) && n_simulations >= 1) {
    n_simulations
  } else {
    NULL
  }
  check_train_controls(max_epochs, batch_size, lr, validation_fraction,
                       patience, n_restarts, clip_grad_norm, n = n_hint,
                       min_val_rows = min_val_rows)
  # See npe(): "resnet"/"mlp"/"linear" need torch (the default is "resnet");
  # only "logistic" and a caller-supplied classifier do not. Fails before the
  # simulator runs rather than after (GitHub #250).
  check_torch_for_estimator(
    classifier, c("resnet", "mlp", "linear"), device,
    what = "This classifier",
    alternative = paste("Alternatively use classifier = \"logistic\" for a",
                        "torch-free baseline."))

  # See npe(): seed the base RNG here, unconditionally of which
  # prepare_simulations() branch runs, so train_restarts()'s train/validation
  # split and minibatch order are reproducible even with pre-computed
  # theta/x (GitHub #213).
  if (!is.null(seed)) set.seed(seed)
  prep <- prepare_simulations(prior, simulator, n_simulations, sim_args,
                              theta, x, standardize, seed, verbose)

  re <- fit_ratio_estimator(
    classifier, prep$theta_z, prep$x_z,
    hidden = hidden, n_blocks = n_blocks, num_atoms = num_atoms,
    embedding = embedding_net, max_epochs = max_epochs,
    batch_size = batch_size, lr = lr, validation_fraction = validation_fraction,
    patience = patience, n_restarts = n_restarts,
    clip_grad_norm = clip_grad_norm, seed = seed, verbose = verbose,
    device = device
  )

  # The fitted classifier goes in the `de` slot: everything that carries a
  # fitted estimator around -- save_npe(), check_fit_alive(), cat_fit_common(),
  # summary() -- reads fit$de, and a ratio estimator travels exactly the same
  # way.
  new_nsbi_fit(re, prior, prep, "nsbi_nre",
               list(classifier = estimator_label(classifier),
                    num_atoms = num_atoms))
}

#' @export
print.nsbi_nre <- function(x, ...) {
  cat("<nsbi_nre> Neural Ratio Estimation fit\n")
  cat(sprintf("  classifier        : %s  (learns p(x | theta) / p(x))\n",
              x$classifier))
  cat(sprintf("  contrastive atoms : %d\n", x$num_atoms))
  cat_fit_common(x, "save_nre", data_suffix = "  per observation")
  cat("  -> log_ratio(fit, theta, x), posterior(fit, x_obs = ...)\n")
  invisible(x)
}

#' Dispatch to the requested ratio estimator
#'
#' The counterpart of `fit_density_estimator()`. Arguments are forwarded
#' explicitly rather than through `...` and `formals()`, because there are only
#' two targets and they take disjoint arguments: the neural classifiers take
#' the whole training block, the closed-form one takes none of it.
#' @keywords internal
fit_ratio_estimator <- function(classifier, theta_z, x_z, hidden, n_blocks,
                                num_atoms, embedding, max_epochs, batch_size,
                                lr, validation_fraction, patience, n_restarts,
                                clip_grad_norm, seed, verbose, device) {
  if (is.function(classifier)) return(classifier(theta_z, x_z))
  if (classifier == "logistic") {
    return(fit_logistic_ratio(theta_z, x_z, num_atoms = num_atoms,
                              verbose = verbose))
  }
  fit_nre_net(theta_z, x_z, classifier = classifier, hidden = hidden,
              n_blocks = n_blocks, num_atoms = num_atoms,
              embedding = embedding, max_epochs = max_epochs,
              batch_size = batch_size, lr = lr,
              validation_fraction = validation_fraction, patience = patience,
              n_restarts = n_restarts, clip_grad_norm = clip_grad_norm,
              seed = seed, verbose = verbose, device = device)
}

# ---- the ratio contract ----------------------------------------------------

#' Log ratio of a fitted ratio estimator
#'
#' The ratio estimators' counterpart to [de_log_prob()][density_estimator]:
#' one number per row, \eqn{\log r(\theta, x)}, in standardized space (which
#' for a ratio is also the original space -- see [nre()]). A one-row `x` is
#' broadcast against a taller `theta`, as `de_log_prob()` does.
#'
#' @param de A fitted ratio estimator.
#' @param theta Standardized parameters.
#' @param x Standardized data.
#' @return A numeric vector with one entry per row.
#' @keywords internal
de_log_ratio <- function(de, theta, x) UseMethod("de_log_ratio")

#' The scorer [cross_iid()] needs, with the ratio's argument order
#'
#' [cross_iid()] calls `score(de, target, condition)` -- data first, parameters
#' second, the order [nle()]'s estimator is trained in. A ratio has no target
#' and no condition, and reads better with the parameter first, so the flip
#' happens here rather than in [de_log_ratio()]'s signature.
#' @keywords internal
nre_score <- function(de, x, theta) de_log_ratio(de, theta, x)

#' Summed log ratio over independent observations, with the observation fixed
#'
#' [de_iid_evaluator()]'s counterpart for a ratio estimator. There is no fast
#' path to specialize: the classifier sees `(theta, x_i)` jointly, so every
#' pair costs a forward pass however the loop is arranged.
#' @inheritParams de_iid_evaluator
#' @return `function(theta)` giving one summed log ratio per row of `theta`.
#' @keywords internal
nre_iid_evaluator <- function(de, x, max_batch = 1e5) {
  iid_evaluator(de, x, max_batch, nre_score)
}

#' The `n_theta x n_obs` matrix of log ratios
#'
#' [de_log_lik_iid()]'s counterpart for a ratio estimator, and like
#' [nre_iid_evaluator()] a plain function rather than a generic: there is no
#' fast path for any classifier to specialize.
#' @inheritParams de_log_lik_iid
#' @return An `n_theta x n_obs` matrix of log ratios.
#' @keywords internal
nre_log_ratio_iid <- function(de, x, theta, max_batch = 1e5) {
  iid_matrix(de, x, theta, max_batch, nre_score)
}

#' @export
surrogate_ops.nsbi_nre <- function(fit) {
  # log_jac = FALSE: the standardization Jacobian cancels between the two
  # densities in the ratio, so there is nothing to correct for. See ?nre.
  list(matrix_fn = nre_log_ratio_iid, evaluator = nre_iid_evaluator,
       log_jac = FALSE)
}

#' Evaluate a learned likelihood-to-evidence ratio
#'
#' `log_ratio()` evaluates the ratio learned by [nre()],
#' \eqn{\log r(\theta, x) = \log p(x \mid \theta) - \log p(x)}.
#'
#' Rows of `x` are treated as **independent observations from the same
#' parameter**, so by default the result sums over them. That sum is the
#' log-likelihood of the whole data set up to an additive constant that does
#' not depend on `theta`, which is all a posterior needs and all NRE learns:
#' the evidence term \eqn{\log p(x)} is never estimated separately, so
#' `log_ratio()` is not a log-likelihood you can compare across observations.
#' Differences between two `theta` at the same `x` are exactly the differences
#' in log-likelihood.
#'
#' @param fit An `nsbi_nre` fit from [nre()].
#' @param theta Parameter values: a numeric vector (one parameter set) or an
#'   `n_theta x dim_theta` matrix.
#' @param x Observed data: a numeric vector (one observation) or an
#'   `n_obs x dim_x` matrix whose rows are independent observations.
#' @param sum_iid Sum the log ratio over the rows of `x` (the default). Set
#'   `FALSE` to get the per-observation values instead.
#' @param max_batch Largest number of `(theta, x)` pairs evaluated in one call
#'   to the classifier. Only affects memory and speed.
#' @param ... Unused, for S3 consistency.
#'
#' @return With `sum_iid = TRUE`, a numeric vector with one entry per row of
#'   `theta`. With `sum_iid = FALSE`, an `n_theta x n_obs` matrix.
#'
#' @seealso [posterior()] to turn the ratio into posterior draws, [log_lik()]
#'   for the [nle()] equivalent that *is* a calibrated log-likelihood.
#'
#' @examples
#' prior <- prior_uniform(c(mu = -3), c(mu = 3))
#' fit <- nre(prior, function(mu) c(y = rnorm(1, mu, 0.5)),
#'            n_simulations = 1000, classifier = "logistic")
#'
#' x_obs <- matrix(rnorm(20, mean = 1, sd = 0.5), ncol = 1)
#' grid <- matrix(seq(-2, 2, length.out = 5), ncol = 1)
#' log_ratio(fit, grid, x_obs)
#' @export
log_ratio <- function(fit, theta, x, ...) UseMethod("log_ratio")

#' @rdname log_ratio
#' @export
log_ratio.nsbi_nre <- function(fit, theta, x, sum_iid = TRUE,
                               max_batch = 1e5, ...) {
  surrogate_score(fit, theta, x, sum_iid, max_batch)
}

# ---- neural classifiers ----------------------------------------------------

#' One pre-activation residual block, as in `nflows`' `ResidualNet`
#'
#' `relu -> linear -> relu -> linear`, added back to the input. The second
#' linear layer starts near zero so the block starts as the identity, the same
#' trick [made_module()] uses on its output heads and the reason a deep stack
#' trains at all from a cold start.
#' @keywords internal
nre_residual_block <- function(features) {
  torch::nn_module(
    classname = "nsbi_nre_block",
    initialize = function() {
      self$linear1 <- torch::nn_linear(features, features)
      self$linear2 <- torch::nn_linear(features, features)
      torch::nn_init_uniform_(self$linear2$weight, -1e-3, 1e-3)
      torch::nn_init_uniform_(self$linear2$bias, -1e-3, 1e-3)
    },
    forward = function(input) {
      h <- self$linear1(torch::nnf_relu(input))
      input + self$linear2(torch::nnf_relu(h))
    }
  )
}

#' Build the classifier torch module: `(theta, x)` -> one logit
#'
#' Three architectures behind one module, matching `sbi`'s `"resnet"`,
#' `"mlp"` and `"linear"` classifiers. `sbi` normalizes between the hidden
#' layers of its MLP (`nn.LayerNorm` by default); this one does not, keeping
#' the trunk the same plain `linear/relu` stack every other estimator in the
#' package uses ([mlp_layers()]).
#' @keywords internal
nre_module <- function(dim_x, dim_theta, classifier, hidden, n_blocks,
                       embedding = NULL) {
  torch::nn_module(
    classname = "nsbi_nre_net",
    initialize = function() {
      self$has_embedding <- !is.null(embedding)
      if (self$has_embedding) {
        self$embedding <- build_embedding_module(embedding, dim_x)
      }
      in_features <- dim_theta + embedding_output_dim(embedding, dim_x)
      if (classifier == "resnet") {
        self$initial <- torch::nn_linear(in_features, hidden)
        self$blocks <- torch::nn_module_list(
          lapply(seq_len(n_blocks), function(i) nre_residual_block(hidden)()))
        self$final <- torch::nn_linear(hidden, 1L)
      } else if (classifier == "mlp") {
        layers <- mlp_layers(c(in_features, rep(hidden, n_blocks)))
        layers[[length(layers) + 1L]] <- torch::nn_linear(hidden, 1L)
        self$trunk <- do.call(torch::nn_sequential, layers)
      } else {
        self$trunk <- torch::nn_linear(in_features, 1L)
      }
    },
    forward = function(theta, x) {
      h <- torch::torch_cat(list(theta, embed_x(self, x)), dim = 2)
      if (classifier != "resnet") return(self$trunk(h))
      h <- self$initial(h)
      for (k in seq_len(n_blocks)) h <- self$blocks[[k]](h)
      self$final(h)
    }
  )
}

#' Per-row logit of the classifier, as a torch tensor
#'
#' This is the log ratio itself: the module's single output for each
#' `(theta, x)` row.
#' @keywords internal
nre_logit_tensor <- function(net, theta, x) net(theta, x)$squeeze(2)

#' Which rows of a minibatch each simulation is scored against
#'
#' Returns a length-`b * k` vector of row indices, `k` per simulation and
#' running simulation-major: the true parameter for row `i` first, then `k - 1`
#' contrasts drawn without replacement from the *other* rows of the batch. That
#' is `sbi`'s atom construction, done with R's RNG rather than
#' `torch.multinomial` so it is seeded by the same `set.seed()` the training
#' loop's split and batch order already are.
#'
#' `deterministic` freezes the draw. The atoms are resampled every time the
#' loss is evaluated, which is the objective's business during training but
#' makes a poor early-stopping signal: the validation loss would move between
#' epochs because the contrasts changed, not because the classifier did.
#' Freezing it for the validation pass costs nothing and makes the two numbers
#' comparable. (`sbi` resamples there too, and pays for it in noisier
#' stopping.)
#' @keywords internal
nre_atom_rows <- function(b, k, deterministic = FALSE) {
  if (deterministic) return(with_fixed_seed(1L, nre_atom_rows(b, k)))
  out <- integer(b * k)
  for (i in seq_len(b)) {
    # Draw from 1..b-1 and step over i, which is cheaper than rejecting.
    contrast <- sample.int(b - 1L, k - 1L)
    out[((i - 1L) * k + 1L):(i * k)] <- c(i, contrast + (contrast >= i))
  }
  out
}

#' Run `expr` with the RNG parked at `seed`, restoring the stream afterwards
#'
#' The caller's random stream is untouched, so a training run is still
#' reproducible from the seed it was given.
#' @keywords internal
with_fixed_seed <- function(seed, expr) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
            add = TRUE)
  }
  set.seed(seed)
  expr
}

#' The atomic NRE-B objective, in the shape [train_conditional_de()] wants
#'
#' The training engine minimizes `-mean(log_prob_fn(net, theta, x))`, so the
#' per-row quantity returned here is the atomic log-softmax: the true
#' parameter's logit minus the log-sum-exp over its atoms. Maximizing it is
#' minimizing the cross-entropy of picking the true parameter out of the
#' `num_atoms` candidates.
#'
#' `num_atoms` is clamped to the batch size, as `sbi` does, so a short final
#' minibatch or a small validation split does not ask for more contrasts than
#' the batch can supply. A batch of one has no contrast at all and contributes
#' nothing, but it still has to contribute a tensor torch built: a constant of
#' its own making has no graph for `backward()` to walk. [minibatches()] keeps
#' single-row batches out of training; a validation split of one row still
#' arrives here.
#' @keywords internal
nre_atomic_log_prob <- function(num_atoms) {
  force(num_atoms)
  function(net, theta, x) {
    b <- theta$shape[1]
    k <- min(num_atoms, b)
    if (k < 2L) return(nre_logit_tensor(net, theta, x) * 0)
    rows <- nre_atom_rows(b, k, deterministic = !net$training)
    logits <- net(theta[rows, , drop = FALSE],
                  x[rep(seq_len(b), each = k), , drop = FALSE])$view(c(b, k))
    logits[, 1] - torch::torch_logsumexp(logits, dim = 2)
  }
}

#' Train a neural ratio estimator on standardized (theta, x)
#'
#' Passes `min_val_rows = 2L` down to [fit_torch_de()], unlike the MDN/MAF/NSF
#' callers. The atomic objective ([nre_atomic_log_prob()]) needs a second row
#' to contrast the true parameter against; with a split of one row on either
#' side it silently returns a constant zero loss every epoch instead of a real
#' signal -- breaking early stopping when it is the validation side (GitHub
#' #188), and training on zero gradient with no error when it is the training
#' side (GitHub #239). `check_train_controls()` enforces `min_val_rows` on
#' both sides of the split for exactly this reason.
#' @keywords internal
fit_nre_net <- function(theta, x, classifier = "resnet", hidden = 50L,
                        n_blocks = 2L, num_atoms = 10L,
                        max_epochs = 2000L, batch_size = 200L, lr = 5e-4,
                        validation_fraction = 0.1, patience = 20L,
                        n_restarts = 1L, clip_grad_norm = 5, embedding = NULL,
                        seed = NULL, verbose = FALSE, device = "cpu") {
  fit_torch_de(
    theta, x,
    build_net_fn = function(dim_x, dim_theta)
      nre_module(dim_x, dim_theta, classifier, hidden, n_blocks, embedding)(),
    log_prob_fn = nre_atomic_log_prob(num_atoms),
    class = "nsbi_re_net",
    arch = list(classifier = classifier, hidden = hidden, n_blocks = n_blocks,
                num_atoms = num_atoms),
    max_epochs = max_epochs, batch_size = batch_size, lr = lr,
    validation_fraction = validation_fraction, patience = patience,
    n_restarts = n_restarts, clip_grad_norm = clip_grad_norm,
    embedding = embedding, seed = seed, verbose = verbose, device = device,
    min_val_rows = 2L
  )
}

#' @export
de_log_ratio.nsbi_re_net <- function(de, theta, x) {
  de_log_prob_torch(de, theta, x, nre_logit_tensor)
}

# ---- closed-form logistic classifier (pure R) ------------------------------

#' Quadratic feature basis for the logistic ratio estimator
#'
#' `[1, z, vech(z z')]` for `z = (theta, x)`. The basis is chosen so the
#' estimator is *exact* for a linear-Gaussian simulator: there
#' \eqn{\log p(x \mid \theta)} is a quadratic form in \eqn{(\theta, x)}, so the
#' log ratio's parameter dependence lies inside this span and the fit is
#' limited only by estimation error (see [fit_logistic_ratio()] for why the
#' evidence term does not spoil that). It is the regression oracle for [nre()]
#' that `"linear_gaussian"` is for [npe()] and [nle()].
#'
#' It costs `1 + d + d(d+1)/2` columns for `d = dim_theta + dim_x`, so it is a
#' baseline for small models, not a substitute for a neural classifier on wide
#' data.
#' @keywords internal
nre_features <- function(theta, x) {
  z <- cbind(theta, x)
  d <- ncol(z)
  idx <- which(upper.tri(matrix(0, d, d), diag = TRUE), arr.ind = TRUE)
  cbind(1, z, z[, idx[, 1], drop = FALSE] * z[, idx[, 2], drop = FALSE])
}

#' Ridge-penalized logistic regression by iteratively reweighted least squares
#'
#' `stats::glm.fit()` would do this, but it warns and wanders off when the two
#' classes are separable, which a sharp ratio makes easy to hit. The ridge term
#' keeps the normal equations solvable and the coefficients finite, the same
#' role it plays in `fit_linear_gaussian()`, and it is measured against each
#' column's own scale for the same reason: under `standardize = FALSE` the
#' quadratic features carry the fourth power of the data's units, so an
#' absolute 1e-6 is either nothing at all or the only thing left. On a
#' simulator whose output has sd 5e-4 the absolute version shrank the fit to
#' noise; the relative one leaves it alone.
#' @keywords internal
irls_logistic <- function(X, y, ridge = 1e-6, max_iter = 100L, tol = 1e-8) {
  beta <- rep(0, ncol(X))
  for (it in seq_len(max_iter)) {
    eta <- as.numeric(X %*% beta)
    mu <- 1 / (1 + exp(-eta))
    w <- pmax(mu * (1 - mu), 1e-8)
    z <- eta + (y - mu) / w
    XtWX <- crossprod(X * w, X)
    diag(XtWX) <- diag(XtWX) + ridge * ridge_scale(diag(XtWX))
    step <- as.numeric(solve(XtWX, crossprod(X * w, z)))
    delta <- max(abs(step - beta))
    beta <- step
    if (delta < tol) break
  }
  beta
}

#' Fit the closed-form logistic ratio estimator on standardized (theta, x)
#'
#' The atomic objective at two atoms has a closed form. For one simulation
#' \eqn{(\theta_i, x_i)} and one contrast \eqn{\theta_j}, the two-way softmax
#' reduces to
#' \eqn{\log \sigma\!\left(f(\theta_i, x_i) - f(\theta_j, x_i)\right)}, and with
#' \eqn{f = w'\varphi} linear in the features that is a logistic regression on
#' the *difference* \eqn{\varphi(\theta_i, x_i) - \varphi(\theta_j, x_i)} with
#' every label positive. So the whole estimator is one ridge-penalized IRLS on
#' an `n * (num_atoms - 1)` by `ncol(Phi)` design, no optimizer and no torch.
#'
#' Working in differences is what makes the fit *exact* for a linear-Gaussian
#' simulator rather than merely close. The differences cancel every term that
#' depends on `x` alone, including the evidence \eqn{\log p(x)}, which is the
#' one part of the log ratio a quadratic basis cannot represent. What is left
#' is the parameter dependence \eqn{\log p(x \mid \theta)}, which for a
#' linear-Gaussian model lies exactly in the span of [nre_features()]. The
#' price is the one the atomic objective always pays: the level of the fitted
#' ratio at a given `x` is arbitrary.
#'
#' Contrasts come from cyclic shifts of the parameter rows rather than a random
#' draw. The rows are independent prior draws in random order already, so a
#' shift is as good a scramble, and it makes the fit deterministic: the same
#' simulations give the same estimator whether or not a seed was set.
#' @keywords internal
fit_logistic_ratio <- function(theta, x, num_atoms = 10L, ridge = 1e-6,
                               verbose = FALSE) {
  theta <- as_theta_matrix(theta)
  x <- as_theta_matrix(x)
  n <- nrow(theta)
  if (n < 2L) {
    stop("The logistic ratio estimator needs at least 2 simulations to build ",
         "a contrasting pair.", call. = FALSE)
  }
  k <- min(as.integer(num_atoms), n)
  true <- nre_features(theta, x)
  # One block of difference rows per contrast shift. Averaging over several
  # shifts is still the two-atom objective -- it is the same expectation, just
  # estimated from more pairs -- so it lowers the variance of the fit without
  # changing what it converges to.
  diffs <- lapply(seq_len(k - 1L), function(s) {
    rotated <- c(seq.int(s + 1L, n), seq_len(s))
    true - nre_features(theta[rotated, , drop = FALSE], x)
  })
  X <- do.call(rbind, diffs)
  beta <- irls_logistic(X, rep(1, nrow(X)), ridge = ridge)
  verbose_cat(verbose, sprintf(
    "[logistic] fitted on %d sims, %d params, %d data dims, %d features\n",
    n, ncol(theta), ncol(x), ncol(true)))
  structure(
    list(beta = beta, dim_theta = ncol(theta), dim_x = ncol(x),
         num_atoms = k),
    class = c("nsbi_re_logistic", "nsbi_de")
  )
}

#' @export
de_log_ratio.nsbi_re_logistic <- function(de, theta, x) {
  theta <- as_theta_matrix(theta, de$dim_theta)
  x <- as_theta_matrix(x, de$dim_x)
  if (nrow(x) == 1L && nrow(theta) > 1L) {
    x <- matrix(x, nrow = nrow(theta), ncol = ncol(x), byrow = TRUE)
  }
  as.numeric(nre_features(theta, x) %*% de$beta)
}

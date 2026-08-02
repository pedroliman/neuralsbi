test_that("check_matrix() names the argument it is complaining about", {
  expect_error(check_matrix(matrix(0, 3, 5), 2L, "theta",
                            "one parameter per column"),
               "`theta` must have 2 columns \\(one parameter per column\\)")
  expect_error(check_matrix(matrix(0, 3, 5), 2L, "x"), "^`x` ")
  expect_error(check_matrix(NULL, 2L, "x_obs"), "`x_obs` is missing")
  expect_error(check_matrix(letters[1:2], 2L, "theta"), "`theta` must be numeric")
  expect_error(check_matrix(list(1, 2), 2L, "theta"), "not a list")
  expect_error(check_matrix(array(0, c(2, 2, 2)), 2L, "x"),
               "3-dimensional array")
  expect_error(check_matrix(numeric(0), 2L, "theta"), "`theta` is empty")
  expect_error(check_matrix(data.frame(a = 1, b = "x"), 2L, "theta"),
               "non-numeric columns: b")
})

test_that("check_matrix() errors on a wrong-length vector instead of recycling", {
  # as_theta_matrix() reshapes this into a 2 x 2 matrix by recycling; two
  # answers from three numbers is the failure this validator exists for.
  expect_error(check_matrix(c(0.1, 0.2, 0.3), 2L, "theta"),
               "length-3 vector")
  expect_error(check_matrix(c(1, 2, 3, 4), 2L, "theta"), "length-4 vector")
  expect_silent(as_theta_matrix(c(1, 2, 3, 4), 2L))
})

test_that("check_matrix() reads a bare vector as one row and keeps names", {
  out <- check_matrix(c(a = 1, b = 2), 2L, "theta")
  expect_equal(dim(out), c(1L, 2L))
  expect_equal(colnames(out), c("a", "b"))
  expect_type(out, "double")

  # With no required width, any bare vector is a single row.
  expect_equal(dim(check_matrix(c(1, 2, 3), arg = "theta")), c(1L, 3L))
})

test_that("check_matrix() passes matrices and data frames through", {
  m <- matrix(1:6, ncol = 2, dimnames = list(NULL, c("a", "b")))
  expect_equal(check_matrix(m, 2L, "theta"), matrix(as.double(1:6), ncol = 2,
    dimnames = list(NULL, c("a", "b"))))
  expect_equal(check_matrix(data.frame(a = 1:3, b = 4:6), 2L, "theta"),
               matrix(as.double(1:6), ncol = 2,
                      dimnames = list(NULL, c("a", "b"))))
  expect_equal(dim(check_matrix(matrix(0, 0, 2), 2L, "theta")), c(0L, 2L))
})

test_that("check_matrix() suggests a transpose when the rows match", {
  expect_error(check_matrix(matrix(0, 2, 100), 2L, "theta"),
               "Did you mean to transpose it")
  msg <- tryCatch(check_matrix(matrix(0, 3, 100), 2L, "theta"),
                  error = conditionMessage)
  expect_false(grepl("transpose", msg))
})

test_that("check_count() takes one whole number at or above the bound", {
  expect_identical(check_count(5, "n_simulations"), 5L)
  expect_identical(check_count(0, "warmup", min = 0L), 0L)

  expect_error(check_count(2.5, "n_rounds"),
               "`n_rounds` must be a single whole number of at least 1, not 2.5")
  expect_error(check_count(0, "n_rounds"), "`n_rounds` must be")
  expect_error(check_count(NA, "n_rounds"), "not NA")
  expect_error(check_count(c(1, 2), "n_rounds"), "a length-2 numeric vector")
  expect_error(check_count("3", "n_rounds"), "a character value")
  expect_error(check_count(NULL, "n_rounds"), "not NULL")
  expect_error(check_count(Inf, "n_rounds"), "not Inf")
  expect_error(check_count(1, "thin", min = 2L, why = "for the reason given"),
               "at least 2 for the reason given")
})

test_that("check_prob() takes one number in the unit interval", {
  expect_identical(check_prob(0.5, "alpha"), 0.5)
  expect_identical(check_prob(0, "alpha", open = FALSE), 0)

  expect_error(check_prob(0, "alpha"),
               "`alpha` must be a single number strictly between 0 and 1, not 0")
  expect_error(check_prob(1, "alpha"), "strictly between 0 and 1")
  expect_error(check_prob(1.5, "alpha", open = FALSE),
               "between 0 and 1 inclusive")
  expect_error(check_prob(NA_real_, "alpha"), "not NA")
  expect_error(check_prob(c(0.1, 0.2), "alpha"), "a length-2 numeric vector")
})

test_that("check_positive() takes one finite number above zero", {
  expect_identical(check_positive(2L, "lr"), 2)

  expect_error(check_positive(0, "lr"),
               "`lr` must be a single positive number, not 0")
  expect_error(check_positive(-1, "lr"), "not -1")
  expect_error(check_positive(Inf, "lr"), "not Inf")
  expect_error(check_positive(NaN, "lr"), "not NaN")
  expect_error(check_positive(NULL, "lr"), "not NULL")
})

test_that("check_prior() checks the class and the dimension", {
  prior <- prior_uniform(c(a = 0, b = 0), c(a = 1, b = 1))
  expect_identical(check_prior(prior, dim = 2L), prior)

  expect_error(check_prior(list(), "prior"),
               "`prior` must be an nsbi_prior object, not list")
  expect_error(check_prior(NULL), "not NULL")
  expect_error(check_prior(list(), "proposal"), "^`proposal` ")
  expect_error(check_prior(prior, dim = 3L),
               "`prior` covers 2 parameters but 3 are expected here")
  expect_error(check_prior(prior, dim = 1L), "but 1 is expected here")
})

test_that("check_finite() names the argument and the first bad entry", {
  expect_silent(check_finite(matrix(1:6, ncol = 2), "theta"))

  m <- matrix(as.double(1:6), ncol = 2)
  m[2, 2] <- NA
  expect_error(check_finite(m, "theta"),
               "`theta` contains 1 non-finite value \\(NA\\), first at row 2, column 2")

  v <- c(1, NaN, Inf, NA)
  expect_error(check_finite(v, "x"),
               "`x` contains 3 non-finite values \\(NA/NaN/Inf\\), first at position 2")
  expect_error(check_finite("a", "x"), "`x` must be numeric")
})

test_that("log_lik() rejects a wrong-length vector by name", {
  fit <- nle(prior_uniform(c(mu = -3), c(mu = 3)),
             function(mu) c(y = stats::rnorm(1, mu, 0.5)), n_simulations = 400,
             density_estimator = "linear_gaussian", seed = 1)

  # dim_theta is 1, so a length-5 vector is not one parameter set. It used to
  # be reshaped into five of them without a word.
  expect_error(log_lik(fit, seq(-2, 2, length.out = 5), matrix(0, 1, 1)),
               "`theta` must have 1 column")
  expect_error(log_lik(fit, 0, matrix(0, 2, 2)), "`x` must have 1 column")
})

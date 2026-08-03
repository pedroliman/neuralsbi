library(testthat)
library(neuralsbi)

# CI sets NEURALSBI_JUNIT to a path so Codecov test analytics has a JUnit XML
# file to read. The default run is unchanged: no env var, no reporter argument,
# so R CMD check still picks its own reporter. CheckReporter rides alongside the
# JUnit one because it is what turns a failure into a non-zero exit.
junit_path <- Sys.getenv("NEURALSBI_JUNIT")

if (nzchar(junit_path)) {
  # test_check() runs the files from tests/testthat, so a relative path would
  # land there rather than where the caller meant. Resolve it first.
  junit_path <- file.path(
    normalizePath(dirname(junit_path), mustWork = FALSE),
    basename(junit_path)
  )
  test_check(
    "neuralsbi",
    reporter = MultiReporter$new(list(
      JunitReporter$new(file = junit_path),
      CheckReporter$new()
    ))
  )
} else {
  test_check("neuralsbi")
}

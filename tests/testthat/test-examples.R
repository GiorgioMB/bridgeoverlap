test_that("the exact eight schools evidence is reproduced", {
  expect_equal(eight_schools_log_evidence(), -31.311347, tolerance = 1e-6)
})

test_that("the conjugate regression evidence matches numerical integration", {
  set.seed(50)
  n <- 6; X <- matrix(rnorm(n), n, 1); y <- 0.5 * X[, 1] + rnorm(n)
  v0 <- 2; a0 <- 2; b0 <- 1.5
  inner <- function(s2) vapply(s2, function(v) {
    # integrate beta out analytically given sigma^2: y ~ N(0, v (I + v0 X X'))
    S <- v * (diag(n) + v0 * tcrossprod(X))
    R <- chol(S); z <- backsolve(R, y, transpose = TRUE)
    exp(-0.5 * sum(z^2) - sum(log(diag(R))) - n / 2 * log(2 * pi) +
          a0 * log(b0) - lgamma(a0) - (a0 + 1) * log(v) - b0 / v)
  }, numeric(1))
  num <- log(integrate(inner, 0, Inf, rel.tol = 1e-10)$value)
  expect_equal(conjugate_log_evidence(X, y, v0, a0, b0), num, tolerance = 1e-7)
})

test_that("the shipped Stan programs exist", {
  for (nm in c("eight_schools_centered_ml", "eight_schools_noncentered_ml",
               "conjugate_regression_ml", "gaussian_mixture_ml")) {
    expect_true(file.exists(bridgeoverlap_stan_file(nm)))
  }
  expect_error(bridgeoverlap_stan_file("nope"), "Unknown")
})

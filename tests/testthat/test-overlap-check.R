a <- function(x, s = 1) dnorm(x, log = TRUE) - dnorm(x, 0, s, log = TRUE)

test_that("the overlap check is a fixed threshold, 10 for independent draws", {
  set.seed(1); s <- sim_gauss(S1 = 3000)
  d <- bridge_diagnostics(s$a_post, s$a_prop, iid = TRUE, eb = FALSE, jackknife = FALSE)
  expect_equal(d$n_ov_threshold, 10)
  expect_identical(d$flags[["overlap"]], !(d$n_ov >= 10))
})

test_that("near-perfect runs are flagged at every sample size and get the finite-sample interval", {
  set.seed(2)
  for (S in c(500, 5000)) {
    d <- bridge_diagnostics(a(rnorm(S), 1.02), a(rnorm(S, 0, 1.02), 1.02), iid = TRUE, jackknife = FALSE)
    expect_gt(d$I, 0.9)
    expect_true(d$flags[["overlap"]])
    expect_true(all(is.finite(d$eb)))
    expect_match(paste(d$actions, collapse = " "), "near-perfect")
  }
})

test_that("Markov chain output uses n_ov / tau with threshold 30", {
  set.seed(3); S <- 4000; rho <- 0.8
  x <- as.numeric(stats::arima.sim(list(ar = rho), S, sd = sqrt(1 - rho^2)))
  la <- function(t) dnorm(t, log = TRUE) - dnorm(t, 0.3, 1.3, log = TRUE)
  d <- bridge_diagnostics(la(x), la(rnorm(S, 0.3, 1.3)), iid = FALSE, jackknife = FALSE)
  expect_equal(d$n_ov_threshold, 30)
  expect_gt(d$tau, 1)
  expect_lt(d$n_ovP, d$n_ov)
  expect_identical(d$flags[["overlap"]], !(d$n_ovP >= 30))
})

test_that("Geyer's estimator recovers the autocorrelation time of an AR(1) chain", {
  set.seed(10)
  rho <- 0.7
  x <- as.numeric(arima.sim(list(ar = rho), n = 50000))
  expect_equal(iat(x), (1 + rho) / (1 - rho), tolerance = 0.12)
  expect_equal(iat(rnorm(20000)), 1, tolerance = 0.1)
})

test_that("multi-chain and within-chain estimators agree on well-mixed chains", {
  set.seed(11)
  x <- unlist(lapply(1:4, function(m) as.numeric(arima.sim(list(ar = 0.5), n = 5000))))
  t_w <- iat(x, rep(5000, 4), "within")
  t_m <- iat(x, rep(5000, 4), "multichain")
  t_a <- iat(x, rep(5000, 4), "concat_ar")
  expect_equal(t_w, 3, tolerance = 0.15)
  expect_equal(t_m, t_w, tolerance = 0.15)
  expect_equal(t_a, t_w, tolerance = 0.15)
  expect_error(iat(x, c(100, 100)), "sum")
})

test_that("chains that disagree inflate the multi-chain estimate", {
  set.seed(12)
  x <- c(rnorm(2000, 0), rnorm(2000, 3))
  expect_gt(iat(x, c(2000, 2000), "multichain"), 10 * iat(x, c(2000, 2000), "within"))
})

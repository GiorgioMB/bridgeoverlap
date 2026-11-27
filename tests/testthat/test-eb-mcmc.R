test_that("the MCMC finite-sample interval contains the estimate, widens as gamma decreases, and covers", {
  set.seed(31); hits <- 0; rho <- 0.8
  la <- function(t) dnorm(t, log = TRUE) - dnorm(t, 0.3, 1.3, log = TRUE)
  for (r in 1:30) {
    x <- as.numeric(stats::arima.sim(list(ar = rho), 3000, sd = sqrt(1 - rho^2)))
    ap <- la(x); aq <- la(rnorm(3000, 0.3, 1.3))
    int <- bridge_fs_interval_mcmc(ap, aq, gamma = 1 - rho)
    lz <- bridge_logml(ap, aq)$logml
    expect_true(int[["lower"]] <= lz && lz <= int[["upper"]])
    hits <- hits + (int[["lower"]] <= 0 && 0 <= int[["upper"]])
  }
  expect_gte(hits, 29)
  wider <- bridge_fs_interval_mcmc(ap, aq, gamma = (1 - rho) / 2)
  expect_gt(wider[["upper"]] - wider[["lower"]], int[["upper"]] - int[["lower"]])
  expect_error(bridge_fs_interval_mcmc(ap, aq, gamma = 0), "gamma")
})

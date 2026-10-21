test_that("the finite-sample interval contains the estimate and is wider than the delta-method interval", {
  set.seed(20)
  s <- sim_gauss(S1 = 3000)
  d <- bridge_diagnostics(s$a_post, s$a_prop, iid = TRUE, jackknife = FALSE, eb_method = "mp")
  expect_true(d$eb[["lower"]] < d$logml && d$logml < d$eb[["upper"]])
  expect_true(d$eb[["lower"]] < s$logZ && s$logZ < d$eb[["upper"]])
  ratio <- diff(d$eb) / (2 * qnorm(0.975) * d$mcse)
  # about 2.3 + 12 / sqrt(n_ov) times as wide for delta = 0.05
  expect_equal(unname(ratio), 2.3 + 12 / sqrt(d$n_ov), tolerance = 0.15)
  expect_equal(bridge_eb_interval(s$a_post, s$a_prop), d$eb, tolerance = 1e-6)
  # the default finite-sample interval (Bennett) is narrower but still wider than the delta interval
  d2 <- bridge_diagnostics(s$a_post, s$a_prop, iid = TRUE, jackknife = FALSE)
  expect_equal(bridge_fs_interval(s$a_post, s$a_prop), d2$eb, tolerance = 1e-6)
  expect_lt(diff(d2$eb), diff(d$eb)); expect_gt(diff(d2$eb), 2 * qnorm(0.975) * d2$mcse)
})

test_that("the finite-sample interval is unbounded when the samples barely overlap", {
  ap <- c(rep(Inf, 20), 0)   # one posterior draw in the overlap region
  aq <- c(rep(-Inf, 20), 0)  # one proposal draw in the overlap region
  int <- bridge_eb_interval(ap, aq)
  expect_equal(unname(int), c(-Inf, Inf))
})

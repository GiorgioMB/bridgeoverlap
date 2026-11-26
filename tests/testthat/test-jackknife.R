test_that("null quantiles are reproducible, sensible and leave the RNG untouched", {
  set.seed(30); r1 <- runif(1)
  set.seed(30); q <- jackknife_null_quantiles(10, probs = c(0.5, 0.99), nsim = 5000, seed = 99); r2 <- runif(1)
  expect_identical(r1, r2)
  expect_identical(q, jackknife_null_quantiles(10, probs = c(0.5, 0.99), nsim = 5000, seed = 99))
  expect_true(all(q["G_max", ] <= (10 - 1) / sqrt(10)))
  expect_true(all(diff(q["I_max", ]) > 0))
})

test_that("blocks respect chain boundaries and G is adjusted to a multiple of M", {
  b <- bridgeoverlap:::.block_ids_post(c(30, 30, 30), 6)
  expect_equal(as.vector(table(b)), rep(15, 6))
  expect_true(all(b[1:30] %in% 1:2) && all(b[31:60] %in% 3:4))
  set.seed(31)
  s <- sim_gauss(S1 = 900)
  expect_message(d <- bridge_diagnostics(s$a_post, s$a_prop, chains = 3, G = 20, nsim = 5000),
                 "using G = 21")
  expect_equal(d$settings$G, 21L)
})

test_that("the failure diagnostics are close to their null behaviour for a correct run", {
  set.seed(32)
  s <- sim_gauss(S1 = 4000, m = 0.2, s = 1.4)
  d <- bridge_diagnostics(s$a_post, s$a_prop, iid = TRUE, nsim = 5000)
  expect_true(is.finite(d$t_T) && abs(d$t_T) < 6)
  expect_true(d$V_jack > 0.1 && d$V_jack < 5)
  expect_false(isTRUE(d$flags[["overlap"]]))
})

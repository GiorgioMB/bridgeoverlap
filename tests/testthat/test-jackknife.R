test_that("blocks respect chain boundaries and G is adjusted to a multiple of M", {
  b <- bridgeoverlap:::.block_ids_post(c(30, 30, 30), 6)
  expect_equal(as.vector(table(b)), rep(15, 6))
  expect_true(all(b[1:30] %in% 1:2) && all(b[31:60] %in% 3:4))
  set.seed(31)
  s <- sim_gauss(S1 = 900)
  expect_message(d <- bridge_diagnostics(s$a_post, s$a_prop, chains = 3, G = 20),
                 "using G = 21")
  expect_equal(d$settings$G, 21L)
})

test_that("the failure diagnostics are close to their null behaviour for a correct run", {
  set.seed(32)
  s <- sim_gauss(S1 = 4000, m = 0.2, s = 1.4)
  d <- bridge_diagnostics(s$a_post, s$a_prop, iid = TRUE)
  expect_true(is.finite(d$V_jack) && d$V_jack < 3)
  expect_true(d$V_jack > 0.1 && d$V_jack < 5)
  expect_false(isTRUE(d$flags[["overlap"]]))
})

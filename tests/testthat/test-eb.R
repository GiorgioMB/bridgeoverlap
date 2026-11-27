test_that("the finite-sample interval contains the estimate, is wider than the delta-method interval and covers", {
  set.seed(20); hits <- 0
  for (r in 1:40) {
    s <- sim_gauss(S1 = 1500)
    fs <- bridge_fs_interval(s$a_post, s$a_prop)
    lz <- bridge_logml(s$a_post, s$a_prop)$logml
    expect_true(fs[["lower"]] <= lz && lz <= fs[["upper"]])
    hits <- hits + (fs[["lower"]] <= s$logZ && s$logZ <= fs[["upper"]])
  }
  expect_gte(hits, 39)
  d <- bridge_diagnostics(s$a_post, s$a_prop, iid = TRUE, jackknife = FALSE)
  expect_true(d$interval[["lower"]] < d$logml && d$logml < d$interval[["upper"]])
  expect_gt(d$interval[["upper"]] - d$interval[["lower"]], 2 * stats::qnorm(0.975) * d$mcse)
})

test_that("the vectorized radius solves Bennett's equation", {
  v <- c(0.5, 3, 40, 900); L <- log(2 / 0.04); x <- .bennett_radius(v, L)
  h <- function(u) (1 + u) * log1p(u) - u
  expect_lt(max(abs(v * h(x / v) - L)), 1e-8)
})

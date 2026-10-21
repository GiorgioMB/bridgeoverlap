# Tests for the plain rule and the finite-sample interval (Bennett).
a <- function(x, s = 1) dnorm(x, log = TRUE) - dnorm(x, 0, s, log = TRUE)

test_that("the plain rule uses the fixed threshold and flags near-perfect runs", {
  set.seed(1); x <- rnorm(2000); y <- rnorm(2000, 0, 1.02)
  d <- bridge_diagnostics(a(x, 1.02), a(y, 1.02), iid = TRUE, eb = FALSE, jackknife = FALSE)
  expect_equal(d$n_ov_threshold, 10)
  expect_true(d$I > 0.9)
  expect_true(d$flags[["overlap"]])                      # n_ov small near perfection
  expect_match(paste(d$actions, collapse = " "), "near-perfect")
})

test_that("the Bennett interval contains the estimate, is narrower than the old one, and covers", {
  set.seed(2); hits <- 0
  for (r in 1:40) {
    x <- rnorm(1500); y <- rnorm(1500, 1.5)
    ap <- dnorm(x, log = TRUE) - dnorm(x, 1.5, log = TRUE); aq <- dnorm(y, log = TRUE) - dnorm(y, 1.5, log = TRUE)
    fs <- bridge_fs_interval(ap, aq); eb <- bridge_eb_interval(ap, aq)
    lz <- bridge_logml(ap, aq)$logml
    expect_true(fs[["lower"]] <= lz && lz <= fs[["upper"]])
    expect_lt(fs[["upper"]] - fs[["lower"]], eb[["upper"]] - eb[["lower"]])
    hits <- hits + (fs[["lower"]] <= 0 && 0 <= fs[["upper"]])
  }
  expect_gte(hits, 39)
})

test_that("the vectorized radius solves Bennett's equation", {
  v <- c(0.5, 3, 40, 900); L <- log(2 / 0.04); x <- bridgeoverlap:::.bennett_radius(v, L)
  h <- function(u) (1 + u) * log1p(u) - u
  expect_lt(max(abs(v * h(x / v) - L)), 1e-8)
})

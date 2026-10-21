test_that("the near-perfect regime is flagged only when N H cannot be certified", {
  a <- function(t) dnorm(t, log = TRUE) - dnorm(t, 0, 1.02, log = TRUE)   # proposal close to the target
  set.seed(70)
  small <- bridge_diagnostics(a(rnorm(15)), a(rnorm(15, 0, 1.02)), iid = TRUE, eb = FALSE, jackknife = FALSE,
                              rule = "gate")
  expect_gt(small$I, 0.5)
  expect_lt(small$NH_lower, small$NH_required)
  expect_true(small$flags[["overlap"]])
  expect_match(small$overlap_reason, "near-perfect")
  expect_true(is.na(small$accuracy_bound))
  set.seed(71)
  large <- bridge_diagnostics(a(rnorm(2000)), a(rnorm(2000, 0, 1.02)), iid = TRUE, eb = FALSE, jackknife = FALSE,
                              rule = "gate")
  expect_gt(large$I, 0.5)
  expect_gte(large$NH_lower, large$NH_required)
  expect_false(large$flags[["overlap"]])
  expect_true(is.finite(large$accuracy_bound))
  expect_equal(large$NH_required, 16 * log(2 / 0.05))
})

test_that("the overlap flag under poor overlap is unchanged", {
  set.seed(72)
  a <- function(t) dnorm(t, log = TRUE) - dnorm(t, 4, 1, log = TRUE)       # mean shift, poor overlap
  d <- bridge_diagnostics(a(rnorm(100)), a(rnorm(100, 4)), iid = TRUE, eb = FALSE, jackknife = FALSE)
  expect_lte(d$I, 0.5)
  expect_equal(d$flags[["overlap"]], d$n_ov < 10)
  if (d$flags[["overlap"]]) expect_match(d$overlap_reason, "n_ov-hat")
})

test_that("the continuous rule flags the tiny near-perfect run and passes the large one", {
  a <- function(t) dnorm(t, log = TRUE) - dnorm(t, 0, 1.02, log = TRUE)
  set.seed(70)
  small <- bridge_diagnostics(a(rnorm(15)), a(rnorm(15, 0, 1.02)), iid = TRUE, eb = FALSE, jackknife = FALSE, rule = "continuous")
  expect_true(small$flags[["overlap"]])
  set.seed(71)
  large <- bridge_diagnostics(a(rnorm(2000)), a(rnorm(2000, 0, 1.02)), iid = TRUE, eb = FALSE, jackknife = FALSE, rule = "continuous")
  expect_false(large$flags[["overlap"]])
  expect_true(is.finite(large$accuracy_bound))
})

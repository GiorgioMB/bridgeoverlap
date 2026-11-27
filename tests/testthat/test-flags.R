a <- function(t, m = 0, s = 1) dnorm(t, log = TRUE) - dnorm(t, m, s, log = TRUE)

test_that("poor overlap is flagged with a rerun factor, good overlap passes", {
  set.seed(6)
  bad <- bridge_diagnostics(a(rnorm(60), 6), a(rnorm(60, 6), 6), iid = TRUE, jackknife = FALSE)
  expect_true(bad$flags[["overlap"]])
  expect_match(paste(bad$actions, collapse = " "), "poor overlap|do not overlap")
  if (is.finite(bad$draws_factor)) expect_gte(bad$draws_factor, 1)
  set.seed(7)
  good <- bridge_diagnostics(a(rnorm(4000), 0.5, 1.3), a(rnorm(4000, 0.5, 1.3), 0.5, 1.3), iid = TRUE, jackknife = FALSE)
  expect_false(good$flags[["overlap"]])
})

test_that("the MCSE never exceeds its deterministic bound (Proposition 4(a))", {
  set.seed(8)
  for (m in c(2, 4, 8)) {
    d <- bridge_diagnostics(a(rnorm(200), m), a(rnorm(200, m), m), iid = TRUE, eb = FALSE, jackknife = FALSE)
    if (is.finite(d$logml)) expect_lte(d$mcse^2, 2 * 200 / 199 * min(1, 1 / d$NH0) + 1e-10)
  }
})

# Tests for the continuous overlap rule, cross-fitting folds, the decision rule,
# the printout and the Bayes-factor combination.

test_that("the continuous rule is the gated rule under poor overlap and continuous in I-hat", {
  set.seed(3)
  a <- function(t, m) dnorm(t, log = TRUE) - dnorm(t, m, log = TRUE)
  x <- rnorm(400); y <- rnorm(400, 3)
  dc <- bridge_diagnostics(a(x, 3), a(y, 3), iid = TRUE, jackknife = FALSE, rule = "continuous")
  dg <- bridge_diagnostics(a(x, 3), a(y, 3), iid = TRUE, jackknife = FALSE, rule = "gate")
  expect_lt(dc$I, 0.5)
  expect_equal(dc$n_ov_threshold, 10)
  expect_identical(dc$flags[["overlap"]], dg$flags[["overlap"]])
  # threshold as a function of I-hat: 10 * min(1, 2 (1 - I))
  thr <- function(I) 10 * max(0, min(1, 2 * (1 - I)))
  expect_equal(thr(0.5 - 1e-9), thr(0.5 + 1e-9), tolerance = 1e-6)
  y2 <- rnorm(400, 0.2, 1.05)
  d2 <- bridge_diagnostics(rule = "continuous", a(x, 0) - dnorm(x, 0.2, 1.05, log = TRUE) + dnorm(x, 0, log = TRUE),
                           dnorm(y2, log = TRUE) - dnorm(y2, 0.2, 1.05, log = TRUE), iid = TRUE, jackknife = FALSE)
  expect_gt(d2$I, 0.5)
  expect_equal(d2$n_ov_threshold, thr(d2$I))
  expect_true(is.finite(d2$NH_lower))
})

test_that("finite-sample quantities are reported in every regime", {
  set.seed(4)
  x <- rnorm(30); y <- rnorm(30, 0.1, 1.02)
  a <- function(t) dnorm(t, log = TRUE) - dnorm(t, 0.1, 1.02, log = TRUE)
  d <- bridge_diagnostics(a(x), a(y), iid = TRUE, jackknife = FALSE)
  expect_true(all(is.finite(d$eb)))
  expect_true(is.finite(d$NH_lower))
  expect_true(is.na(d$accuracy_bound) || d$accuracy_bound > 0)
})

test_that("folds: within-fold variances, fold-aware blocks and EB radii", {
  set.seed(5)
  S <- 800
  # two folds with different proposals g1 = N(0.3, 1.2^2), g2 = N(-0.2, 1.1^2), target N(0, 1)
  x1 <- rnorm(S / 2); x2 <- rnorm(S / 2)
  y1 <- rnorm(S / 2, 0.3, 1.2); y2 <- rnorm(S / 2, -0.2, 1.1)
  l1 <- function(t) dnorm(t, log = TRUE) - dnorm(t, 0.3, 1.2, log = TRUE)
  l2 <- function(t) dnorm(t, log = TRUE) - dnorm(t, -0.2, 1.1, log = TRUE)
  ap <- c(l1(x1), l2(x2)); aq <- c(l1(y1), l2(y2))
  folds <- list(post = rep(1:2, each = S / 2), prop = rep(1:2, each = S / 2))
  d <- bridge_diagnostics(ap, aq, chains = 2, folds = folds, iid = TRUE, G = 10)
  dp <- bridge_diagnostics(ap, aq, chains = 2, iid = TRUE, G = 10)
  expect_equal(d$logml, dp$logml)                 # same root: pooled log ratios
  lr <- 0
  pp <- plogis(d$logml - ap - lr); pq <- plogis(d$logml - aq - lr)
  n_within <- sum(sapply(1:2, function(k) S / 2 * var(pq[folds$prop == k]) + S / 2 * var(pp[folds$post == k])))
  expect_equal(d$n_ov, n_within)
  expect_lte(d$n_ov, dp$n_ov + 1e-8)               # pooled variances include the between-fold term
  expect_equal(d$settings$K, 2L)
  expect_true(d$eb[["lower"]] < d$logml && d$logml < d$eb[["upper"]])
  expect_gt(d$eb[["upper"]] - d$eb[["lower"]], dp$eb[["upper"]] - dp$eb[["lower"]])  # delta / K per fold
  expect_true(all(is.finite(d$jackknife$c_b)))
  expect_error(bridge_diagnostics(ap, aq, folds = list(post = rep(1:2, each = S / 2), prop = rep(2:1, each = S / 2)), iid = TRUE),
               "contiguous")
})

test_that("decision rule maps flags to actions", {
  set.seed(6)
  a <- function(t) dnorm(t, log = TRUE) - dnorm(t, 6, log = TRUE)
  d <- bridge_diagnostics(a(rnorm(60)), a(rnorm(60, 6)), iid = TRUE, jackknife = FALSE)
  expect_true(d$flags[["overlap"]])
  expect_match(d$decision, "do not report")
  expect_true(any(grepl("rerun from scratch|overlap", d$actions)))
  if (is.finite(d$draws_factor)) expect_gte(d$draws_factor, 1)
  rep <- bridge_report(d)
  expect_true(all(c("quantity", "value", "rule", "status") %in% names(rep)))
  expect_output(print(d), "Decision")
})

test_that("Bayes factors combine two runs", {
  set.seed(7)
  a <- function(t, lz) lz + dnorm(t, log = TRUE) - dnorm(t, 0.5, 1.3, log = TRUE)
  d1 <- bridge_diagnostics(a(rnorm(2000), 1), a(rnorm(2000, 0.5, 1.3), 1), iid = TRUE, nsim = 2e4)
  d2 <- bridge_diagnostics(a(rnorm(2000), -0.5), a(rnorm(2000, 0.5, 1.3), -0.5), iid = TRUE, nsim = 2e4)
  bf <- bridge_bayes_factor(d1, d2)
  expect_equal(bf$log_bf, d1$logml - d2$logml)
  expect_equal(bf$mcse, sqrt(d1$mcse^2 + d2$mcse^2))
  expect_true(bf$eb[["lower"]] <= bf$log_bf && bf$log_bf <= bf$eb[["upper"]])
  expect_output(print(bf), "log Bayes factor")
})

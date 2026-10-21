test_that("the root solver agrees with uniroot and W is non-increasing", {
  set.seed(1)
  s <- sim_gauss()
  r <- bridge_logml(s$a_post, s$a_prop)
  expect_true(r$converged)
  u <- uniroot(function(cc) bridge_estimating_function(cc, s$a_post, s$a_prop),
               c(-20, 20), tol = 1e-14)$root
  expect_equal(r$logml, u, tolerance = 1e-10)
  W <- bridge_estimating_function(seq(-5, 8, by = 0.25), s$a_post, s$a_prop)
  expect_true(all(diff(W) <= 0))
})

test_that("the estimate converges from far-away starting values", {
  set.seed(2)
  s <- sim_gauss()
  ref <- bridge_logml(s$a_post, s$a_prop)$logml
  for (c0 in c(-200, 200)) {
    expect_equal(bridge_logml(s$a_post, s$a_prop, c0 = c0)$logml, ref, tolerance = 1e-10)
  }
})

test_that("the MCSE identity and the contraction factor hold at the root", {
  set.seed(3)
  s <- sim_gauss(S1 = 1500, S2 = 1000)
  d <- bridge_diagnostics(s$a_post, s$a_prop, iid = TRUE, jackknife = FALSE, eb = FALSE)
  S1 <- 1500; S2 <- 1000; N <- S1 + S2; s1 <- S1 / N; s2 <- S2 / N
  ch <- d$logml; Z <- exp(ch); lr <- log(s1 / s2)
  pp <- plogis(ch - s$a_post - lr); pq <- plogis(ch - s$a_prop - lr)
  expect_equal(sum(pp), sum(1 - pq), tolerance = 1e-8)
  # delta-method MCSE of Fruhwirth-Schnatter / Micaletto & Vehtari (2025)
  Ni <- exp(s$a_prop) / (s1 * exp(s$a_prop) + s2 * Z)
  Dj <- 1 / (s1 * exp(s$a_post) + s2 * Z)
  re2 <- var(Ni) / (S2 * mean(Ni)^2) + var(Dj) / (S1 * mean(Dj)^2)
  expect_equal(d$mcse, sqrt(re2), tolerance = 1e-8)
  expect_equal(d$n_ov, (d$NH0 * d$mcse)^2, tolerance = 1e-8)
  # contraction factor T = phi'(c-hat), phi(x) = log T(exp(x))
  Tmap <- function(z) mean(exp(s$a_prop) / (s1 * exp(s$a_prop) + s2 * z)) /
    mean(1 / (s1 * exp(s$a_post) + s2 * z))
  h <- 1e-5
  dphi <- (log(Tmap(exp(ch + h))) - log(Tmap(exp(ch - h)))) / (2 * h)
  expect_equal(d$T, dphi, tolerance = 1e-6)
})

test_that("the two-point family gives Delta = log(n_prop / n_post) and the conventions at +-Inf", {
  eps <- 0.1; logZ <- 0.4; S <- 60
  lf <- function(t) ifelse(t == 2, -Inf, logZ + log(c(1 - eps, eps, 1)[t + 1]))
  lg <- function(t) ifelse(t == 0, -Inf, log(c(1, eps, 1 - eps)[t + 1]))
  for (seed in 1:10) {
    set.seed(seed)
    xs <- sample(0:2, S, TRUE, prob = c(1 - eps, eps, 0))
    ys <- sample(0:2, S, TRUE, prob = c(0, eps, 1 - eps))
    est <- bridge_logml(lf(xs) - lg(xs), lf(ys) - lg(ys))$logml
    np <- sum(ys == 1); nq <- sum(xs == 1)
    expected <- if (np == 0) -Inf else if (nq == 0) Inf else logZ + log(np / nq)
    expect_equal(est, expected, tolerance = 1e-9)
  }
  expect_equal(bridge_logml(rep(Inf, 5), c(0, -Inf))$logml, Inf)
  expect_equal(bridge_logml(c(0, Inf), rep(-Inf, 5))$logml, -Inf)
})

test_that("diagnostics are invariant to a common shift of the log ratios", {
  set.seed(4)
  s <- sim_gauss(S1 = 800)
  d1 <- bridge_diagnostics(s$a_post, s$a_prop, iid = TRUE, nsim = 5000)
  d2 <- bridge_diagnostics(s$a_post + 50, s$a_prop + 50, iid = TRUE, nsim = 5000)
  expect_equal(d2$logml - d1$logml, 50, tolerance = 1e-8)
  for (k in c("mcse", "n_ov", "I", "NH0", "T", "t_T", "V_jack", "I_max", "G_max")) {
    expect_equal(d2[[k]], d1[[k]], tolerance = 1e-6)
  }
  expect_equal(unname(d2$eb - d1$eb), c(50, 50), tolerance = 1e-6)
})

test_that("NA log ratios are replaced with a warning", {
  set.seed(5)
  s <- sim_gauss(S1 = 200)
  ap <- s$a_post; ap[3] <- NA
  expect_warning(bridge_logml(ap, s$a_prop), "NA/NaN")
})

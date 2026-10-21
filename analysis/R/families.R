# Families F1-F3 of Section 7.1 ------------------------------------------------------
#
# In each family the log ratio log p - log g depends on the draw only through a
# scalar u: u = |theta|^2 for the radial families F1 and F3 and u = theta for
# F2. Sampling u directly is therefore exact and much cheaper than sampling
# theta. All targets are normalized densities multiplied by exp(logZ).
#
#   F1  p = N(0, I_d),     g = N(0, sigma^2 I_d)   (scale mismatch)
#   F2  p = N(0, 1),       g = N(mu_g, 1)          (mean shift)
#   F3  p = t_3(0, I_d),   g = N(0, 3 I_d)         (heavy-tailed target; g/f has
#                                                   infinite variance under g)

family_F1 <- function(d, sigma) {
  list(
    id = "F1", label = sprintf("F1 (d=%d, sigma=%g)", d, sigma), d = d,
    lr = function(u) d * log(sigma) - u / 2 * (1 - 1 / sigma^2),
    rpost = function(n) stats::rchisq(n, d),
    rprop = function(n) sigma^2 * stats::rchisq(n, d),
    dpost = function(u) stats::dchisq(u, d),
    dprop = function(u) stats::dchisq(u / sigma^2, d) / sigma^2,
    qpost = function(p) stats::qchisq(p, d),
    qprop = function(p) sigma^2 * stats::qchisq(p, d),
    lower = 0
  )
}

family_F2 <- function(mu) {
  list(
    id = "F2", label = sprintf("F2 (mu_g=%g)", mu), d = 1L, mu = mu,
    lr = function(u) -mu * u + mu^2 / 2,
    rpost = function(n) stats::rnorm(n),
    rprop = function(n) stats::rnorm(n, mu),
    dpost = function(u) stats::dnorm(u),
    dprop = function(u) stats::dnorm(u, mu),
    qpost = function(p) stats::qnorm(p),
    qprop = function(p) mu + stats::qnorm(p),
    lower = -Inf
  )
}

family_F3 <- function(d, nu = 3, var_g = 3) {
  const <- lgamma((nu + d) / 2) - lgamma(nu / 2) - d / 2 * log(nu * pi) +
    d / 2 * log(2 * pi * var_g)
  list(
    id = "F3", label = sprintf("F3 (d=%d)", d), d = d,
    lr = function(u) const - (nu + d) / 2 * log1p(u / nu) + u / (2 * var_g),
    rpost = function(n) d * stats::rf(n, d, nu),
    rprop = function(n) var_g * stats::rchisq(n, d),
    dpost = function(u) stats::df(u / d, d, nu) / d,
    dprop = function(u) stats::dchisq(u / var_g, d) / var_g,
    qpost = function(p) d * stats::qf(p, d, nu),
    qprop = function(p) var_g * stats::qchisq(p, d),
    lower = 0
  )
}

families_universality <- function() {
  list(family_F1(10, 1.5), family_F1(40, 1.3),
       family_F2(2), family_F2(4), family_F2(6),
       family_F3(20), family_F3(100))
}

# Adaptive quadrature split at break points (sharp features are resolved).
integrate_pieces <- function(f, breaks, lower) {
  b <- sort(unique(c(lower, breaks[is.finite(breaks) & breaks > lower], Inf)))
  sum(vapply(seq_len(length(b) - 1L), function(i) {
    stats::integrate(f, b[i], b[i + 1L], rel.tol = 1e-11, abs.tol = 0,
                     subdivisions = 2000L, stop.on.error = FALSE)$value
  }, numeric(1)))
}

# Population overlap functionals (Definition 1) for s1 = S1/N.
overlap_population <- function(fam, s1 = 0.5) {
  s2 <- 1 - s1
  lr0 <- log(s1 / s2)
  pistar <- function(u) stats::plogis(-fam$lr(u) - lr0)
  probs <- c(1e-12, 1e-6, 1e-3, 0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.98, 1 - 1e-3, 1 - 1e-6, 1 - 1e-12)
  half <- tryCatch(stats::uniroot(function(u) fam$lr(u) + lr0,
                                  range(c(fam$qpost(probs), fam$qprop(probs))))$root,
                   error = function(e) NULL)
  breaks <- c(fam$qpost(probs), fam$qprop(probs), half)
  E <- function(dens, h) integrate_pieces(function(u) h(u) * dens(u), breaks, fam$lower)
  Ep1 <- E(fam$dpost, pistar)
  Ep2 <- E(fam$dpost, function(u) pistar(u)^2)
  Eg1 <- E(fam$dprop, pistar)
  Eg2 <- E(fam$dprop, function(u) pistar(u)^2)
  H <- s1 * Ep1
  I <- H / (s1 * s2)
  list(H = H, H_check = s2 * (1 - Eg1), I = I, beta = 1 - I,
       var_p = Ep2 - Ep1^2, var_g = Eg2 - Eg1^2, s1 = s1, pistar = pistar)
}

# Effective overlap and sigma_N for S1 = S2 = S.
overlap_at <- function(pop, S) {
  N <- 2 * S
  list(n_ov = N * pop$H * pop$beta,
       n_ov_check = S * (pop$var_p + pop$var_g),  # Lemma 2(e)
       sigma_N = sqrt(pop$beta / (N * pop$H)))
}

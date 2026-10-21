# Eight schools (Section 7.5) ------------------------------------------------------------
#
# y_j ~ N(theta_j, sigma_j^2), theta_j ~ N(mu, tau^2), mu ~ N(0, 5^2),
# tau ~ half-Cauchy(0, 5). Unconstrained spaces: (theta, mu, log tau) for the
# centred and (eta, mu, log tau), theta = mu + tau eta, for the non-centred
# parameterization. log Z = -31.3113 (bridgeoverlap::eight_schools_log_evidence()).

es <- local({
  d <- eight_schools_data()
  list(y = d$y, sigma = d$sigma, J = d$J, logZ = eight_schools_log_evidence(d))
})
es_names <- c(paste0("theta", 1:8), "mu", "log_tau")

# log p(y | tau) with theta and mu integrated out (Sherman-Morrison), vectorized in tau.
es_log_marg_tau <- function(tau) {
  vapply(tau, function(t) {
    dj <- es$sigma^2 + t^2
    s <- sum(1 / dj); tt <- sum(es$y / dj); q <- sum(es$y^2 / dj)
    -0.5 * (8 * log(2 * pi) + sum(log(dj)) + log1p(25 * s) + q - 25 * tt^2 / (1 + 25 * s))
  }, numeric(1))
}

# Inverse-CDF table for the marginal posterior of log tau on a fine grid.
es_tau_table <- local({
  tab <- NULL
  function() {
    if (is.null(tab)) {
      lt <- seq(log(1e-8), log(1e4), length.out = 200001)
      tau <- exp(lt)
      ld <- es_log_marg_tau(tau) + stats::dcauchy(tau, 0, 5, log = TRUE) + log(2) + lt
      dens <- exp(ld - max(ld))
      cdf <- c(0, cumsum((dens[-1] + dens[-length(dens)]) / 2 * diff(lt)))
      tab <<- list(lt = lt, cdf = cdf / cdf[length(cdf)])
    }
    tab
  }
})

# Unnormalized log density of log tau given y (half-Cauchy(0, 5) prior, Jacobian of log).
es_log_dens_lt <- function(lt) es_log_marg_tau(exp(lt)) + stats::dcauchy(exp(lt), 0, 5, log = TRUE) + log(2) + lt

# Exact draws of log tau by rejection from a Student-t_4 envelope. The target decays
# exponentially in both tails and the envelope polynomially, so the ratio vanishes
# at both ends; its maximum is located on a fine grid and inflated by 5%.
es_lt_rejection <- local({
  env <- NULL
  function(n) {
    if (is.null(env)) {
      lt <- seq(log(1e-8), log(1e4), length.out = 100001)
      ld <- es_log_dens_lt(lt)
      m <- lt[which.max(ld)]
      w <- exp(ld - max(ld)); s <- sqrt(sum(w * (lt - m)^2) / sum(w))
      le <- stats::dt((lt - m) / s, 4, log = TRUE) - log(s)
      env <<- list(m = m, s = s, logM = max(ld - le) + log(1.05))
    }
    out <- numeric(0)
    while (length(out) < n) {
      k <- ceiling(1.6 * (n - length(out))) + 10L
      x <- env$m + env$s * stats::rt(k, 4)
      la <- es_log_dens_lt(x) - (stats::dt((x - env$m) / env$s, 4, log = TRUE) - log(env$s)) - env$logM
      if (any(la > 0)) stop("rejection envelope violated")
      out <- c(out, x[log(stats::runif(k)) < la])
    }
    out[seq_len(n)]
  }
})

# Exact posterior draws in the centred unconstrained space (vectorized). log tau is
# drawn by inverting the CDF on a grid of 200001 points ("grid")
# or by rejection ("rejection"); theta and mu are then drawn from their Gaussian
# conditionals.
es_exact_draws <- function(n, lt_method = c("grid", "rejection")) {
  lt_method <- match.arg(lt_method)
  if (lt_method == "grid") {
    tab <- es_tau_table()
    keep <- !duplicated(tab$cdf)
    lt <- stats::approx(tab$cdf[keep], tab$lt[keep], xout = stats::runif(n))$y
  } else {
    lt <- es_lt_rejection(n)
  }
  tau <- exp(lt)
  dj <- outer(tau^2, es$sigma^2, "+")
  prec_mu <- 1 / 25 + rowSums(1 / dj)
  mu <- stats::rnorm(n, as.vector((1 / dj) %*% es$y) / prec_mu, 1 / sqrt(prec_mu))
  prec_th <- outer(1 / tau^2, 1 / es$sigma^2, "+")
  mean_th <- (matrix(es$y / es$sigma^2, n, 8, byrow = TRUE) + mu / tau^2) / prec_th
  th <- matrix(stats::rnorm(n * 8, mean_th, 1 / sqrt(prec_th)), n, 8)
  out <- cbind(th, mu, lt)
  colnames(out) <- es_names
  out
}

es_to_noncentered <- function(U) {
  tau <- exp(U[, 10])
  out <- U
  out[, 1:8] <- (U[, 1:8] - U[, 9]) / tau
  colnames(out) <- c(paste0("eta", 1:8), "mu", "log_tau")
  out
}

# Unconstrained log densities including the Jacobian of log tau and every constant.
es_lp_centered <- function(u, data = NULL) {
  th <- u[1:8]; mu <- u[9]; lt <- u[10]; tau <- exp(lt)
  sum(stats::dnorm(es$y, th, es$sigma, log = TRUE)) +
    sum(stats::dnorm(th, mu, tau, log = TRUE)) + stats::dnorm(mu, 0, 5, log = TRUE) +
    log(2) + stats::dcauchy(tau, 0, 5, log = TRUE) + lt
}

es_lp_noncentered <- function(u, data = NULL) {
  eta <- u[1:8]; mu <- u[9]; lt <- u[10]; tau <- exp(lt)
  sum(stats::dnorm(es$y, mu + tau * eta, es$sigma, log = TRUE)) +
    sum(stats::dnorm(eta, log = TRUE)) + stats::dnorm(mu, 0, 5, log = TRUE) +
    log(2) + stats::dcauchy(tau, 0, 5, log = TRUE) + lt
}

# Row-wise version for the random-walk Metropolis sampler.
es_lp_centered_rows <- function(U) {
  tau <- exp(U[, 10])
  rowSums(stats::dnorm(matrix(es$y, nrow(U), 8, byrow = TRUE), U[, 1:8],
                       matrix(es$sigma, nrow(U), 8, byrow = TRUE), log = TRUE)) +
    rowSums(stats::dnorm(U[, 1:8], U[, 9], tau, log = TRUE)) +
    stats::dnorm(U[, 9], 0, 5, log = TRUE) + log(2) +
    stats::dcauchy(tau, 0, 5, log = TRUE) + U[, 10]
}

# Random-walk Metropolis on the centred parameterization: M chains of length n
# started at exact draws, Gaussian increments with covariance 2.38^2 Sigma / 10.
es_rwm <- function(n, M, Sigma) {
  L <- t(chol(2.38^2 * Sigma / 10))
  X <- es_exact_draws(M)
  lp <- es_lp_centered_rows(X)
  draws <- array(NA_real_, c(n, M, 10))
  acc <- 0
  for (t in seq_len(n)) {
    Y <- X + matrix(stats::rnorm(M * 10), M, 10) %*% t(L)
    lpy <- es_lp_centered_rows(Y)
    ok <- log(stats::runif(M)) < lpy - lp
    ok[is.na(ok)] <- FALSE
    X[ok, ] <- Y[ok, ]
    lp[ok] <- lpy[ok]
    acc <- acc + sum(ok)
    draws[t, , ] <- X
  }
  list(chains = lapply(seq_len(M), function(m) {
    out <- draws[, m, ]
    colnames(out) <- es_names
    out
  }), acceptance = acc / (n * M))
}

es_bounds <- function(names) {
  list(lb = stats::setNames(rep(-Inf, length(names)), names),
       ub = stats::setNames(rep(Inf, length(names)), names))
}

# Gaussian proposal fitted by moments, its draws and log density (as bridgesampling's
# method = "normal" on unconstrained parameters).
es_gfit <- function(X) list(m = colMeans(X), V = stats::cov(X))
es_gdraw <- function(n, fit) {
  L <- t(chol(fit$V)); d <- length(fit$m)
  sweep(matrix(stats::rnorm(n * d), n, d) %*% t(L), 2, fit$m, "+")
}
es_glog <- function(X, fit) {
  d <- length(fit$m); R <- chol(fit$V)
  Z <- backsolve(R, t(X) - fit$m, transpose = TRUE)
  -0.5 * colSums(Z^2) - sum(log(diag(R))) - d / 2 * log(2 * pi)
}

# Normal mixture with label switching (posteriordb normal_2-normal_mixture) -------------
#
# y_n ~ theta N(mu_1, 1) + (1 - theta) N(mu_2, 1), theta ~ U(0, 1), mu_k ~ N(0, 10^2).
# The posterior is invariant under (theta, mu_1, mu_2) -> (1 - theta, mu_2, mu_1); with the
# posteriordb data (components at -10 and 10) the two labellings are separated by a region
# of negligible mass, so log Z = log 2 + log of the integral over one labelling, computed by
# the trapezoid rule on a Laplace-standardized grid in (logit theta, mu_1, mu_2).
mixture_data <- function() {
  f <- file.path(BO_ROOT, "data", "normal_2.json.zip")
  ex <- tempfile(); dir.create(ex)
  js <- utils::unzip(f, exdir = ex)
  jsonlite::read_json(js[grepl("json$", js)][1], simplifyVector = TRUE)
}
# log f on the unconstrained scale u = (logit theta, mu_1, mu_2), with the Jacobian; vectorized over rows.
mixture_lf <- function(U, y) {
  th <- stats::plogis(U[, 1])
  ll <- vapply(seq_len(nrow(U)), function(i) {
    a <- log(th[i]) + stats::dnorm(y, U[i, 2], 1, log = TRUE)
    b <- log1p(-th[i]) + stats::dnorm(y, U[i, 3], 1, log = TRUE)
    m <- pmax(a, b); sum(m + log(exp(a - m) + exp(b - m)))
  }, numeric(1))
  ll + stats::dnorm(U[, 2], 0, 10, log = TRUE) + stats::dnorm(U[, 3], 0, 10, log = TRUE) + log(th) + log1p(-th)
}
mixture_log_evidence <- function(y, h = 0.25, zmax = 9) {
  f <- function(u) -mixture_lf(matrix(u, 1), y)
  start <- c(stats::qlogis(mean(y < 0)), mean(y[y < 0]), mean(y[y > 0]))
  opt <- stats::optim(start, f, method = "BFGS", hessian = TRUE)
  Lc <- t(chol(solve(opt$hessian)))
  zs <- seq(-zmax, zmax, by = h)
  Zg <- as.matrix(expand.grid(zs, zs, zs))
  lv <- mixture_lf(sweep(Zg %*% t(Lc), 2, opt$par, "+"), y)
  m <- max(lv)
  list(logZ = log(2) + m + log(sum(exp(lv - m))) + 3 * log(h) + sum(log(diag(Lc))),
       edge = max(lv[apply(abs(Zg) >= zmax - 1e-9, 1, any)]) - m, mode = opt$par)
}
# One-component model y ~ N(mu, 1), mu ~ N(0, 10^2): closed form.
single_log_evidence <- function(y, v0 = 100) {
  n <- length(y)
  -0.5 * (n * log(2 * pi) + log1p(v0 * n) + sum(y^2) - v0 * sum(y)^2 / (1 + v0 * n))
}

# Frühwirth-Schnatter (2004): a proposal invariant under relabelling. The posterior draws are
# relabelled to mu_1 < mu_2, a Gaussian is fitted to the first half of every chain, and the
# proposal is the equal mixture of it and its image under the permutation
# P(u) = (-u_1, u_3, u_2) (|det P| = 1). f and g are both invariant, so the responsibilities,
# and hence the estimator, do not depend on which labelling each chain explored.
mixture_symmetric_bridge <- function(fit, y) {
  ud <- posterior::as_draws_array(fit$unconstrain_draws(format = "draws_array"))
  n_iter <- posterior::niterations(ud); M <- posterior::nchains(ud)
  arr <- unclass(ud)
  P <- function(U) cbind(-U[, 1], U[, 3], U[, 2])
  relab <- function(U) { s <- U[, 2] > U[, 3]; U[s, ] <- P(U[s, , drop = FALSE]); U }
  h <- n_iter %/% 2L
  fitp <- lapply(seq_len(M), function(m) relab(matrix(arr[seq_len(h), m, ], h)))
  bridge <- do.call(rbind, lapply(seq_len(M), function(m) matrix(arr[h + seq_len(n_iter - h), m, ], n_iter - h)))
  g <- es_gfit(do.call(rbind, fitp))
  lg <- function(U) { a <- es_glog(U, g); b <- es_glog(P(U), g); m <- pmax(a, b); m + log(0.5 * exp(a - m) + 0.5 * exp(b - m)) }
  Y <- es_gdraw(nrow(bridge), g)
  flip <- stats::runif(nrow(Y)) < 0.5
  Y[flip, ] <- P(Y[flip, , drop = FALSE])
  lf <- function(U) vapply(seq_len(nrow(U)), function(i) fit$log_prob(U[i, ], jacobian = TRUE), numeric(1))
  list(a_post = lf(bridge) - lg(bridge), a_prop = lf(Y) - lg(Y), chains = rep(seq_len(M), each = n_iter - h))
}

# Radon models of posteriordb with exact evidence ------------------------------------
#
# y = F gamma + Z_j eta_j + eps within county j, eta_j ~ N(0, D), D = diag(sigma_r^2),
# eps ~ N(0, sigma_y^2 I), gamma ~ N(0, 10^2 I), sigma_y, sigma_r ~ half-N(0, 1). The
# county effects and gamma are integrated out in closed form (Woodbury identities on
# per-county sufficient statistics), which leaves an integral over the 1 + q log-scales
# u = log(sigma). It is evaluated by the trapezoid rule on a grid in Laplace-standardized
# coordinates z = L^{-1}(u - u_hat); the integrand is analytic and decays at least
# exponentially, so the rule converges geometrically. radon_check_logml() halves the step
# and widens the grid to confirm the digits.
#
# Models (posteriordb names; Stan programs in inst/stan/*_ml.stan):
#   pooled  F = [1, floor], no county effects                      (radon_pooled)
#   vi      F = [1, floor], Z = [1]                                 (radon_variable_intercept_*)
#   hi      F = [1, log_uppm, floor], Z = [1]                       (radon_hierarchical_intercept_*)
#   vis     F = [1, floor], Z = [1, floor]                          (radon_variable_intercept_slope_*)

radon_data <- function(which = c("mn", "all")) {
  which <- match.arg(which)
  f <- file.path(BO_ROOT, "data", sprintf("radon_%s.json.zip", which))
  ex <- tempfile(); dir.create(ex)
  js <- utils::unzip(f, exdir = ex)
  jsonlite::read_json(js[grepl("json$", js)][1], simplifyVector = TRUE)
}

.radon_design <- function(dat, model) {
  one <- rep(1, dat$N)
  Fm <- switch(model, pooled = cbind(one, dat$floor_measure), vi = cbind(one, dat$floor_measure),
               hi = cbind(one, dat$log_uppm, dat$floor_measure), vis = cbind(one, dat$floor_measure))
  Zm <- switch(model, pooled = NULL, vi = cbind(one), hi = cbind(one), vis = cbind(one, dat$floor_measure))
  list(F = Fm, Z = Zm, y = dat$log_radon, g = dat$county_idx, J = dat$J, n = dat$N)
}

# Per-county sufficient statistics.
.radon_stats <- function(des) {
  g <- factor(des$g, levels = seq_len(des$J))
  sums <- function(M) rowsum(M, g, reorder = TRUE)
  p <- ncol(des$F); q <- if (is.null(des$Z)) 0L else ncol(des$Z)
  st <- list(p = p, q = q, n = des$n, nj = as.numeric(table(g)), yy = as.numeric(sums(cbind(des$y^2))),
             Fy = sums(des$F * des$y), FF = sums(do.call(cbind, lapply(seq_len(p), function(a) des$F[, a] * des$F))))
  if (q > 0) {
    st$Zy <- sums(des$Z * des$y)
    st$ZZ <- sums(do.call(cbind, lapply(seq_len(q), function(a) des$Z[, a] * des$Z)))
    st$ZF <- sums(do.call(cbind, lapply(seq_len(q), function(a) des$Z[, a] * des$F)))   # columns (a, b): Z_a F_b
  }
  st
}

# log p(y | sigma_y, sigma_r) for a matrix SG of scales (rows: points; columns sigma_y,
# sigma_r...). County effects are removed by per-county Woodbury identities, gamma by a
# closed-form p x p Gaussian integral (p <= 3), all vectorized over points.
.radon_loglik <- function(SG, st, V0 = 100) {
  p <- st$p; q <- st$q; J <- length(st$nj); G <- nrow(SG)
  s2 <- SG[, 1]^2
  L <- log(s2) %o% st$nj                                     # G x J: n_j log sigma_y^2
  Q <- (1 / s2) %o% st$yy
  A <- array(0, c(G, p, p)); b <- matrix(0, G, p)
  for (a in seq_len(p)) {
    b[, a] <- (1 / s2) * sum(st$Fy[, a])
    for (c2 in seq_len(p)) A[, a, c2] <- (1 / s2) * sum(st$FF[, (a - 1) * p + c2])
  }
  if (q > 0) {
    d2 <- SG[, -1, drop = FALSE]^2
    ZZ <- function(a, c2) matrix(st$ZZ[, (a - 1) * q + c2], G, J, byrow = TRUE)
    ZF <- function(a, c2) matrix(st$ZF[, (a - 1) * p + c2], G, J, byrow = TRUE)   # Z_a F_c2
    Zy <- function(a) matrix(st$Zy[, a], G, J, byrow = TRUE)
    # Sigma_j = s2 I + Z_j D Z_j'; with K_j = (D^{-1} + Z_j'Z_j / s2)^{-1} / s2^2,
    # Sigma_j^{-1} = I / s2 - Z_j K_j Z_j' and log det Sigma_j = n_j log s2 + log det(I + D Z_j'Z_j / s2).
    if (q == 1) {
      dd <- d2[, 1]
      m <- 1 + dd * ZZ(1, 1) / s2                            # det(I + D Z'Z / s2)
      k11 <- dd / (s2^2 * m)                                 # K
      L <- L + log(m)
      zy <- Zy(1)
      Q <- Q - k11 * zy^2
      for (a in seq_len(p)) {
        za <- ZF(1, a)
        b[, a] <- b[, a] - rowSums(k11 * za * zy)
        for (c2 in seq_len(p)) A[, a, c2] <- A[, a, c2] - rowSums(k11 * za * ZF(1, c2))
      }
    } else {
      d1 <- d2[, 1]; dB <- d2[, 2]
      z11 <- ZZ(1, 1); z22 <- ZZ(2, 2); z12 <- ZZ(1, 2)
      # B = s2 I + D Z'Z; K = B^{-1} D / s2 (symmetric) and det(I + D Z'Z / s2) = det(B) / s2^2,
      # both finite at sigma_r = 0, where the grid has its half-weight end points.
      detB <- (s2 + d1 * z11) * (s2 + dB * z22) - d1 * dB * z12^2
      L <- L + log(detB) - 2 * log(s2)
      k11 <- d1 * (s2 + dB * z22) / (s2 * detB); k22 <- dB * (s2 + d1 * z11) / (s2 * detB)
      k12 <- -d1 * dB * z12 / (s2 * detB)
      zy1 <- Zy(1); zy2 <- Zy(2)
      Q <- Q - (k11 * zy1^2 + 2 * k12 * zy1 * zy2 + k22 * zy2^2)
      for (a in seq_len(p)) {
        za1 <- ZF(1, a); za2 <- ZF(2, a)
        b[, a] <- b[, a] - rowSums(k11 * za1 * zy1 + k12 * (za1 * zy2 + za2 * zy1) + k22 * za2 * zy2)
        for (c2 in seq_len(p)) {
          zc1 <- ZF(1, c2); zc2 <- ZF(2, c2)
          A[, a, c2] <- A[, a, c2] - rowSums(k11 * za1 * zc1 + k12 * (za1 * zc2 + za2 * zc1) + k22 * za2 * zc2)
        }
      }
    }
  }
  # gamma ~ N(0, V0 I_p): P = A + I / V0, contribution -1/2 {p log V0 + log det P - b' P^{-1} b}
  P <- A; for (a in seq_len(p)) P[, a, a] <- P[, a, a] + 1 / V0
  if (p == 2) {
    det <- P[, 1, 1] * P[, 2, 2] - P[, 1, 2]^2
    quad <- (P[, 2, 2] * b[, 1]^2 - 2 * P[, 1, 2] * b[, 1] * b[, 2] + P[, 1, 1] * b[, 2]^2) / det
  } else if (p == 3) {
    p11 <- P[, 1, 1]; p22 <- P[, 2, 2]; p33 <- P[, 3, 3]; p12 <- P[, 1, 2]; p13 <- P[, 1, 3]; p23 <- P[, 2, 3]
    C11 <- p22 * p33 - p23^2; C12 <- -(p12 * p33 - p13 * p23); C13 <- p12 * p23 - p13 * p22
    C22 <- p11 * p33 - p13^2; C23 <- -(p11 * p23 - p12 * p13); C33 <- p11 * p22 - p12^2
    det <- p11 * C11 + p12 * C12 + p13 * C13
    quad <- (C11 * b[, 1]^2 + C22 * b[, 2]^2 + C33 * b[, 3]^2 + 2 * C12 * b[, 1] * b[, 2] +
               2 * C13 * b[, 1] * b[, 3] + 2 * C23 * b[, 2] * b[, 3]) / det
  } else stop("p must be 2 or 3")
  -0.5 * (st$n * log(2 * pi) + rowSums(L) + rowSums(Q) + p * log(V0) + log(det) - quad)
}

# log of the integrand on the linear scale: log p(y | sigma) + sum_k log{2 phi(sigma_k)}.
.radon_logint <- function(SG, st) .radon_loglik(SG, st) + rowSums(log(2) + stats::dnorm(SG, log = TRUE))

# Trapezoid rule on a tensor grid in (sigma_y, sigma_r) on the linear scale. The integrand
# is even in every sigma_r (it depends on sigma_r^2), so the rule on [0, b] with half weight
# at 0 is spectrally accurate; sigma_y is integrated around its mode. Ranges are set where
# the log integrand has dropped by `drop` along each axis through the mode.
radon_log_evidence <- function(model, dat, n_y = 61L, n_r = 81L, drop = 40, chunk = 20000L) {
  st <- .radon_stats(.radon_design(dat, model))
  k <- 1L + st$q
  f <- function(s) -.radon_logint(matrix(s, 1), st)
  opt <- stats::optim(rep(0.5, k), f, method = "L-BFGS-B", lower = rep(1e-6, k), upper = rep(5, k))
  top <- -opt$value
  edge_of <- function(j, dir) {             # distance along axis j where the log integrand drops by `drop`
    g <- function(t) { s <- opt$par; s[j] <- s[j] + dir * t; if (s[j] <= 0) return(drop + 1); top - (-f(s)) - drop }
    hi <- 0.01; while (g(hi) < 0 && hi < 20) hi <- 2 * hi
    if (g(hi) < 0) return(hi)
    stats::uniroot(g, c(0, hi))$root
  }
  axes <- lapply(seq_len(k), function(j) {
    up <- opt$par[j] + 1.3 * edge_of(j, +1)
    if (j == 1L) seq(max(1e-6, opt$par[1] - 1.3 * edge_of(1, -1)), up, length.out = n_y) else seq(0, up, length.out = n_r)
  })
  w <- lapply(axes, function(x) { h <- diff(x)[1]; ww <- rep(h, length(x)); ww[c(1, length(x))] <- h / 2; ww })
  grid <- as.matrix(expand.grid(axes)); lw <- log(Reduce(`%o%`, w))
  lv <- numeric(nrow(grid))
  for (s in seq(1L, nrow(grid), by = chunk)) {
    i <- s:min(nrow(grid), s + chunk - 1L)
    lv[i] <- .radon_logint(grid[i, , drop = FALSE], st)
  }
  lv[!is.finite(lv)] <- -Inf
  m <- max(lv)
  far <- apply(grid, 1, function(g) any(g >= vapply(axes, max, 0) - 1e-12) || g[1] <= min(axes[[1]]) + 1e-12)
  list(logZ = m + log(sum(exp(lv + as.vector(lw) - m))), mode = opt$par, edge = max(lv[far]) - m, k = k,
       points = nrow(grid))
}

radon_check_logml <- function(model, dat) {
  a <- radon_log_evidence(model, dat)
  b <- radon_log_evidence(model, dat, n_y = 91L, n_r = 121L, drop = 50)
  c(logZ = b$logZ, diff_step = abs(a$logZ - b$logZ), edge = max(a$edge, b$edge))
}

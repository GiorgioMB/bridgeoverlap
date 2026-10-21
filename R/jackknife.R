# Delete-a-group jackknife -----------------------------------------
#
# Both samples are partitioned into G contiguous blocks aligned with the
# sampling order; posterior blocks never straddle chain boundaries (G is a
# multiple of the number of chains M and each chain holds G/M blocks). Deleting
# block b from both samples leaves s1 and s2 unchanged, so the leave-block-out
# estimate is the root of W - W_b with the same responsibility function.

.split_sizes <- function(n, k) diff(round(seq(0, n, length.out = k + 1L)))

.block_ids_post <- function(chain_lengths, G) {
  per <- G %/% length(chain_lengths)
  unlist(lapply(seq_along(chain_lengths), function(m) {
    rep((m - 1L) * per + seq_len(per), times = .split_sizes(chain_lengths[m], per))
  }), use.names = FALSE)
}

.block_ids_prop <- function(S2, G) rep(seq_len(G), times = .split_sizes(S2, G))

# Cross-fitting: G / K blocks per fold. Posterior blocks are formed within the
# chains of each fold (chains are nested in folds and ordered fold by fold), and
# the proposal draws of fold k are split into the G / K blocks numbered like the
# posterior blocks of that fold, so that deleting a block removes posterior and
# proposal draws of one fold in proportion and leaves s1, s2 unchanged.
.block_ids_prop_folds <- function(fold_prop, G) {
  K <- max(fold_prop)
  per <- G %/% K
  out <- integer(length(fold_prop))
  for (k in seq_len(K)) {
    idx <- which(fold_prop == k)
    out[idx] <- (k - 1L) * per + rep(seq_len(per), times = .split_sizes(length(idx), per))
  }
  out
}

# Contraction factor T = 1 - H-hat / H0-hat at the root c.
.calibration_stat <- function(a_post, a_prop, c, lr) {
  NH0 <- sum(.resp(a_post, c, lr))
  NH <- sum(.resp_var(a_post, c, lr)) + sum(.resp_var(a_prop, c, lr))
  1 - NH / NH0
}

.jackknife <- function(a_post, a_prop, lr, chat, block_post, block_prop, G,
                       tol, maxit) {
  c_b <- T_b <- rep(NA_real_, G)
  for (b in seq_len(G)) {
    ap <- a_post[block_post != b]
    aq <- a_prop[block_prop != b]
    r <- .bridge_root(ap, aq, lr, c0 = chat, tol = tol, maxit = maxit)
    c_b[b] <- r$c
    if (is.finite(r$c)) T_b[b] <- .calibration_stat(ap, aq, r$c, lr)
  }
  list(c_b = c_b, T_b = T_b)
}

.jack_cache <- new.env(parent = emptyenv())

#' Null quantiles of the jackknife influence statistics
#'
#' Simulates the pivotal limits, for a fixed number of blocks `G`, of the
#' maximal standardized block influence \eqn{\mathcal I_{\max}} and of the
#' Grubbs-type statistic \eqn{\mathcal G_{\max}} under the fixed-`G` regime, in which the
#' numbers of draws grow with `G` fixed. With
#' \eqn{Z_1,\dots,Z_G} i.i.d. standard normal and
#' \eqn{D = \max_b |Z_b - \bar Z|}, the limits are
#' \eqn{\mathcal I_{\max} = D\sqrt G/(G-1)} (the standardized block effects have
#' variance \eqn{1/(G-1)}, not \eqn{1/G}) and
#' \eqn{\mathcal G_{\max} = D / s_Z}, the Grubbs statistic. Results are cached
#' for the session, and the caller's random number stream is left untouched.
#'
#' @param G Number of blocks.
#' @param probs Probabilities at which to return quantiles.
#' @param nsim Number of simulated replicates.
#' @param seed Seed used for the simulation.
#' @return A matrix with rows `I_max` and `G_max` and one column per element of
#'   `probs`.
#' @examples
#' jackknife_null_quantiles(20, probs = 1 - 0.05 / 4, nsim = 2e4)
#' @export
jackknife_null_quantiles <- function(G, probs = c(0.95, 1 - 0.05 / 4),
                                     nsim = 1e5L, seed = 20260917L) {
  G <- as.integer(G)
  if (G < 3L) stop("`G` must be at least 3.", call. = FALSE)
  key <- paste(G, nsim, seed, sep = ":")
  sims <- .jack_cache[[key]]
  if (is.null(sims)) {
    sims <- .with_seed(seed, {
      I_max <- G_max <- numeric(nsim)
      done <- 0L
      while (done < nsim) {
        n <- min(20000L, nsim - done)
        Z <- matrix(stats::rnorm(n * G), n, G)
        Zc <- Z - rowMeans(Z)
        A <- abs(Zc)
        D <- A[cbind(seq_len(n), max.col(A, ties.method = "first"))]
        s <- sqrt(rowSums(Zc^2) / (G - 1L))
        idx <- done + seq_len(n)
        I_max[idx] <- D * sqrt(G) / (G - 1L)
        G_max[idx] <- D / s
        done <- done + n
      }
      list(I_max = I_max, G_max = G_max)
    })
    assign(key, sims, envir = .jack_cache)
  }
  out <- rbind(I_max = stats::quantile(sims$I_max, probs, names = FALSE),
               G_max = stats::quantile(sims$G_max, probs, names = FALSE))
  colnames(out) <- format(probs, digits = 6)
  out
}

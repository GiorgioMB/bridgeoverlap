# Integration with bridgesampling ---------------------------------------------------
#
# A `bridge` object from bridgesampling::bridge_sampler(method = "normal",
# repetitions = 1) stores
#   q11 = log f at the posterior draws used in the iterative scheme,
#   q12 = log g at the same draws,
#   q21 = log f at the proposal draws,  q22 = log g at the proposal draws,
# so the log ratios needed by Algorithm 1 are a_post = q11 - q12 and
# a_prop = q21 - q22. No further evaluations are required.

.extract_log_ratios <- function(x) {
  if (!inherits(x, "bridge")) {
    stop("Expected a `bridge` object from bridgesampling::bridge_sampler().", call. = FALSE)
  }
  if (!identical(x$method, "normal")) {
    stop("Only method = \"normal\" is supported: warp3 objects do not store the ",
         "warp matrix needed to reconstruct the log ratios.", call. = FALSE)
  }
  if (is.null(x$q11) || is.null(x$q21)) {
    stop("The bridge object does not contain q11/q21 (use repetitions = 1).", call. = FALSE)
  }
  list(a_post = x$q11 - x$q12, a_prop = x$q21 - x$q22)
}

#' @rdname bridge_diagnostics
#' @export
bridge_diagnostics.bridge <- function(x, chains = attr(x, "bridgeoverlap_chains"),
                                      iid = isTRUE(attr(x, "bridgeoverlap_iid")),
                                      folds = attr(x, "bridgeoverlap_folds"), ...) {
  a <- .extract_log_ratios(x)
  out <- bridge_diagnostics.default(a$a_post, a$a_prop, chains = chains, iid = iid,
                                    folds = folds, ...)
  out$bridgesampling <- list(
    logml = x$logml, niter = x$niter,
    mcse_logml = if (is.null(x$mcse_logml)) NA_real_ else x$mcse_logml,
    difference = x$logml - out$logml
  )
  out
}

#' @rdname bridge_diagnostics
#' @export
bridge_diagnostics.bridge_list <- function(x, ...) {
  stop("Diagnostics need a single bridge run: call bridge_sampler() with repetitions = 1.",
       call. = FALSE)
}

#' Chain identifiers of the posterior draws used in the bridge
#'
#' [bridgesampling::bridge_sampler()] splits the posterior draws into a part
#' that fits the proposal and a part that enters the iterative scheme. Which
#' draws end up in the second part, and hence their chain structure, depends
#' on the method: the `stanfit`, `stanreg` and `mcmc.list` methods use the
#' second half of *each* chain (`split = "per_chain"`), whereas the `matrix`
#' method, and the `CmdStanMCMC` method of \pkg{bridgesampling} 1.2-1 that
#' calls it, use the second half of the *stacked* draws (`split = "stacked"`).
#'
#' @param n_iter Number of post-warmup iterations per chain.
#' @param n_chains Number of chains.
#' @param split `"per_chain"` or `"stacked"`, see Details.
#' @return Integer vector of chain ids, one per posterior draw in the bridge,
#'   to be passed as `chains` to [bridge_diagnostics()].
#' @export
#' @examples
#' bridge_chain_ids(1000, 4, "per_chain")   # 4 chains of 500 draws
#' table(bridge_chain_ids(1000, 4, "stacked"))  # chains 3 and 4, 1000 draws each
bridge_chain_ids <- function(n_iter, n_chains, split = c("per_chain", "stacked")) {
  split <- match.arg(split)
  if (split == "per_chain") {
    keep <- !(seq_len(n_iter) %in% seq_len(round(n_iter / 2)))
    rep(seq_len(n_chains), each = sum(keep))
  } else {
    ids <- rep(seq_len(n_chains), each = n_iter)
    nr <- length(ids)
    ids[!(seq_len(nr) %in% seq_len(round(nr / 2)))]
  }
}

#' Two-fold cross-fitted bridge sampling
#'
#' Splits the posterior draws into two groups A and B. With an even number of
#' chains (the default, `fold = "auto"`), A and B are the first and the second
#' half of the chains; otherwise, or with `fold = "halves"`, they are the first
#' and second halves of every chain. Fold 1 bridges B with a proposal
#' \eqn{g_A} fitted on A, fold 2 bridges A with \eqn{g_B} fitted on B (two calls
#' of [bridgesampling::bridge_sampler()]), and the two estimating functions are
#' added. Folding by chains keeps the draws of each fold independent of the
#' draws its proposal is fitted on, which makes every fold's estimating function
#' exactly unbiased at \eqn{\log Z}; halves of the same
#' chain are dependent, so folding by halves is exact only for independent
#' draws. With chains stuck in different modes, folding by chains also makes the
#' two proposals miss each other's draws, which the overlap certificate flags. Each fold draws as many proposal draws
#' as it has posterior draws (the default of \pkg{bridgesampling}), so the
#' allocation is proportional, \eqn{s_1 = s_2 = 1/2} in each fold and overall,
#' and every fold's estimating function has mean zero at \eqn{\log Z}
#' conditionally on the draws its proposal was fitted on. Responsibilities
#' depend on a draw only through its log ratio, so the cross-fitted estimating
#' function is that of the pooled log ratios and its root is the ordinary
#' estimate computed from them.
#'
#' The returned object records the chains (each half-chain is a chain of its
#' fold) and the folds, so that [bridge_diagnostics()] computes the variances,
#' the autocorrelation times, the jackknife blocks and the finite-sample bounds
#' within folds. Use a number of blocks `G` that is a
#' multiple of the number of sequences, e.g. `G = 16` for four chains.
#'
#' Fitting the proposal on the same draws that enter the bridge biases
#' \eqn{\log\hat Z} at second order, by an amount that grows with the number of
#' fitted parameters; splitting removes the bias but
#' leaves half of the posterior draws out of the bridge. Cross-fitting removes
#' the leading bias term while using every draw in both roles.
#'
#' @param samples An `mcmc.list` (one element per chain) or a matrix (a single
#'   chain) of posterior draws, as accepted by
#'   [bridgesampling::bridge_sampler()].
#' @param fold `"auto"` (by chains when their number is even, otherwise by
#'   halves of every chain), `"chains"` or `"halves"`.
#' @param ... Further arguments for [bridgesampling::bridge_sampler()], e.g.
#'   `log_posterior`, `data`, `lb`, `ub`.
#' @param use_neff Passed to [bridgesampling::bridge_sampler()]; `FALSE` (the
#'   default) gives the optimal bridge estimator.
#' @param silent Passed to [bridgesampling::bridge_sampler()].
#' @return A `bridge` object with pooled `q11`, `q12`, `q21`, `q22`, the
#'   cross-fitted `logml`, and the chain structure in the attribute
#'   `"bridgeoverlap_chains"`. The two underlying runs are kept in the attribute
#'   `"bridgeoverlap_runs"`.
#' @export
bridge_crossfit <- function(samples, ..., fold = c("auto", "chains", "halves"),
                            use_neff = FALSE, silent = TRUE) {
  fold <- match.arg(fold)
  chains <- if (inherits(samples, "mcmc.list")) {
    lapply(samples, as.matrix)
  } else if (is.matrix(samples)) {
    list(samples)
  } else {
    stop("`samples` must be an mcmc.list or a matrix.", call. = FALSE)
  }
  M <- length(chains)
  if (fold == "auto") fold <- if (M >= 2L && M %% 2L == 0L) "chains" else "halves"
  if (fold == "chains" && (M < 2L || M %% 2L != 0L)) {
    stop("Folding by chains needs an even number of chains (at least two).", call. = FALSE)
  }
  as_input <- function(first, second) {
    mats <- Map(rbind, first, second)
    if (length(mats) == 1L) mats[[1L]] else coda::mcmc.list(lapply(mats, coda::mcmc))
  }
  if (fold == "chains") {
    # Fold 1 bridges the chains of group B with the proposal fitted on group A, fold 2 the
    # reverse. Pairing chain i of one group with chain i of the other makes the per-chain
    # split of bridgesampling (first half fits, second half bridges) do exactly this.
    n <- min(vapply(chains, nrow, integer(1)))
    if (n < 4L) stop("Chains are too short for cross-fitting.", call. = FALSE)
    chains <- lapply(chains, function(m) m[seq_len(n), , drop = FALSE])
    A <- chains[seq_len(M / 2L)]; B <- chains[M / 2L + seq_len(M / 2L)]
    h <- n
  } else {
    if (M > 1L && !silent) message("Folding by halves of chains: exact only for independent draws.")
    n <- min(vapply(chains, nrow, integer(1)))
    n <- n - n %% 2L  # equal halves
    if (n < 8L) stop("Chains are too short for cross-fitting.", call. = FALSE)
    h <- n / 2L
    A <- lapply(chains, function(m) m[seq_len(h), , drop = FALSE])
    B <- lapply(chains, function(m) m[h + seq_len(h), , drop = FALSE])
  }
  run1 <- bridgesampling::bridge_sampler(as_input(A, B), ..., use_neff = use_neff,
                                         silent = silent, repetitions = 1)
  run2 <- bridgesampling::bridge_sampler(as_input(B, A), ..., use_neff = use_neff,
                                         silent = silent, repetitions = 1)
  out <- list(
    q11 = c(run1$q11, run2$q11), q12 = c(run1$q12, run2$q12),
    q21 = c(run1$q21, run2$q21), q22 = c(run1$q22, run2$q22),
    niter = c(run1$niter, run2$niter), method = "normal"
  )
  out$logml <- bridge_logml(out$q11 - out$q12, out$q21 - out$q22)$logml
  class(out) <- "bridge"
  Mf <- length(A)                     # sequences per fold
  attr(out, "bridgeoverlap_chains") <- c(rep(seq_len(Mf), each = h), Mf + rep(seq_len(Mf), each = h))
  attr(out, "bridgeoverlap_folds") <- list(post = rep(1:2, c(length(run1$q11), length(run2$q11))),
                                           prop = rep(1:2, c(length(run1$q21), length(run2$q21))))
  attr(out, "bridgeoverlap_fold_by") <- fold
  attr(out, "bridgeoverlap_runs") <- list(run1, run2)
  out
}

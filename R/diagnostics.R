#' Bridge sampling with diagnostics
#'
#' Computes the optimal bridge sampling estimate of \eqn{\log Z}, its delta-method MCSE
#' (with a within-chain long-run variance for Markov chain output) and two checks of the
#' MCSE: the overlap check, \eqn{\hat n_{ov} \ge 10} for independent draws and
#' \eqn{\hat n^P_{ov} \ge 30} for Markov chain output, and the variance check, which
#' compares the delete-a-group jackknife variance with the squared MCSE. It also returns a
#' finite-sample confidence interval, valid at every sample size: for independent draws
#' always, and for Markov chain output when a lower bound `gamma` on the spectral gap of
#' the sampler is given. No additional evaluations of the posterior or the proposal are
#' needed.
#'
#' @param x Either an object of class `bridge` returned by
#'   [bridgesampling::bridge_sampler()], or a numeric vector of log ratios
#'   \eqn{\log\{f/g\}} at the posterior draws.
#' @param a_prop Numeric vector of log ratios at the proposal draws (only when `x` is
#'   numeric).
#' @param chains Chain structure of the posterior draws: `NULL` (a single sequence), the
#'   number of chains of equal length, or a vector of chain ids of length `length(x)`.
#' @param iid Logical: are the posterior draws independent? If `FALSE`, the MCSE uses a
#'   within-chain estimate of the integrated autocorrelation time.
#' @param gamma For Markov chain output, an optional lower bound on the spectral gap of
#'   the sampler; if given, the finite-sample interval is that of
#'   [bridge_fs_interval_mcmc()].
#' @param G Number of jackknife blocks; a multiple of the number of chains and folds.
#' @param alpha Level of the variance check.
#' @param delta Error probability of the finite-sample interval.
#' @param tau_method Estimator of the integrated autocorrelation time of the posterior
#'   responsibilities; see [iat()].
#' @param thresholds List with the overlap thresholds `n_ov` (default 10, independent
#'   draws) and `n_ovP` (default 30, Markov chain output).
#' @param folds For cross-fitted output: a list with integer vectors `post` and `prop`
#'   giving the fold of every draw (see [bridge_crossfit()]).
#' @param interval Logical: compute the finite-sample interval?
#' @param jackknife Logical: compute the variance check?
#' @param tol,maxit Tolerance and iteration limit of the root solver.
#' @param keep_responsibilities Logical: store the fitted responsibilities?
#' @param ... Passed to methods.
#' @return An object of class `bridge_diagnostics`: a list with the estimate `logml`, the
#'   `mcse`, the estimated effective overlaps `n_ov` and `n_ovP`, the autocorrelation time
#'   `tau`, the jackknife variance ratio `V_jack`, the finite-sample `interval`, the
#'   logical `flags` (`overlap`, `variance`), the `decision`, the `actions` and further
#'   quantities used by the analysis scripts.
#' @details
#' * **Overlap check.** The MCSE is reliable when the estimated effective overlap is
#'   large. Below the threshold, report the finite-sample interval. Under poor overlap, a
#'   fresh run with more draws raises \eqn{\hat n_{ov}} in proportion; near a perfect
#'   proposal it does not, and the finite-sample interval is the answer.
#' * **Variance check.** \eqn{\mathcal V_{jack}}, the jackknife variance over \eqn{G}
#'   blocks divided by the squared MCSE, is approximately \eqn{\chi^2_{G-1}/(G-1)} when the
#'   MCSE is right. With blocks aligned with the chains, misestimated autocorrelation, a
#'   dominating block and chains that disagree all make it large.
#'
#' The estimate uses \eqn{s_1 = S_1/N}. It coincides with the `logml` of
#' [bridgesampling::bridge_sampler()] when that is called with `use_neff = FALSE`; with the
#' default `use_neff = TRUE`, \pkg{bridgesampling} replaces \eqn{S_1} by a median effective
#' sample size inside the bridge function, which changes the estimator slightly.
#' @examples
#' set.seed(1)
#' S <- 2000
#' x <- rnorm(S); y <- rnorm(S, 0.5, 1.4)            # posterior and proposal draws
#' a <- function(t) dnorm(t, log = TRUE) - dnorm(t, 0.5, 1.4, log = TRUE)  # log Z = 0
#' d <- bridge_diagnostics(a(x), a(y), iid = TRUE)
#' d
#' bridge_report(d)
#' @export
bridge_diagnostics <- function(x, ...) UseMethod("bridge_diagnostics")

#' @rdname bridge_diagnostics
#' @export
bridge_diagnostics.default <- function(x, a_prop, chains = NULL, iid = FALSE, gamma = NULL,
                                       G = 20L, alpha = 0.05, delta = 0.05,
                                       tau_method = c("within", "multichain",
                                                      "concat_geyer", "concat_ar"),
                                       thresholds = list(), folds = NULL,
                                       interval = TRUE, jackknife = TRUE,
                                       tol = 1e-12, maxit = 500L,
                                       keep_responsibilities = FALSE, ...) {
  tau_method <- match.arg(tau_method)
  .check_scalar(alpha, "alpha", 0, 1)
  .check_scalar(delta, "delta", 0, 1)
  if (!is.null(gamma) && (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) ||
                          gamma <= 0 || gamma > 1))
    stop("`gamma` must be NULL or a number in (0, 1].", call. = FALSE)
  thr <- utils::modifyList(list(n_ov = 10, n_ovP = 30), thresholds)
  a_post <- .clean_log_ratios(x, "posterior")
  a_prop <- .clean_log_ratios(a_prop, "proposal")
  S1 <- length(a_post); S2 <- length(a_prop)
  if (S1 < 2L || S2 < 2L) stop("Both samples need at least two draws.", call. = FALSE)
  N <- S1 + S2; s1 <- S1 / N; s2 <- S2 / N; lr <- log(S1 / S2)
  # without a chain structure, every fold of cross-fitted output is its own sequence
  if (is.null(chains) && !is.null(folds) && length(folds$post) == S1) chains <- folds$post
  chain_lengths <- .chain_lengths(chains, S1)
  M <- length(chain_lengths)
  chain_id <- rep.int(seq_len(M), chain_lengths)
  fd <- .check_folds(folds, S1, S2, chain_id)
  K <- fd$K
  G <- if (jackknife) .resolve_G(G, M, chain_lengths, S2, K) else NA_integer_
  root <- .bridge_root(a_post, a_prop, lr, tol = tol, maxit = maxit)
  chat <- root$c
  n0 <- if (iid) thr$n_ov else thr$n_ovP
  settings <- list(S1 = S1, S2 = S2, N = N, s1 = s1, s2 = s2, M = M, chain_lengths = chain_lengths,
                   G = G, K = K, alpha = alpha, delta = delta, iid = iid, gamma = gamma,
                   n0 = n0, interval = interval, jackknife = jackknife,
                   tau_method = if (iid) "iid" else tau_method)
  if (!is.finite(chat)) return(.no_overlap_result(chat, root, settings))

  # Responsibilities and overlap, pooled over folds.
  pp <- .resp(a_post, chat, lr); pq <- .resp(a_prop, chat, lr)
  NH0 <- sum(pp)
  I_hat <- (sum(.resp_var(a_post, chat, lr)) + sum(.resp_var(a_prop, chat, lr))) / (N * s1 * s2)

  # Estimated effective overlap and MCSE, with variances within folds.
  fs <- do.call(rbind, lapply(seq_len(K), function(k) {
    ip <- fd$post == k; iq <- fd$prop == k
    tau <- if (iid) 1 else iat(pp[ip], rle(chain_id[ip])$lengths, method = tau_method)
    c(S1 = sum(ip), S2 = sum(iq), V_post = stats::var(pp[ip]), V_prop = stats::var(pq[iq]), tau = tau)
  }))
  prop_term <- sum(fs[, "S2"] * fs[, "V_prop"])
  post_term <- sum(fs[, "S1"] * fs[, "V_post"])
  n_ov <- prop_term + post_term
  n_P <- prop_term + sum(fs[, "S1"] * fs[, "V_post"] * fs[, "tau"])
  tau <- if (iid) 1 else if (post_term > 0) (n_P - prop_term) / post_term else NA_real_
  mcse <- sqrt(n_P) / NH0
  n_ovP <- if (n_P > 0) n_ov^2 / n_P else NA_real_

  # Finite-sample interval: independent draws, or Markov chain output with a gap bound.
  int <- c(lower = NA_real_, upper = NA_real_)
  if (interval && iid) {
    int <- .fs_interval(a_post, a_prop, lr, chat, delta,
                        fold_post = if (K > 1L) fd$post, fold_prop = if (K > 1L) fd$prop)
  } else if (interval && !is.null(gamma) && K == 1L) {
    int <- .hull_interval(function(cc) .fs_F_mcmc(cc, a_post, a_prop, lr, delta, gamma),
                          chat, a_post, a_prop, lr, 150L)
  }

  # Variance check: delete-a-group jackknife, with blocks within chains and folds.
  V_jack <- NA_real_; c_b <- NULL
  if (jackknife) {
    bpost <- .block_ids_post(chain_lengths, G)
    bprop <- if (K == 1L) .block_ids_prop(S2, G) else .block_ids_prop_folds(fd$prop, G)
    c_b <- .jackknife(a_post, a_prop, lr, chat, bpost, bprop, G, tol, maxit)
    if (all(is.finite(c_b))) V_jack <- (G - 1) / G * sum((c_b - mean(c_b))^2) / mcse^2
  }
  v0 <- if (jackknife) stats::qchisq(1 - alpha, df = G - 1) / (G - 1) else NA_real_

  stat <- if (iid) n_ov else n_ovP
  flags <- c(overlap = !isTRUE(stat >= n0),
             variance = if (is.finite(V_jack)) V_jack > v0 else NA)
  dec <- .decide(flags, stat = stat, n0 = n0, I_hat = I_hat, iid = iid,
                 have_interval = all(is.finite(int)))
  out <- list(
    logml = chat, mcse = mcse, n_ov = n_ov, n_ovP = n_ovP, tau = tau, V_jack = V_jack,
    interval = int, flags = flags, decision = dec$decision, actions = dec$actions,
    draws_factor = dec$draws_factor, critical = c(n0 = n0, V_jack = v0),
    I = I_hat, NH0 = NH0, n_P = n_P, mcse_iid = sqrt(n_ov) / NH0,
    V_post = if (K == 1L) fs[1, "V_post"] else post_term / S1,
    V_prop = if (K == 1L) fs[1, "V_prop"] else prop_term / S2,
    folds = if (K > 1L) fs else NULL, jackknife = if (jackknife) list(c_b = c_b) else NULL,
    solver = root[c("iterations", "converged", "status")], settings = settings
  )
  if (keep_responsibilities) out$responsibilities <- list(posterior = pp, proposal = pq)
  class(out) <- "bridge_diagnostics"
  out
}

# Validate the fold structure of cross-fitted output. Folds are integer labels
# 1..K, contiguous in both samples, with chains nested in folds and the proposal
# draws allocated in proportion to the posterior draws (S2k / S1k = S2 / S1),
# which keeps every fold's estimating function unbiased at log Z.
.check_folds <- function(folds, S1, S2, chain_id) {
  if (is.null(folds)) return(list(K = 1L, post = rep(1L, S1), prop = rep(1L, S2)))
  if (!is.list(folds) || is.null(folds$post) || is.null(folds$prop)) {
    stop("`folds` must be a list with integer vectors `post` and `prop`.", call. = FALSE)
  }
  fp <- as.integer(folds$post); fq <- as.integer(folds$prop)
  if (length(fp) != S1 || length(fq) != S2) stop("`folds` vectors must have lengths S1 and S2.", call. = FALSE)
  K <- max(fp)
  if (!setequal(unique(fp), seq_len(K)) || !setequal(unique(fq), seq_len(K))) {
    stop("Fold labels must be 1, ..., K in both samples.", call. = FALSE)
  }
  if (any(diff(fp) < 0) || any(diff(fq) < 0)) stop("Folds must be contiguous and ordered.", call. = FALSE)
  if (any(tapply(fp, chain_id, function(v) length(unique(v))) > 1L)) {
    stop("Every chain must lie in a single fold.", call. = FALSE)
  }
  ratio <- tabulate(fq, K) / tabulate(fp, K)
  if (any(abs(ratio / (S2 / S1) - 1) > 0.02)) {
    warning("Proposal draws are not allocated in proportion to the posterior draws of each fold ",
            "(S2k / S1k should equal S2 / S1); the cross-fitted estimating function is then biased.",
            call. = FALSE)
  }
  list(K = K, post = fp, prop = fq)
}

# Decision rule: report log Z-hat +/- z MCSE when both checks pass; otherwise map each
# flag to an action. Under poor overlap n_ov grows linearly in N at a fixed proposal, so
# the factor of more draws needed is n0 / n_ov; near a perfect proposal it does not grow.
.decide <- function(flags, stat, n0, I_hat, iid, have_interval) {
  factor <- if (isTRUE(flags[["overlap"]]) && is.finite(stat) && stat > 0) max(1, n0 / stat) else NA_real_
  use <- if (have_interval) "report the finite-sample interval" else if (iid)
    "the finite-sample interval is unbounded" else
    "report the finite-sample interval given a lower bound on the spectral gap (argument gamma)"
  act <- character()
  if (isTRUE(flags[["overlap"]])) {
    act <- c(act, if (isTRUE(I_hat > 0.5)) {
      paste0("near-perfect proposal: the MCSE cannot be checked here; ", use, ", which is short in this regime")
    } else if (is.finite(factor)) {
      sprintf("poor overlap: %s, or rerun from scratch with both samples at least %s times larger (do not extend this run), or improve the proposal",
              use, .fmt(factor, 2))
    } else "the samples (almost) do not overlap: improve the proposal before adding draws")
  }
  if (isTRUE(flags[["variance"]])) {
    act <- c(act, paste("the MCSE is too small (the jackknife and long-run variances disagree):",
                        "run longer chains, estimate the long-run variance within chains, and look for",
                        "a dominating block or chains in different modes"))
  }
  list(decision = if (length(act)) "do not report log Z-hat +/- z MCSE" else "report log Z-hat +/- z MCSE",
       actions = act, draws_factor = factor)
}

.chain_lengths <- function(chains, S1) {
  if (is.null(chains)) return(S1)
  if (length(chains) == 1L && is.numeric(chains)) {
    M <- as.integer(chains)
    if (M < 1L || S1 %% M != 0L) {
      stop("`chains` gives ", M, " chains but the ", S1,
           " posterior draws cannot be split into equal chains; pass chain ids.",
           call. = FALSE)
    }
    return(rep(S1 %/% M, M))
  }
  if (length(chains) != S1) {
    stop("`chains` must be NULL, the number of chains, or a vector of length S1.",
         call. = FALSE)
  }
  r <- rle(as.character(chains))
  if (anyDuplicated(r$values)) {
    stop("Chain ids must form contiguous runs (concatenate the chains).", call. = FALSE)
  }
  r$lengths
}

.resolve_G <- function(G, M, chain_lengths, S2, K = 1L) {
  G <- as.integer(G)
  if (G %% M != 0L) {
    G_new <- M * max(1L, as.integer(round(G / M)))
    if (G_new < 3L) G_new <- M * ceiling(3 / M)
    message("G = ", G, " is not a multiple of the number of chains (", M,
            "); using G = ", G_new, ".")
    G <- G_new
  }
  if (G < 3L) stop("At least 3 jackknife blocks are needed.", call. = FALSE)
  if (G %% K != 0L) stop("G must be a multiple of the number of folds.", call. = FALSE)
  if (min(chain_lengths) < G %/% M || S2 < G) {
    stop("Too few draws for G = ", G, " jackknife blocks.", call. = FALSE)
  }
  G
}

.no_overlap_result <- function(chat, root, settings) {
  na <- NA_real_
  out <- list(
    logml = chat, mcse = na, n_ov = 0, n_ovP = 0, tau = na, V_jack = na,
    interval = c(lower = -Inf, upper = Inf), flags = c(overlap = TRUE, variance = NA),
    decision = "do not report log Z-hat +/- z MCSE",
    actions = "the samples do not overlap: improve the proposal before adding draws",
    draws_factor = na, critical = c(n0 = settings$n0, V_jack = na),
    I = 0, NH0 = 0, n_P = 0, mcse_iid = na, V_post = na, V_prop = na, folds = NULL, jackknife = NULL,
    solver = root[c("iterations", "converged", "status")], settings = settings
  )
  class(out) <- "bridge_diagnostics"
  out
}

#' @export
print.bridge_diagnostics <- function(x, digits = 4, ...) {
  st <- x$settings
  f <- function(v, d = digits) if (is.null(v) || is.na(v)) "NA" else format(signif(v, d), scientific = FALSE)
  pf <- function(flag) if (is.na(flag)) "" else if (flag) "  FLAG" else "  pass"
  cat("Bridge sampling (",
      if (st$iid) "independent draws" else sprintf("MCMC, %d chain%s", st$M, if (st$M > 1) "s" else ""),
      "; S1 = ", st$S1, ", S2 = ", st$S2, if (isTRUE(st$K > 1L)) sprintf(", %d folds", st$K) else "", ")\n", sep = "")
  if (!is.finite(x$logml)) {
    cat("  log Z-hat = ", x$logml, ": the samples do not overlap (", x$solver$status, ").\n", sep = "")
    return(invisible(x))
  }
  cat(sprintf("  log Z-hat        %s   MCSE %s\n", f(x$logml, 8), f(x$mcse)))
  cat(sprintf("  Overlap check    %s = %s (>= %s)%s%s\n", if (st$iid) "n_ov-hat" else "n_ovP-hat",
              f(if (st$iid) x$n_ov else x$n_ovP), f(st$n0, 3),
              if (st$iid) "" else sprintf(", tau-hat = %s", f(x$tau, 3)), pf(x$flags[["overlap"]])))
  if (isTRUE(st$jackknife)) {
    cat(sprintf("  Variance check   V_jack = %s (<= %s, G = %d)%s\n", f(x$V_jack, 3), f(x$critical[["V_jack"]], 3),
                st$G, pf(x$flags[["variance"]])))
  }
  if (all(is.finite(x$interval)) || st$iid) {
    cat(sprintf("  %g%% finite-sample interval [%s, %s]\n", 100 * (1 - st$delta),
                f(x$interval[["lower"]], 8), f(x$interval[["upper"]], 8)))
  }
  cat("  Decision: ", x$decision, "\n", sep = "")
  for (a in x$actions) cat("    - ", a, "\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.bridge_diagnostics <- function(x, row.names = NULL, optional = FALSE, ...) {
  st <- x$settings
  data.frame(
    logml = x$logml, mcse = x$mcse, n_ov = x$n_ov, n_ovP = x$n_ovP, tau = x$tau, V_jack = x$V_jack,
    interval_lower = x$interval[["lower"]], interval_upper = x$interval[["upper"]],
    flag_overlap = x$flags[["overlap"]], flag_variance = x$flags[["variance"]],
    pass = !any(x$flags %in% TRUE), draws_factor = x$draws_factor,
    I = x$I, NH0 = x$NH0, n_P = x$n_P, mcse_iid = x$mcse_iid,
    S1 = st$S1, S2 = st$S2, M = st$M, G = st$G, K = st$K %||% 1L, iid = st$iid,
    tau_method = st$tau_method, row.names = row.names, stringsAsFactors = FALSE
  )
}

#' @export
summary.bridge_diagnostics <- function(object, ...) as.data.frame(object)

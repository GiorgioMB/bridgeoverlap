#' Bridge sampling with diagnostics
#'
#' Computes the optimal bridge sampling estimate of \eqn{\log Z} together with
#' the recommended diagnostics: the delta-method
#' MCSE (with a long-run variance for Markov chain output), the effective
#' overlap certificate \eqn{\hat n_{ov}} with \eqn{\hat I} and
#' \eqn{N\hat H_0}, the MCMC certificate \eqn{\hat n^P_{ov}}, the finite-sample
#' interval (independent draws), the calibration test for
#' sampler failure and the delete-a-group jackknife comparisons
#' \eqn{\mathcal V_{jack}}, \eqn{\mathcal I_{\max}} and \eqn{\mathcal G_{\max}}.
#' No additional evaluations of the posterior or the proposal are needed.
#'
#' @param x Either an object of class `bridge` returned by
#'   [bridgesampling::bridge_sampler()] with `method = "normal"` and
#'   `repetitions = 1`, or a numeric vector of log ratios
#'   \eqn{\log f - \log g} at the posterior draws, in sampling order with
#'   chains concatenated.
#' @param a_prop Numeric vector of log ratios at the proposal draws (only when
#'   `x` is numeric).
#' @param chains Chain structure of the posterior draws: `NULL` (a single
#'   sequence), the number of equal-length concatenated chains, or a vector of
#'   chain identifiers of length \eqn{S_1} with contiguous runs. For `bridge`
#'   objects created by [bridge_sampler_cmdstanr()] or [bridge_crossfit()] it is
#'   read from the object; see [bridge_chain_ids()] for other objects.
#' @param iid Logical: are the posterior draws independent (exact draws)? If
#'   `TRUE`, \eqn{\hat\tau = 1}, the finite-sample interval and the
#'   near-perfect accuracy bound are computed, and the overlap threshold is
#'   `thresholds$n_ov`; otherwise \eqn{\hat\tau} is estimated with `tau_method`
#'   and the threshold is `thresholds$n_ovP`.
#' @param G Number of jackknife blocks; a multiple of the number of chains
#'   (adjusted if necessary).
#' @param alpha Overall level of the four failure tests, each run at
#'   `alpha / 4`.
#' @param delta Error probability for the finite-sample interval and bounds.
#' @param tau_method Estimator of the integrated autocorrelation time of the
#'   posterior responsibilities; see [iat()]. Ignored when `iid = TRUE`.
#' @param thresholds List with the overlap thresholds `n_ov` (default 10, for
#'   independent draws), `n_ovP` (default 30, for MCMC) and `I_gate` (default
#'   1/2), the value of \eqn{\hat I} above which the threshold is relaxed.
#' @param eb_method Finite-sample interval: `"bennett"` (default; from Bennett's inequality) or `"mp"` (an empirical Bernstein interval).
#' @param rule How the certificate is read. `"plain"` (default): flag if the
#'   certificate is below `threshold`. `"continuous"`: flag if
#'   the certificate is below `threshold * min(1, (1 - I-hat) / (1 - I_gate))`,
#'   so that the threshold is constant under poor overlap and shrinks linearly
#'   in the near-perfect regime, where it amounts to requiring about
#'   `2 * threshold` ambiguous draws; the decision changes continuously with
#'   \eqn{\hat I}. `"gate"`: the gated rule, which reads the
#'   certificate only when \eqn{\hat I \le} `I_gate` and otherwise requires the
#'   lower confidence bound for \eqn{NH} to reach \eqn{16\log(2/\delta)}.
#' @param folds For cross-fitted output: a list with integer
#'   vectors `post` (length \eqn{S_1}) and `prop` (length \eqn{S_2}) giving
#'   the fold of every draw. Variances, autocorrelation times, jackknife blocks
#'   and the finite-sample interval are then computed within folds. Read from
#'   the object for [bridge_crossfit()] output.
#' @param nsim Number of simulations for the null quantiles of the jackknife
#'   influence statistics ([jackknife_null_quantiles()]).
#' @param eb Logical: compute the finite-sample interval?
#'   Defaults to `iid`; the interval is not valid for Markov chain output.
#' @param jackknife Logical: compute the jackknife diagnostics and the
#'   calibration test? Setting `FALSE` saves time in large simulation studies.
#' @param tol,maxit Tolerance and iteration limit of the root solver.
#' @param keep_responsibilities Logical: store the fitted responsibilities?
#' @param ... Passed to methods.
#'
#' @return An object of class `bridge_diagnostics`: a list whose scalar
#'   elements include `logml` (the estimate \eqn{\hat c}), `mcse`
#'   (\eqn{\hat\sigma}, or \eqn{\hat\sigma_P} for MCMC), `mcse_iid` (ignoring
#'   autocorrelation), `mcse_log1p` (the \eqn{\{\log(1+\hat\sigma^2)\}^{1/2}}
#'   convention of \pkg{bridgesampling}), `n_ov`, `I`, `NH0`, `NH`, `n_P`,
#'   `n_ovP`, `tau`, `stuck_chains` (chains whose responsibilities are
#'   constant, i.e. the sampler never moved), `T`, `t_T`, `p_T`, `V_jack`,
#'   `I_max`, `G_max`, the
#'   finite-sample interval `eb`, the lower confidence bound `NH_lower` for
#'   \eqn{NH} and the absolute `accuracy_bound` (independent draws), together
#'   with `flags`, `thresholds`, `jackknife` and `settings`. Use
#'   [as.data.frame()] for a one-row summary.
#'
#' @details
#' Interpretation:
#' * **Overlap.** When \eqn{\hat I \le 1/2} (poor overlap), \eqn{\hat n_{ov}}
#'   counts the ambiguous draws and controls the accuracy of the estimate and
#'   of its MCSE; the run is flagged if \eqn{\hat n_{ov} < 10} (independent
#'   draws) or \eqn{\hat n^P_{ov} < 30} (MCMC). When \eqn{\hat I > 1/2} the
#'   proposal is nearly exact and a small \eqn{\hat n_{ov}} is harmless; the run
#'   is then flagged only if the lower confidence bound `NH_lower` for
#'   \eqn{NH} is below \eqn{16\log(2/\delta)}, in which case
#'   the accuracy of the estimate cannot be certified. Otherwise, for independent
#'   draws, an absolute accuracy bound is reported.
#' * **Calibration test.** \eqn{\mathcal T = 1 - \hat H/\hat H_0}, the
#'   derivative of the fixed-point map at the solution, studentized with the
#'   jackknife and compared with \eqn{t_{G-1}}. It detects samplers whose error
#'   varies with the responsibilities on the overlap region (e.g. mode
#'   imbalance). No statistic of the bridge output can detect errors that are
#'   constant on the overlap region.
#' * **Jackknife.** \eqn{\mathcal V_{jack}} compares a batch-type variance with
#'   the MCSE (flags an underestimated autocorrelation), \eqn{\mathcal I_{\max}}
#'   and \eqn{\mathcal G_{\max}} flag a single dominating block.
#'
#' The estimate uses \eqn{s_1 = S_1/N}. It coincides with the `logml` of
#' [bridgesampling::bridge_sampler()] when that is called with
#' `use_neff = FALSE`; with the default `use_neff = TRUE`, \pkg{bridgesampling}
#' replaces \eqn{S_1} by a median effective sample size inside the bridge
#' function, which changes the estimator slightly.
#'
#' @examples
#' set.seed(1)
#' S <- 2000
#' x <- rnorm(S); y <- rnorm(S, 0.5, 1.4)            # posterior and proposal draws
#' a <- function(t) dnorm(t, log = TRUE) - dnorm(t, 0.5, 1.4, log = TRUE)  # log Z = 0
#' d <- bridge_diagnostics(a(x), a(y), iid = TRUE, nsim = 2e4)
#' d
#' as.data.frame(d)
#' @export
bridge_diagnostics <- function(x, ...) UseMethod("bridge_diagnostics")

#' @rdname bridge_diagnostics
#' @export
bridge_diagnostics.default <- function(x, a_prop, chains = NULL, iid = FALSE,
                                       G = 20L, alpha = 0.05, delta = 0.05,
                                       tau_method = c("within", "multichain",
                                                      "concat_geyer", "concat_ar"),
                                       thresholds = list(),
                                       rule = c("plain", "continuous", "gate"),
                                       eb_method = c("bennett", "mp"),
                                       folds = NULL, nsim = 1e5L,
                                       eb = iid, jackknife = TRUE,
                                       tol = 1e-12, maxit = 500L,
                                       keep_responsibilities = FALSE, ...) {
  tau_method <- match.arg(tau_method)
  rule <- match.arg(rule); eb_method <- match.arg(eb_method)
  .check_scalar(alpha, "alpha", 0, 1)
  .check_scalar(delta, "delta", 0, 1)
  thr <- utils::modifyList(list(n_ov = 10, n_ovP = 30, I_gate = 0.5), thresholds)
  a_post <- .clean_log_ratios(x, "posterior")
  a_prop <- .clean_log_ratios(a_prop, "proposal")
  S1 <- length(a_post)
  S2 <- length(a_prop)
  if (S1 < 2L || S2 < 2L) stop("Both samples need at least two draws.", call. = FALSE)
  N <- S1 + S2
  s1 <- S1 / N
  s2 <- S2 / N
  lr <- log(S1 / S2)
  # without a chain structure, every fold of cross-fitted output is its own sequence
  if (is.null(chains) && !is.null(folds) && length(folds$post) == S1) chains <- folds$post
  chain_lengths <- .chain_lengths(chains, S1)
  M <- length(chain_lengths)
  chain_id <- rep.int(seq_len(M), chain_lengths)
  fd <- .check_folds(folds, S1, S2, chain_id)
  K <- fd$K
  G <- if (jackknife) .resolve_G(G, M, chain_lengths, S2, K) else NA_integer_

  # Step 2: the estimate (root of the pooled estimating function).
  root <- .bridge_root(a_post, a_prop, lr, tol = tol, maxit = maxit)
  chat <- root$c
  settings <- list(S1 = S1, S2 = S2, N = N, s1 = s1, s2 = s2, M = M,
                   chain_lengths = chain_lengths, G = G, K = K, alpha = alpha,
                   delta = delta, iid = iid, eb = eb, jackknife = jackknife,
                   rule = rule, tau_method = if (iid) "iid" else tau_method)
  if (!is.finite(chat)) return(.no_overlap_result(chat, root, settings, thr))

  # Step 3: responsibilities, overlap and contraction factor (pooled over folds).
  pp <- .resp(a_post, chat, lr)
  pq <- .resp(a_prop, chat, lr)
  NH0 <- sum(pp)
  NH <- sum(.resp_var(a_post, chat, lr)) + sum(.resp_var(a_prop, chat, lr))
  I_hat <- (NH / N) / (s1 * s2)
  T_stat <- 1 - NH / NH0

  # Step 4: certificate and MCSE. Variances are computed within folds: with
  # cross-fitting each fold has its own responsibility function.
  fold_stats <- lapply(seq_len(K), function(k) {
    ip <- fd$post == k; iq <- fd$prop == k
    V_post <- stats::var(pp[ip]); V_prop <- stats::var(pq[iq])
    tau <- if (iid) 1 else iat(pp[ip], rle(chain_id[ip])$lengths, method = tau_method)
    c(S1 = sum(ip), S2 = sum(iq), V_post = V_post, V_prop = V_prop, tau = tau,
      a = sum(pp[ip]), b = sum(1 - pq[iq]))
  })
  fs <- do.call(rbind, fold_stats)
  V_post <- if (K == 1L) fs[1, "V_post"] else sum(fs[, "S1"] * fs[, "V_post"]) / S1
  V_prop <- if (K == 1L) fs[1, "V_prop"] else sum(fs[, "S2"] * fs[, "V_prop"]) / S2
  post_term <- sum(fs[, "S1"] * fs[, "V_post"])
  n_ov <- sum(fs[, "S2"] * fs[, "V_prop"]) + post_term
  n_P <- sum(fs[, "S2"] * fs[, "V_prop"]) + sum(fs[, "S1"] * fs[, "V_post"] * fs[, "tau"])
  tau <- if (iid) 1 else if (post_term > 0) sum(fs[, "S1"] * fs[, "V_post"] * fs[, "tau"]) / post_term else NA_real_
  # chains whose responsibilities never change: the sampler did not move
  stuck_chains <- if (iid) 0L else sum(vapply(
    split(pp, chain_id),
    function(v) length(v) > 1L && stats::var(v) == 0, logical(1)))
  mcse <- sqrt(n_P) / NH0
  mcse_iid <- sqrt(n_ov) / NH0
  n_ovP <- if (n_P > 0) n_ov^2 / n_P else NA_real_

  # Step 5: finite-sample results, reported in every regime. Lower confidence bound
  # for N H, per fold at level delta / K; for independent draws the
  # accuracy bound and the finite-sample interval.
  lg <- log(2 * K / delta)
  NH_lower_k <- pmax(sqrt(pmin(fs[, "a"], fs[, "b"])) - sqrt(2 / 3 * lg), 0)^2
  NH_lower <- sum(NH_lower_k)
  NH_required <- 16 * lg
  accuracy_bound <- if (iid && all(NH_lower_k >= NH_required)) max(2 * sqrt(lg / NH_lower_k)) else NA_real_
  eb_int <- c(lower = NA_real_, upper = NA_real_)
  if (eb) eb_int <- if (eb_method == "bennett") {
    .fs_interval(a_post, a_prop, lr, chat, delta, fold_post = if (K > 1L) fd$post, fold_prop = if (K > 1L) fd$prop)
  } else {
    .eb_interval(a_post, a_prop, lr, chat, delta, fold_post = if (K > 1L) fd$post, fold_prop = if (K > 1L) fd$prop)
  }

  # Steps 6-7: jackknife, with blocks within chains and within folds.
  c_b <- T_b <- NA_real_
  if (jackknife) {
    bpost <- .block_ids_post(chain_lengths, G)
    bprop <- if (K == 1L) .block_ids_prop(S2, G) else .block_ids_prop_folds(fd$prop, G)
    jk <- .jackknife(a_post, a_prop, lr, chat, bpost, bprop, G, tol, maxit)
    c_b <- jk$c_b
    T_b <- jk$T_b
  }
  if (jackknife && all(is.finite(c_b)) && all(is.finite(T_b))) {
    c_bar <- mean(c_b)
    V_jack <- (G - 1) / G * sum((c_b - c_bar)^2) / mcse^2
    I_max <- max(abs(c_b - chat)) / mcse
    Q <- G * chat - (G - 1) * c_b
    G_max <- if (stats::sd(Q) > 0) max(abs(Q - mean(Q))) / stats::sd(Q) else NA_real_
    se_T <- sqrt((G - 1) / G * sum((T_b - mean(T_b))^2))
    t_T <- if (se_T > 0) T_stat / se_T else NA_real_
  } else {
    V_jack <- I_max <- G_max <- t_T <- se_T <- NA_real_
  }
  p_T <- if (is.finite(t_T)) 2 * stats::pt(-abs(t_T), df = G - 1) else NA_real_

  # Step 8: the overlap check. Plain rule (default): a fixed threshold
  # on the certificate, so near-perfect runs are flagged and get the finite-sample
  # interval. The continuous rule (threshold shrinking linearly with 1 - I-hat above
  # the gate) and the gated rule are kept for comparison.
  crit <- list(t_T = NA_real_, V_jack = NA_real_, I_max = NA_real_, G_max = NA_real_)
  if (jackknife) {
    q_jack <- jackknife_null_quantiles(G, probs = 1 - alpha / 4, nsim = nsim)
    crit <- list(
      t_T = stats::qt(1 - alpha / 8, df = G - 1),
      V_jack = stats::qchisq(1 - alpha / 4, df = G - 1) / (G - 1),
      I_max = q_jack["I_max", 1],
      G_max = q_jack["G_max", 1]
    )
  }
  base <- if (iid) thr$n_ov else thr$n_ovP
  stat <- if (iid) n_ov else n_ovP
  weight <- max(0, min(1, (1 - I_hat) / (1 - thr$I_gate)))
  if (rule == "plain") {                 # a fixed threshold; near perfection no relaxation can be certified
    n_ov_threshold <- base
    overlap_flag <- !isTRUE(stat >= base)
  } else if (rule == "continuous") {
    n_ov_threshold <- base * weight
    overlap_flag <- !isTRUE(stat >= n_ov_threshold)
  } else {
    n_ov_threshold <- if (I_hat <= thr$I_gate) base else NA_real_
    overlap_flag <- if (I_hat <= thr$I_gate) !isTRUE(stat >= base) else NH_lower < NH_required
  }
  certificate <- if (iid) "n_ov" else "n_ovP"
  overlap_reason <- if (!overlap_flag) "" else if (rule != "gate" || I_hat <= thr$I_gate) {
    sprintf("%s-hat = %s < %s", certificate, .fmt(stat, 3), .fmt(n_ov_threshold, 3))
  } else {
    "near-perfect regime not certified: lower bound for N H below 16 log(2/delta)"
  }
  test <- function(stat, cv) if (jackknife && is.finite(stat)) stat > cv else NA
  flags <- c(
    overlap = overlap_flag,
    calibration = test(abs(t_T), crit$t_T),
    jackknife_variance = test(V_jack, crit$V_jack),
    influential_block_imax = test(I_max, crit$I_max),
    influential_block_grubbs = test(G_max, crit$G_max)
  )
  dec <- .decide(flags, stat = stat, threshold = n_ov_threshold, stuck = stuck_chains,
                 iid = iid, eb_int = eb_int, I_hat = I_hat)

  out <- list(
    logml = chat, mcse = mcse, mcse_iid = mcse_iid,
    mcse_log1p = sqrt(log1p(mcse^2)),
    mcse_mw = sqrt(max(1 / NH0 - 1 / (N * s1 * s2), 0)),
    n_ov = n_ov, I = I_hat, NH0 = NH0, NH = NH, n_P = n_P, n_ovP = n_ovP,
    tau = tau, stuck_chains = stuck_chains, V_post = V_post, V_prop = V_prop,
    T = T_stat, t_T = t_T, p_T = p_T, V_jack = V_jack, I_max = I_max,
    G_max = G_max, eb = eb_int, NH_lower = NH_lower, NH_required = NH_required,
    accuracy_bound = accuracy_bound, n_ov_threshold = n_ov_threshold,
    overlap_reason = overlap_reason,
    flags = flags, decision = dec$decision, actions = dec$actions,
    draws_factor = dec$draws_factor,
    thresholds = c(thr, crit),
    regime = if (I_hat <= thr$I_gate) "poor overlap" else "near-perfect proposal",
    folds = if (K > 1L) fs else NULL,
    jackknife = if (jackknife) list(c_b = c_b, T_b = T_b, se_T = se_T) else NULL,
    solver = root[c("iterations", "converged", "status")],
    settings = settings
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

# Decision rule: pass if no flag is raised; otherwise map each flag
# to an action. The overlap action quantifies "more draws in proportion to 1/H":
# n_ov grows linearly in N at a fixed proposal, so the factor is threshold / n_ov.
.decide <- function(flags, stat, threshold, stuck, iid, eb_int, I_hat = NA_real_) {
  raised <- names(flags)[flags %in% TRUE]
  if (isTRUE(stuck > 0)) raised <- c(raised, "stuck_chains")
  factor <- if (isTRUE(flags[["overlap"]]) && is.finite(stat) && stat > 0 && is.finite(threshold)) {
    max(1, threshold / stat)
  } else NA_real_
  act <- c(
    overlap = if (isTRUE(I_hat > 0.5)) {
      paste0("near-perfect proposal: no statistic of the bridge output can certify the normal interval here, ",
             "and more draws do not raise n_ov; ", if (iid) "report the finite-sample interval, which is short in this regime" else
             "thin the chain to a lag with a mixing bound and report the finite-sample interval of the thinned draws")
    } else if (is.finite(factor)) {
      sprintf("rerun from scratch with both samples at least %s times larger (the draws needed scale as 1/H; do not extend this run until it passes), or improve the proposal (warp-III, heavier tails, mixture)%s",
              .fmt(factor, 2), if (iid) "; meanwhile report the finite-sample interval" else "")
    } else "the samples (almost) do not overlap: improve the proposal before adding draws",
    calibration = "recheck the sampler (dispersed initial values, compare chains, mode balance) and refit the proposal by cross-fitting",
    jackknife_variance = "the MCSE is too small: run longer chains and estimate the long-run variance within chains",
    influential_block_imax = "one block dominates: inspect it for a transient or a rare-mode visit; run longer chains or recheck the sampler",
    influential_block_grubbs = "one block dominates: inspect it for a transient or a rare-mode visit; run longer chains or recheck the sampler",
    stuck_chains = "a chain never moved: recheck the sampler before bridge sampling"
  )
  list(decision = if (length(raised)) "do not report the normal interval" else "report log Z-hat +/- z MCSE",
       actions = unique(unname(act[raised])), draws_factor = factor)
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

.no_overlap_result <- function(chat, root, settings, thr) {
  na <- NA_real_
  out <- list(
    logml = chat, mcse = na, mcse_iid = na, mcse_log1p = na, mcse_mw = na,
    n_ov = 0, I = 0, NH0 = 0, NH = 0, n_P = 0, n_ovP = 0, tau = na, stuck_chains = NA_integer_,
    V_post = na, V_prop = na, T = na, t_T = na, p_T = na, V_jack = na,
    I_max = na, G_max = na,
    eb = c(lower = if (settings$iid) -Inf else na, upper = if (settings$iid) Inf else na),
    NH_lower = na, NH_required = na, accuracy_bound = na, n_ov_threshold = na,
    overlap_reason = "the samples do not overlap",
    flags = c(overlap = TRUE, calibration = NA, jackknife_variance = NA,
              influential_block_imax = NA, influential_block_grubbs = NA),
    decision = "do not report the normal interval",
    actions = "the samples do not overlap: improve the proposal before adding draws",
    draws_factor = NA_real_,
    thresholds = thr, regime = "no overlap", folds = NULL,
    jackknife = NULL, solver = root[c("iterations", "converged", "status")],
    settings = settings
  )
  class(out) <- "bridge_diagnostics"
  out
}

#' @export
print.bridge_diagnostics <- function(x, digits = 4, ...) {
  st <- x$settings
  f <- function(v, d = digits) if (is.null(v) || is.na(v)) "NA" else format(signif(v, d), scientific = FALSE)
  cat("Bridge sampling diagnostics (",
      if (st$iid) "independent draws" else sprintf("MCMC, %d chain%s", st$M, if (st$M > 1) "s" else ""),
      "; S1 = ", st$S1, ", S2 = ", st$S2,
      if (isTRUE(st$jackknife)) paste0(", G = ", st$G) else "", ")\n", sep = "")
  if (!is.finite(x$logml)) {
    cat("  log Z-hat = ", x$logml, ": the samples do not overlap (", x$solver$status, ").\n", sep = "")
    return(invisible(x))
  }
  cat(sprintf("  log Z-hat          %s   MCSE %s%s\n", f(x$logml, 8), f(x$mcse),
              if (st$iid) "" else sprintf("   (ignoring autocorrelation: %s)", f(x$mcse_iid))))
  if (isTRUE(st$eb)) {
    cat(sprintf("  %g%% finite-sample interval  [%s, %s]\n", 100 * (1 - st$delta),
                f(x$eb[["lower"]], 8), f(x$eb[["upper"]], 8)))
  }
  cat(sprintf("  Overlap            n_ov-hat %s   I-hat %s   N*H0-hat %s   (%s; threshold %s)\n",
              f(x$n_ov), f(x$I, 3), f(x$NH0), x$regime, f(x$n_ov_threshold, 3)))
  if (!st$iid) {
    cat(sprintf("  MCMC               tau-hat %s (%s)   n_ovP-hat %s\n",
                f(x$tau), st$tau_method, f(x$n_ovP)))
    if (isTRUE(x$stuck_chains > 0)) {
      cat(sprintf("  WARNING: %d of %d chains have constant responsibilities (the sampler did not move).\n",
                  x$stuck_chains, st$M))
    }
  }
  cat(sprintf("  Calibration test   T %s   t_T %s   p %s\n", f(x$T), f(x$t_T, 3), f(x$p_T, 3)))
  if (isTRUE(st$jackknife)) {
    cat(sprintf("  Jackknife          V_jack %s   I_max %s   G_max %s\n",
                f(x$V_jack, 3), f(x$I_max, 3), f(x$G_max, 3)))
  }
  if (isTRUE(st$K > 1L)) cat(sprintf("  Cross-fitted       %d folds; variances, blocks and bounds within folds\n", st$K))
  if (is.finite(x$NH_lower %||% NA)) {
    cat(sprintf("  Finite-sample      N*H >= %s (prob. >= %g)%s\n", f(x$NH_lower), 1 - st$delta,
                if (st$iid && is.finite(x$accuracy_bound))
                  sprintf(";  |log Z-hat - log Z| <= %s  (prob. >= %g)", f(x$accuracy_bound), 1 - 2 * st$delta)
                else if (st$iid) sprintf(";  accuracy bound needs N*H >= %s", f(x$NH_required)) else ""))
  }
  fl <- x$flags
  lab <- c(overlap = "overlap certificate", calibration = "calibration test (sampler failure)",
           jackknife_variance = "jackknife variance (MCSE too small)",
           influential_block_imax = "influential block (I_max)",
           influential_block_grubbs = "influential block (G_max)")
  if (isTRUE(fl[["overlap"]]) && nzchar(x$overlap_reason %||% "")) lab[["overlap"]] <- sprintf("overlap (%s)", x$overlap_reason)
  raised <- names(fl)[fl %in% TRUE]
  if (length(raised)) {
    cat(sprintf("  FLAGS (alpha = %g, failure tests at alpha/4): %s\n", st$alpha,
                paste(lab[raised], collapse = "; ")))
  } else {
    cat(sprintf("  No flags raised (alpha = %g, failure tests at alpha/4).\n", st$alpha))
  }
  cat("  Decision: ", x$decision %||% "", "\n", sep = "")
  for (a in x$actions %||% character()) cat("    - ", a, "\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.bridge_diagnostics <- function(x, row.names = NULL, optional = FALSE, ...) {
  st <- x$settings
  data.frame(
    logml = x$logml, mcse = x$mcse, mcse_iid = x$mcse_iid, mcse_log1p = x$mcse_log1p,
    n_ov = x$n_ov, I = x$I, NH0 = x$NH0, NH = x$NH, n_P = x$n_P, n_ovP = x$n_ovP,
    tau = x$tau, stuck_chains = x$stuck_chains, T = x$T, t_T = x$t_T, p_T = x$p_T, V_jack = x$V_jack,
    I_max = x$I_max, G_max = x$G_max,
    eb_lower = x$eb[["lower"]], eb_upper = x$eb[["upper"]],
    NH_lower = x$NH_lower, accuracy_bound = x$accuracy_bound,
    n_ov_threshold = x$n_ov_threshold %||% NA_real_, draws_factor = x$draws_factor %||% NA_real_,
    pass = !any(x$flags %in% TRUE) && !isTRUE(x$stuck_chains > 0),
    overlap_reason = x$overlap_reason %||% "",
    flag_overlap = x$flags[["overlap"]],
    flag_calibration = x$flags[["calibration"]],
    flag_jackknife_variance = x$flags[["jackknife_variance"]],
    flag_influential_imax = x$flags[["influential_block_imax"]],
    flag_influential_grubbs = x$flags[["influential_block_grubbs"]],
    flag_any_failure = if (all(is.na(x$flags[-1L]))) NA else any(x$flags[-1L] %in% TRUE),
    S1 = st$S1, S2 = st$S2, M = st$M, G = st$G, K = st$K %||% 1L, iid = st$iid,
    rule = st$rule %||% "gate", tau_method = st$tau_method,
    row.names = row.names, stringsAsFactors = FALSE
  )
}

#' @export
summary.bridge_diagnostics <- function(object, ...) as.data.frame(object)

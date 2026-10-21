# Helpers for cmdstanr ----------------------------------------------------------
#
# Bridge sampling needs log f at arbitrary points of the unconstrained space,
# i.e. the model's log_prob with the Jacobian adjustment. cmdstanr provides it
# through "model methods", which must be compiled together with the model.
# cmdstanr skips compilation (and therefore the model methods) when an
# up-to-date executable already exists, so cmdstan_model_bridge() checks for
# the methods and forces a recompilation if they are missing.
#
# Note: cmdstanr's log_prob() drops the constants of `~` statements (it calls
# log_prob_propto). For log Z to be the marginal likelihood, write the model
# block with `target += *_lpdf(...)` and include the normalizing constants of
# truncated priors (e.g. + log(2) for a half-normal or half-Cauchy prior). For
# calibration studies against a reference computed with the same program, the
# dropped constants are irrelevant because they shift log Z by a constant.

.require_cmdstanr <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required: install.packages('cmdstanr', ",
         "repos = c('https://stan-dev.r-universe.dev', getOption('repos')))",
         call. = FALSE)
  }
}

.has_model_methods <- function(mod) {
  env <- tryCatch(mod$.__enclos_env__$private$model_methods_env_,
                  error = function(e) NULL)
  !is.null(env) && length(env$hpp_code_) > 0 && is.function(env$model_ptr)
}

#' Compile a Stan program with the model methods needed for bridge sampling
#'
#' Wrapper around [cmdstanr::cmdstan_model()] with
#' `compile_model_methods = TRUE` that forces a recompilation when an existing
#' executable would otherwise be reused without the model methods (which makes
#' `$log_prob()` and `$unconstrain_draws()` unavailable). Programs in a
#' directory that is not writable are compiled from a copy in
#' `tools::R_user_dir("bridgeoverlap", "cache")`. If the environment variable
#' `BO_STAN_BUILD_DIR` is set, programs are always compiled from a copy in that
#' directory: concurrent processes (e.g. cluster array tasks) must each use their
#' own, because every new session rebuilds the executable.
#'
#' @param stan_file Path to the Stan program.
#' @param ... Further arguments passed to [cmdstanr::cmdstan_model()] and
#'   `$compile()`, e.g. `cpp_options` or `include_paths`.
#' @param quiet Passed to [cmdstanr::cmdstan_model()].
#' @return A `CmdStanModel` object.
#' @export
cmdstan_model_bridge <- function(stan_file, ..., quiet = TRUE) {
  .require_cmdstanr()
  # cmdstanr writes the executable next to the .stan file; programs in a
  # read-only library (e.g. those shipped with this package) are copied to the
  # user cache directory first.
  build <- Sys.getenv("BO_STAN_BUILD_DIR", "")
  if (nzchar(build) || file.access(dirname(normalizePath(stan_file)), 2L) != 0L) {
    cache <- if (nzchar(build)) build else tools::R_user_dir("bridgeoverlap", "cache")
    dir.create(cache, showWarnings = FALSE, recursive = TRUE)
    target <- file.path(cache, basename(stan_file))
    if (!file.exists(target) || !identical(readLines(target), readLines(stan_file))) {
      file.copy(stan_file, target, overwrite = TRUE)
    }
    stan_file <- target
  }
  mod <- cmdstanr::cmdstan_model(stan_file, compile_model_methods = TRUE,
                                 quiet = quiet, ...)
  if (!.has_model_methods(mod)) {
    mod$compile(compile_model_methods = TRUE, force_recompile = TRUE,
                quiet = quiet, ...)
  }
  mod
}

.cmdstanr_log_posterior <- function(pars, data) {
  out <- tryCatch(data$log_prob(unname(pars), jacobian = TRUE),
                  error = function(e) -Inf)
  if (is.na(out)) -Inf else out
}

#' Bridge sampling for a cmdstanr fit
#'
#' Runs [bridgesampling::bridge_sampler()] on the unconstrained draws of a
#' `CmdStanMCMC` fit, with `log f` evaluated by the model's `$log_prob()`
#' (including the Jacobian), and records the chain structure of the posterior
#' draws that enter the bridge so that [bridge_diagnostics()] can form blocks
#' and long-run variances within chains.
#'
#' @param fit A `CmdStanMCMC` object whose model was compiled with model
#'   methods, e.g. by [cmdstan_model_bridge()].
#' @param split How the draws are split between fitting the proposal and the
#'   bridge: `"per_chain"` (the first half of every chain fits the proposal, as
#'   in the `stanfit` method of \pkg{bridgesampling}; the default) or
#'   `"stacked"` (the first half of the stacked draws, as in the `CmdStanMCMC`
#'   method of \pkg{bridgesampling} >= 1.2-1).
#' @param crossfit Logical: use two-fold cross-fitting ([bridge_crossfit()])
#'   instead of a single split? Uses all draws for both roles.
#' @param use_neff Passed to [bridgesampling::bridge_sampler()]. `FALSE` (the
#'   default) gives the optimal bridge estimator (\eqn{s_1 = S_1/N}).
#' @param silent Passed to [bridgesampling::bridge_sampler()].
#' @param ... Further arguments for [bridgesampling::bridge_sampler()], e.g.
#'   `maxiter`.
#' @return A `bridge` object (see [bridgesampling::bridge_sampler()]) with the
#'   attribute `"bridgeoverlap_chains"`; pass it to [bridge_diagnostics()].
#' @examples
#' \dontrun{
#' file <- bridgeoverlap_stan_file("eight_schools_noncentered_ml")
#' mod <- cmdstan_model_bridge(file)
#' fit <- mod$sample(data = eight_schools_data(), refresh = 0, seed = 1)
#' bs <- bridge_sampler_cmdstanr(fit)
#' bridge_diagnostics(bs)        # exact log Z = -31.3113
#' }
#' @export
bridge_sampler_cmdstanr <- function(fit, split = c("per_chain", "stacked"),
                                    crossfit = FALSE, use_neff = FALSE,
                                    silent = TRUE, ...) {
  .require_cmdstanr()
  split <- match.arg(split)
  if (!inherits(fit, "CmdStanMCMC")) {
    stop("`fit` must be a CmdStanMCMC object.", call. = FALSE)
  }
  tryCatch(fit$init_model_methods(), error = function(e) {
    stop("The model methods of this fit are not available (", conditionMessage(e),
         "). Compile the model with cmdstan_model_bridge(), or with ",
         "cmdstanr::cmdstan_model(..., compile_model_methods = TRUE, force_recompile = TRUE).",
         call. = FALSE)
  })
  ud <- posterior::as_draws_array(fit$unconstrain_draws(format = "draws_array"))
  n_iter <- posterior::niterations(ud)
  n_chains <- posterior::nchains(ud)
  vars <- posterior::variables(ud)
  arr <- unclass(ud)
  chain_mats <- lapply(seq_len(n_chains), function(m) {
    mat <- arr[, m, , drop = FALSE]
    dim(mat) <- c(n_iter, length(vars))
    colnames(mat) <- vars
    mat
  })
  lb <- stats::setNames(rep(-Inf, length(vars)), vars)
  ub <- stats::setNames(rep(Inf, length(vars)), vars)
  common <- list(log_posterior = .cmdstanr_log_posterior, data = fit, lb = lb, ub = ub)
  if (crossfit) {
    samples <- coda::mcmc.list(lapply(chain_mats, coda::mcmc))
    out <- do.call(bridge_crossfit, c(list(samples), common,
                                      list(use_neff = use_neff, silent = silent), list(...)))
  } else if (split == "per_chain") {
    samples <- coda::mcmc.list(lapply(chain_mats, coda::mcmc))
    out <- do.call(bridgesampling::bridge_sampler,
                   c(list(samples), common,
                     list(use_neff = use_neff, silent = silent, repetitions = 1), list(...)))
    attr(out, "bridgeoverlap_chains") <- bridge_chain_ids(n_iter, n_chains, "per_chain")
  } else {
    samples <- do.call(rbind, chain_mats)
    out <- do.call(bridgesampling::bridge_sampler,
                   c(list(samples), common,
                     list(use_neff = use_neff, silent = silent, repetitions = 1), list(...)))
    attr(out, "bridgeoverlap_chains") <- bridge_chain_ids(n_iter, n_chains, "stacked")
  }
  attr(out, "bridgeoverlap_iid") <- FALSE
  out
}

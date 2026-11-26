make_gauss_bridge <- function(seed = 40, use_neff = FALSE, chains = 1) {
  set.seed(seed)
  d <- 3; n <- 2000; logZ <- -1.2
  lp <- function(pars, data) logZ + sum(dnorm(pars, 1, 2, log = TRUE))
  draw <- function() { m <- matrix(rnorm(n * d, 1, 2), n, d); colnames(m) <- paste0("b", 1:d); m }
  lb <- setNames(rep(-Inf, d), paste0("b", 1:d)); ub <- -lb
  samples <- if (chains == 1) draw() else coda::mcmc.list(lapply(1:chains, function(i) coda::mcmc(draw())))
  list(bs = bridgesampling::bridge_sampler(samples, log_posterior = lp, data = NULL, lb = lb, ub = ub,
                                           use_neff = use_neff, silent = TRUE),
       logZ = logZ, samples = samples, lp = lp, lb = lb, ub = ub, n = n)
}

test_that("the estimate equals bridgesampling's logml when use_neff = FALSE", {
  g <- make_gauss_bridge()
  d <- bridge_diagnostics(g$bs, iid = TRUE, nsim = 5000)
  expect_equal(d$logml, g$bs$logml, tolerance = 1e-8)
  expect_equal(d$bridgesampling$difference, 0, tolerance = 1e-8)
  expect_lt(abs(d$logml - g$logZ), 5 * d$mcse)
})

test_that("chain ids follow bridgesampling's splitting rules", {
  expect_equal(bridge_chain_ids(1000, 4, "per_chain"), rep(1:4, each = 500))
  expect_equal(bridge_chain_ids(1000, 4, "stacked"), rep(3:4, each = 1000))
  expect_equal(length(bridge_chain_ids(999, 3, "stacked")), 3 * 999 - round(3 * 999 / 2))
  g <- make_gauss_bridge(chains = 4)
  expect_length(g$bs$q11, 4 * 1000)
  d <- bridge_diagnostics(g$bs, chains = bridge_chain_ids(g$n, 4, "per_chain"), nsim = 5000)
  expect_equal(d$settings$M, 4L)
})

test_that("cross-fitting returns a bridge object usable by bridge_diagnostics", {
  g <- make_gauss_bridge(chains = 2)
  # default: two chains are folded by chains (each chain is bridged whole)
  cf <- bridge_crossfit(g$samples, log_posterior = g$lp, data = NULL, lb = g$lb, ub = g$ub)
  expect_s3_class(cf, "bridge")
  expect_identical(attr(cf, "bridgeoverlap_fold_by"), "chains")
  expect_length(cf$q11, 2 * g$n)
  d <- suppressMessages(bridge_diagnostics(cf, G = 20, nsim = 5000))
  expect_equal(d$settings$M, 2L)
  expect_equal(d$settings$K, 2L)
  expect_equal(d$logml, cf$logml, tolerance = 1e-10)
  expect_lt(abs(d$logml - g$logZ), 5 * d$mcse)
  # halves of every chain (exact only for independent draws)
  ch <- bridge_crossfit(g$samples, log_posterior = g$lp, data = NULL, lb = g$lb, ub = g$ub, fold = "halves")
  expect_identical(attr(ch, "bridgeoverlap_fold_by"), "halves")
  dh <- suppressMessages(bridge_diagnostics(ch, nsim = 5000))
  expect_equal(dh$settings$M, 4L)
  expect_lt(abs(dh$logml - g$logZ), 5 * dh$mcse)
})

test_that("folding by chains needs an even number of chains", {
  g <- make_gauss_bridge(chains = 3)
  expect_error(bridge_crossfit(g$samples, log_posterior = g$lp, data = NULL, lb = g$lb, ub = g$ub, fold = "chains"),
               "even number")
  cf <- bridge_crossfit(g$samples, log_posterior = g$lp, data = NULL, lb = g$lb, ub = g$ub)
  expect_identical(attr(cf, "bridgeoverlap_fold_by"), "halves")
})

test_that("unsupported objects give informative errors", {
  g <- make_gauss_bridge()
  w <- g$bs; w$method <- "warp3"
  expect_error(bridge_diagnostics(w), "warp3")
  bl <- structure(list(logml = c(1, 2)), class = "bridge_list")
  expect_error(bridge_diagnostics(bl), "repetitions = 1")
})

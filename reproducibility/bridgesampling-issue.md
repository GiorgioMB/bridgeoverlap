# Normalization of the posterior term in `error_measures()` and `mcse_logml`

**Package:** bridgesampling 1.1-2 (CRAN) and 1.2-1 (GitHub, `master`).
**Patch:** `bridgesampling-mcse.patch` (two lines).
**Status:** submitted as pull request 50, https://github.com/quentingronau/bridgesampling/pull/50
(opened 16 October 2026; open as of 21 October 2026). The pull request carries this patch; the
autocorrelation estimated across concatenated chains (`use_ess = TRUE`) is raised there as a remark.

## 1. `error_measures()` multiplies the spectral density by `var(f2)` twice

In `R/error_measures.R`,

```r
rho_f2 <- spectrum0.ar(f2)$spec
term1 <- 1 / N2 * var(f1) / mean(f1)^2
term2 <- rho_f2 / N1 * var(f2) / mean(f2)^2
```

`coda::spectrum0.ar(x)$spec` is the spectral density at frequency zero, i.e. the long-run
variance of `x` (it already equals `tau * var(x)`, with `tau` the integrated autocorrelation
time; `coda::effectiveSize(x)` is `length(x) * var(x) / spectrum0.ar(x)$spec`). Multiplying by
`var(f2)` again makes the posterior term too small by the factor `var(f2)`. With the default
`s1 = s2 = 1/2`, `f2` lies in `[0, 2]`, so the factor is at most 1, and it is small whenever the
two samples overlap poorly. The relative error reported by `error_measures()` is therefore too
small for MCMC and for independent draws alike. Correct:

```r
term2 <- rho_f2 / N1 / mean(f2)^2
```

**Reproducible example** (normalized target, so the exact answer is `log Z = 0`; four AR(1)
chains standing in for MCMC output):

```r
library(bridgesampling); library(coda)
set.seed(2026)
ar1 <- function(n, rho) { x <- numeric(n); x[1] <- rnorm(1)
  for (t in 2:n) x[t] <- rho * x[t - 1] + sqrt(1 - rho^2) * rnorm(1); x }
lp <- function(p, data) dnorm(p, log = TRUE)
res <- t(replicate(300, {
  ch <- mcmc.list(lapply(1:4, function(m) mcmc(matrix(ar1(1000, 0.9), dimnames = list(NULL, "x")))))
  b <- bridge_sampler(ch, log_posterior = lp, data = NULL, lb = c(x = -Inf), ub = c(x = Inf),
                      silent = TRUE, use_neff = FALSE)
  c(err = b$logml, cv = error_measures(b)$cv)
}))
sd(res[, "err"])                                     # 0.0054
median(res[, "cv"])                                  # 0.0011
mean(abs(res[, "err"]) <= 1.96 * res[, "cv"])        # 0.37 instead of 0.95
```

With the patch, the reported error is of the order of the actual standard deviation.

## 2. `mcse_logml` divides the posterior term by the number of proposal draws

In `R/bridge_sampler_internals.R`,

```r
var_r <- (mean_numi^2) / (mean_deni^2) * (var_numi / (mean_numi)^2 + var_deni / mean_deni^2)
var_r <- var_r / length(numi)
```

`numi` has one element per proposal draw (`n.2`) and `deni` one per posterior draw (`n.1`), so
the posterior term should be divided by `length(deni)`. The two coincide only when
`n.1 == n.2`, which is the default; any other allocation of proposal draws gives a wrong MCSE.
Correct:

```r
var_r <- (mean_numi^2) / (mean_deni^2) *
  (var_numi / (mean_numi)^2 / length(numi) + var_deni / mean_deni^2 / length(deni))
```

A related remark: with `use_ess = TRUE` the effective sample size is computed by
`coda::effectiveSize()` on the concatenation of the chains, whose junctions are spurious
transitions; estimating the autocorrelation within chains (or with the multi-chain estimator of
Vehtari et al. 2021) is safer.

Found with the R package `bridgeoverlap` (https://github.com/GiorgioMB/bridgeoverlap).

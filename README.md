# bridgeoverlap

Overlap certificates and failure diagnostics for bridge sampling in R, built on
[`bridgesampling`](https://github.com/quentingronau/bridgesampling) and
[`cmdstanr`](https://mc-stan.org/cmdstanr/).

The optimal bridge sampling estimator of Meng and Wong (1996) is the intercept of a logistic
regression that separates posterior draws from proposal draws. Its estimating function is a
sum of bounded terms, so the accuracy of the estimate and of its Monte Carlo standard error
(MCSE) is governed by the **overlap** of the two samples. From the output that 
`bridgesampling::bridge_sampler()` already produces, and at no extra cost, `bridgeoverlap` 
reports the estimate, its MCSE, an effective-overlap certificate, a finite-sample confidence 
interval, a calibration test for sampler failure and delete-a-group jackknife diagnostics, 
and maps every flag to an action.



## Installation

```r
# install.packages("remotes")
remotes::install_github("GiorgioMB/bridgeoverlap")

# optional, for Stan models
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

Works with `bridgesampling` 1.1-2 (CRAN) and 1.2-1 (GitHub).

## Quick start with cmdstanr

```r
library(bridgeoverlap)
mod <- cmdstan_model_bridge(bridgeoverlap_stan_file("eight_schools_noncentered_ml"))
fit <- mod$sample(data = eight_schools_data(), seed = 11, refresh = 0)
bs  <- bridge_sampler_cmdstanr(fit)     # bridgesampling on the unconstrained draws
bridge_diagnostics(bs)
eight_schools_log_evidence()            # -31.31135
```

```
Bridge sampling diagnostics (MCMC, 4 chains; S1 = 2000, S2 = 2000, G = 20)
  log Z-hat          -31.32138   MCSE 0.01072   (ignoring autocorrelation: 0.01014)
  Overlap            n_ov-hat 85.51   I-hat 0.907   N*H0-hat 912   (near-perfect proposal; threshold 5.59)
  MCMC               tau-hat 1.321 (within)   n_ovP-hat 76.44
  Calibration test   T 0.005749   t_T 1.87   p 0.0765
  Jackknife          V_jack 1.15   I_max 0.49   G_max 1.94
  Finite-sample      N*H >= 819.8 (prob. >= 0.95)
  No flags raised (alpha = 0.05, failure tests at alpha/4).
  Decision: report log Z-hat +/- z MCSE
```

`cmdstan_model_bridge()` compiles the model with the model methods that
bridge sampling needs (`$log_prob()`, `$unconstrain_draws()`). cmdstanr
silently skips them when an up-to-date executable exists, so the helper detects
that case and recompiles.

Two-fold cross-fitting uses every draw both to fit the proposal and in the
bridge. In the run above it lowered the MCSE from 0.0109 to 0.0079.

```r
bridge_diagnostics(bridge_sampler_cmdstanr(fit, crossfit = TRUE), G = 16)
```

When a flag is raised, `bridge_diagnostics()` prints the action (for the overlap flag, the
factor by which a fresh, independent run should be larger; do not extend the flagged run), and `bridge_report()` returns the printout as a
table. For a Bayes factor, run each model once and combine:

```r
bridge_bayes_factor(bridge_diagnostics(bs1), bridge_diagnostics(bs2))
```

## Quick start with any sampler

Any `bridge` object with `method = "normal"` and `repetitions = 1` works. Tell
the diagnostics how the posterior draws in the bridge are arranged in chains.

```r
library(bridgesampling)
bs <- bridge_sampler(samples_mcmc_list, log_posterior = lp, data = dat,
                     lb = lb, ub = ub, use_neff = FALSE)
bridge_diagnostics(bs, chains = bridge_chain_ids(n_iter, n_chains, "per_chain"))

# exact (independent) posterior draws: adds the finite-sample interval
bridge_diagnostics(bs, iid = TRUE)

# or directly from log ratios a = log f - log g at posterior and proposal draws
bridge_diagnostics(a_post, a_prop, chains = 4)
```

## What is reported and how to read it

| Output | Meaning | Rule |
|---|---|---|
| `logml`, `mcse` | Estimate and delta-method MCSE. For MCMC the posterior term uses the long-run variance, estimated within chains. |  |
| `n_ov`, `I`, `NH0` | Effective-overlap certificate, overlap ratio Î and the number of ambiguous draws. | Flag n̂_ov < 10 (independent draws; plain rule, the default). Near perfection (Î → 1) every run is flagged and the finite-sample interval is the answer. |
| `n_ovP`, `tau` | MCMC certificate n̂_ov² / n̂_P and the autocorrelation time of the responsibilities. | Flag n̂_ov^P < 30 (MCMC). |
| `eb` | Finite-sample interval from Bennett's inequality, valid at every sample size. | Independent draws only. |
| `NH_lower`, `accuracy_bound` | Lower confidence bound for N·H and, for independent draws, the absolute accuracy bound. | Reported in every regime; the bound needs `NH_lower` ≥ 16 log(2/δ). |
| `decision`, `actions`, `draws_factor` | Decision rule. | Report the normal interval only if no flag is raised. |
| `T`, `t_T`, `p_T` | Calibration test for sampler failure. | Flag if \|t_T\| > t_{G-1, 1-α/8}. |
| `V_jack` | Jackknife variance divided by the squared MCSE: an MCSE that is too small. | Flag if above χ²_{G-1, 1-α/4} / (G-1). |
| `I_max`, `G_max` | A single dominating block, such as a transient or a stuck chain. | Simulated pivotal quantiles at 1-α/4. |
| `stuck_chains` | Chains whose responsibilities never change, i.e. the sampler did not move. | Printed as a warning. |

The four failure tests run at level α/4 each, so a correct run raises any of
them with probability at most α. Use `as.data.frame()` for one row per run.

## Things to know

- **Normalizing constants.** cmdstanr's `$log_prob()` drops the constants of
  `~` statements. For log Z to be the marginal likelihood:
  - write the model block with `target += *_lpdf(...)`;
  - normalize truncated priors, for example `+ log2()` for a half-Cauchy (see
    `inst/stan/`).

  For calibration studies against a reference computed with the same program,
  the dropped constants don't matter.
- **`use_neff`.** By default `bridgesampling` puts a median effective sample
  size into the bridge function. The optimal estimator uses s1 = S1/N, which
  is `use_neff = FALSE`. This is the default in `bridge_sampler_cmdstanr()`, and
  then `bridge_diagnostics()` reproduces `bridgesampling`'s `logml` to about
  1e-15.
- **Near-perfect proposals.** When Î > 1/2 a small n̂_ov is harmless, so the threshold shrinks linearly with 1 − Î and amounts to requiring about
  twice the threshold in ambiguous draws. 
- **Cross-fitting.** `bridge_crossfit()` / `bridge_sampler_cmdstanr(fit, crossfit = TRUE)`
  record the folds; `bridge_diagnostics()` then computes variances, jackknife blocks and
  finite-sample bounds within folds. With an even number of chains the
  folds are groups of chains, which keeps every fold's estimating equation exactly unbiased
  under MCMC; with one chain or an odd number, halves of every chain (exact only for
  independent draws). Use `G` a multiple of the number of sequences.
- **Blind spot.** No statistic of the bridge output can detect a sampler error
  that is constant on the overlap region. The example is a mode
  that the posterior draws and the proposal both miss: see
  `analysis/12_stan_mixture.R`.
- **Supported methods.** Only `method = "normal"` is supported. Warp-III
  objects don't store the warp matrix needed to rebuild the log ratios.

## Reproducing the experiments

`analysis/` contains every experiment: the independent-draw and Markov chain studies, the
fitted-proposal and sampler-failure studies, the eight schools studies, the search over the
problem space, the computer-assisted constants, and the Stan benchmark on `posteriordb`
posteriors with exactly known evidence. Every table is written as a LaTeX `tabular`
(`tab_*.tex`), every figure as a PDF (`fig_*.pdf`), and every reported number as a LaTeX macro
(`results.tex`), all in `generated/` (or the directory named by `BO_OUT`):

```sh
Rscript analysis/run_all.R --stan    # everything (hours on one core; see hpc/ for a cluster)
Rscript analysis/run_all.R           # everything that needs no Stan
BO_QUICK=true BO_SCALE=0.01 Rscript analysis/run_all.R --stan   # smoke test, minutes
```

Replication counts are fixed in `analysis/config.R` and every replication has its own seed, so
a rerun reproduces every number. `generated/MANIFEST.txt` lists the outputs and whether the run
is final. Set `BO_TEX` to a colon-separated list of LaTeX files that use the macros, and the
last step also checks that each macro they use is defined. See
[`analysis/README.md`](analysis/README.md) for the map from outputs to scripts and
[`hpc/README.md`](hpc/README.md) for the cluster workflow.

**`bridgesampling` MCSEs.** `reproducibility/bridgesampling-mcse.patch` corrects two
normalization errors in the MCSEs of `bridgesampling`; it is submitted upstream as
[pull request 50](https://github.com/quentingronau/bridgesampling/pull/50).

## Tests

```r
testthat::test_local()                                # fast, no Stan needed
Sys.setenv(BRIDGEOVERLAP_TEST_CMDSTAN = "true")       # also compile and run Stan
testthat::test_local()
```

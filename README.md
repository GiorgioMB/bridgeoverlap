# bridgeoverlap

Overlap checks and failure diagnostics for the optimal bridge sampling estimator of
Meng and Wong (1996), as described in

> G. Micaletto (2026). *Overlap, not moments: sample size and standard errors for the
> optimal bridge sampling estimator.*

The optimal bridge estimator of log Z is the intercept of a logistic regression that
separates posterior from proposal draws. Its error is the sign of a sum of bounded terms,
so its accuracy is governed by the overlap of the two samples rather than by moments of
the importance ratios. From the fitted responsibilities, at no extra cost, the package
computes:

* the delta-method Monte Carlo standard error (MCSE), with a within-chain long-run
  variance for Markov chain output;
* the estimated effective overlap `n_ov` (and `n_ovP` = `n_ov / tau` for MCMC output),
  and the **overlap check** `n_ov >= 10` (`n_ovP >= 30` for MCMC);
* a finite-sample confidence interval, valid at every sample size for independent draws
  (Bennett's inequality applied to the estimating function);
* a calibration test for sampler failure and delete-a-group jackknife diagnostics for
  misestimated autocorrelation and dominating blocks;
* a decision rule that maps every flag to an action.

## Installation

```r
# install.packages("remotes")
remotes::install_github("GiorgioMB/bridgeoverlap")
```

## Quick start

```r
library(bridgesampling)
library(bridgeoverlap)

bs <- bridge_sampler(fit, ...)     # any bridgesampling object
d  <- bridge_diagnostics(bs)       # estimate, MCSE, n_ov, flags and the decision
d
bridge_report(d)                   # the same as a table

# from log ratios log{f/g} of posterior and proposal draws
d <- bridge_diagnostics(a_post, a_prop, chains = 4, iid = FALSE)
bridge_fs_interval(a_post, a_prop)  # finite-sample interval (independent draws)

# Bayes factor from two independent runs
bridge_bayes_factor(d1, d2)
```

`bridge_sampler_cmdstanr()` runs bridge sampling on a `cmdstanr` fit, and
`bridge_crossfit()` runs two-fold cross-fitting with folds made of whole chains.

## The decision rule

1. Check the sampler first: R-hat, effective sample sizes, divergences.
2. **Overlap flag** (`n_ov < 10`, or `n_ovP < 30`): report the finite-sample interval
   instead of `log Z-hat ± z MCSE`. Under poor overlap, improve the proposal or rerun from
   scratch with both sample sizes multiplied by the printed factor; near a perfect
   proposal, the finite-sample interval is the answer.
3. **Failure flags** (calibration test, jackknife variance ratio, influential blocks,
   stuck chains): recheck the sampler, estimate the long-run variance within chains,
   inspect the dominating block, fit the proposal by cross-fitting.
4. **No flag**: report `log Z-hat ± z MCSE` with `n_ov`, `I-hat` and the finite-sample
   interval.

No statistic of the bridge output can detect sampler errors that are constant on the
overlap region, such as label switching in mixtures; posteriors with known symmetries
need a proposal invariant under them.

## Reproducing the paper

Every number, table and figure of the paper and its supplement is produced by the
scripts in [`analysis/`](analysis/README.md): `Rscript analysis/run_all.R` runs every
study that needs no Stan, and `Rscript analysis/run_all.R --stan` adds the Stan
benchmark. Cluster scripts are in [`hpc/`](hpc/README.md).

`reproducibility/bridgesampling-mcse.patch` corrects two normalization errors in the
MCSE of `bridgesampling` (`error_measures()` and `mcse_logml`).

## License

MIT. Cite with `citation("bridgeoverlap")`.

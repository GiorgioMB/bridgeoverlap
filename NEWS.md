# bridgeoverlap 1.1.0

* `bridge_diagnostics()` reports the estimate, the MCSE, two checks and one interval:
  the overlap check (`n_ov`, or `n_ovP` for MCMC output) and the variance check (`V_jack`),
  with blocks aligned with the chains. The calibration test, the block-influence statistics
  `I_max` and `G_max`, the stuck-chain check, the lower bound for `N H` and the accuracy
  bound are no longer computed, and `jackknife_null_quantiles()` is removed.
* New argument `gamma`: with a lower bound on the spectral gap, the finite-sample interval
  is also reported for Markov chain output (`bridge_fs_interval_mcmc()`).
* `bridge_fs_interval_mcmc()` estimates the variance of both samples (empirical Bernstein
  for the proposal draws, a chain Bernstein bound for the posterior draws), which narrows
  the interval by a factor of two to three.
* Analysis: Pareto k-hat comparison (S24), multiple bridge sampling (S26), the
  finite-sample interval for Markov chain output (S25) and the collapse figure (S23).

# bridgeoverlap 1.0.0

First release, accompanying Micaletto (2026), *Overlap, not moments: sample size and
standard errors for the optimal bridge sampling estimator*.

* `bridge_diagnostics()`: estimate, MCSE, estimated effective overlap, overlap check,
  finite-sample interval and accuracy bound, calibration test, jackknife diagnostics
  and the decision rule (Algorithm 1 of the paper).
* `bridge_crossfit()`: two-fold cross-fitted bridge sampling (Algorithm S1).
* `bridge_fs_interval()` (independent draws) and `bridge_fs_interval_mcmc()` (Markov chain
  output with a known lower bound on the spectral gap), `bridge_report()`, `bridge_bayes_factor()`,
  `bridge_sampler_cmdstanr()` and helpers.

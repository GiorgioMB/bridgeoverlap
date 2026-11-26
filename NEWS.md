# bridgeoverlap 1.0.0

First release, accompanying Micaletto (2026), *Overlap, not moments: sample size and
standard errors for the optimal bridge sampling estimator*.

* `bridge_diagnostics()`: estimate, MCSE, estimated effective overlap, overlap check,
  finite-sample interval and accuracy bound, calibration test, jackknife diagnostics
  and the decision rule (Algorithm 1 of the paper).
* `bridge_crossfit()`: two-fold cross-fitted bridge sampling (Algorithm S1).
* `bridge_fs_interval()`, `bridge_report()`, `bridge_bayes_factor()`,
  `bridge_sampler_cmdstanr()` and helpers.

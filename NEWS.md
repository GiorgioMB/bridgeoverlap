# bridgeoverlap 0.2.0

First public release.

* `bridge_diagnostics()`: the optimal bridge sampling estimate (safeguarded Newton on the log
  scale), its delta-method MCSE with a within-chain long-run variance, the effective-overlap
  certificate, a lower confidence bound for `N H`, the finite-sample interval and accuracy
  bound for independent draws, the calibration test, the delete-a-group jackknife statistics,
  and a decision that maps every flag to an action. The overlap flag uses a fixed threshold
  (`rule = "plain"`); `rule = "continuous"` and `rule = "gate"` are alternatives.
* `bridge_fs_interval()`: the finite-sample interval from Bennett's inequality with a
  high-probability variance bound; `bridge_eb_interval()`: an empirical Bernstein interval.
* `bridge_report()`: the diagnostics of one run as a table; `bridge_bayes_factor()`: the
  diagnostics of a Bayes factor from two independent runs.
* `bridge_sampler_cmdstanr()` and `bridge_crossfit()`: bridge sampling of `cmdstanr` fits, with
  two-fold cross-fitting folded by chains.

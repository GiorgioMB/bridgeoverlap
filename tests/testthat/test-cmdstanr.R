# These tests compile a Stan model and are therefore opt-in:
#   Sys.setenv(BRIDGEOVERLAP_TEST_CMDSTAN = "true")
test_that("bridge sampling a cmdstanr fit recovers the exact eight schools evidence", {
  skip_if_not(identical(Sys.getenv("BRIDGEOVERLAP_TEST_CMDSTAN"), "true"),
              "set BRIDGEOVERLAP_TEST_CMDSTAN=true to run the cmdstanr tests")
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)))
  mod <- cmdstan_model_bridge(bridgeoverlap_stan_file("eight_schools_noncentered_ml"))
  fit <- mod$sample(data = eight_schools_data(), chains = 4, iter_warmup = 500,
                    iter_sampling = 1000, refresh = 0, seed = 1, show_messages = FALSE)
  bs <- bridge_sampler_cmdstanr(fit)
  d <- bridge_diagnostics(bs, nsim = 5000)
  expect_equal(d$settings$M, 4L)
  expect_lt(abs(d$logml - eight_schools_log_evidence()), 5 * d$mcse)
  cf <- bridge_diagnostics(bridge_sampler_cmdstanr(fit, crossfit = TRUE), G = 16, nsim = 5000)
  expect_equal(cf$settings$S1, 4000L)
  expect_lt(abs(cf$logml - eight_schools_log_evidence()), 5 * cf$mcse)
})

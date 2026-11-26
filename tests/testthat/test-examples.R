test_that("the exact eight schools evidence is reproduced", {
  expect_equal(eight_schools_log_evidence(), -31.311347, tolerance = 1e-6)
})


test_that("the shipped Stan programs exist", {
  for (nm in c("eight_schools_centered_ml", "eight_schools_noncentered_ml")) {
    expect_true(file.exists(bridgeoverlap_stan_file(nm)))
  }
  expect_error(bridgeoverlap_stan_file("nope"), "Unknown")
})

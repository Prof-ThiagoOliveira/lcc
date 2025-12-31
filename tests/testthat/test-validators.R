test_that("length and numeric validators enforce constraints", {
  expect_true(lcc:::validate_equal_length(1:3, 4:6))
  expect_error(lcc:::validate_equal_length(1:3, 1:2), class = "lcc_error_input")

  expect_true(lcc:::validate_numeric_no_na(c(0.1, 0.2)))
  expect_error(lcc:::validate_numeric_no_na("a"), class = "lcc_error_input")
  expect_error(lcc:::validate_numeric_no_na(c(1, NA)), class = "lcc_error_input")
  expect_error(lcc:::validate_numeric_no_na(c(1, Inf)), class = "lcc_error_input")
})

test_that("non-degenerate variance validator warns once variance collapses", {
  old <- getOption("lcc.show.warnings")
  on.exit(options(lcc.show.warnings = old), add = TRUE)
  options(lcc.show.warnings = TRUE)

  expect_true(lcc:::validate_non_degenerate_var(c(0, 1, 2)))
  expect_warning(res <- lcc:::validate_non_degenerate_var(rep(1, 3)), class = "lcc_warning")
  expect_false(res)
})

test_that("fisher helpers clamp extreme correlations", {
  expect_equal(lcc:::safe_clamp_r(c(-2, -0.5, 0.9, 2)), c(-1 + 1e-8, -0.5, 0.9, 1 - 1e-8))
  transformed <- lcc:::safe_fisher(c(-0.99, 0, 0.99))
  expect_type(transformed, "double")
  inverted <- lcc:::safe_fisher_inv(transformed)
  expect_equal(round(inverted, 6), round(c(-0.99, 0, 0.99), 6))
})

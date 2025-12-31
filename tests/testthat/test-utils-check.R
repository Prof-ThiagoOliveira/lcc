test_that("check helpers accept valid inputs", {
  expect_invisible(lcc:::check_flag(TRUE))
  expect_invisible(lcc:::check_scalar_integer(5L, lower = 1L, upper = 10L))
  expect_invisible(lcc:::check_scalar_numeric(0.5, lower = 0, upper = 1))
  expect_invisible(lcc:::check_character_scalar("method"))
  expect_invisible(lcc:::check_choice("alpha", c("alpha", "beta")))

  df <- data.frame(a = factor(c("x", "y")), b = c(1, 2), stringsAsFactors = TRUE)
  expect_invisible(lcc:::check_has_columns(df, c("a", "b")))
  expect_invisible(lcc:::check_is_factor_col(df, "a"))
  expect_invisible(lcc:::check_is_numeric_col(df, "b"))
  expect_invisible(lcc:::check_polynomial_degrees(2L, 1L))
  expect_invisible(lcc:::check_REML_flag(FALSE))
  expect_invisible(lcc:::check_num_core(1L))
  expect_invisible(lcc:::check_gs("x", df, "a"))
})

test_that("check helpers signal informative errors", {
  expect_error(lcc:::check_flag("TRUE"), class = "lcc_error_input")
  expect_error(lcc:::check_scalar_integer(2.5), class = "lcc_error_input")
  expect_error(lcc:::check_scalar_integer(0L, lower = 1L), class = "lcc_error_input")
  expect_error(lcc:::check_scalar_integer(11L, upper = 10L), class = "lcc_error_input")
  expect_error(lcc:::check_scalar_numeric("a"), class = "lcc_error_input")
  expect_error(lcc:::check_scalar_numeric(-1, lower = 0), class = "lcc_error_input")
  expect_error(lcc:::check_scalar_numeric(2, upper = 1), class = "lcc_error_input")
  expect_error(lcc:::check_character_scalar(TRUE), class = "lcc_error_input")
  expect_error(lcc:::check_choice("gamma", c("alpha", "beta")), class = "rlang_error")

  df <- data.frame(a = factor("x"), b = "non-numeric", stringsAsFactors = TRUE)
  expect_error(lcc:::check_has_columns(df, c("a", "missing")), class = "lcc_error_input")
  expect_error(lcc:::check_is_factor_col(data.frame(a = 1), "a", label = "A"), class = "lcc_error_input")
  expect_error(lcc:::check_is_numeric_col(data.frame(a = "text"), "a", label = "A"), class = "lcc_error_input")
  expect_error(lcc:::check_polynomial_degrees(1L, 2L), class = "lcc_error_input")
  expect_error(lcc:::check_REML_flag(NA), class = "lcc_error_input")
  expect_error(lcc:::check_num_core(0L), class = "lcc_error_input")

  df_methods <- data.frame(m = factor(c("A", "B")))
  expect_error(lcc:::check_gs("C", df_methods, "m"), class = "lcc_error_input")
})

test_that("warn_general honours lcc.show.warnings option", {
  old <- getOption("lcc.show.warnings")
  on.exit(options(lcc.show.warnings = old), add = TRUE)

  options(lcc.show.warnings = FALSE)
  expect_null(lcc:::warn_general("message suppressed"))

  options(lcc.show.warnings = TRUE)
  expect_warning(lcc:::warn_general("message visible"), class = "lcc_warning")
})

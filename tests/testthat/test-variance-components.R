test_that("extract_random_effects_cov returns covariance matrix and rank", {
  mat <- matrix(c(2, 0.5, 0.5, 1), nrow = 2, ncol = 2)
  res <- with_mocked_bindings(
    lcc:::extract_random_effects_cov("model"),
    getVarCov = function(model, type) {
      expect_identical(type, "random.effects")
      mat
    },
    .package = "nlme"
  )
  expect_equal(res$G, mat)
  expect_identical(res$n_re, 2L)
})

test_that("extract_random_effects_cov warns on heterogeneous groups", {
  mat1 <- matrix(c(2, 0.6, 0.6, 1.5), nrow = 2)
  mat2 <- mat1 + 5e-6
  old <- getOption("lcc.show.warnings")
  on.exit(options(lcc.show.warnings = old), add = TRUE)
  options(lcc.show.warnings = TRUE)

  res <- with_mocked_bindings(
    expect_warning(
      lcc:::extract_random_effects_cov("model"),
      class = "lcc_warning"
    ),
    getVarCov = function(model, type) list(mat1, mat2),
    .package = "nlme"
  )
  expect_equal(res$G, mat1)
  expect_identical(res$n_re, 2L)
})

test_that("extract_random_effects_cov catches invalid matrices", {
  mat_rect <- matrix(seq_len(6), nrow = 2, ncol = 3)
  expect_error(
    with_mocked_bindings(
      lcc:::extract_random_effects_cov("model"),
      getVarCov = function(model, type) mat_rect,
      .package = "nlme"
    ),
    class = "lcc_error_internal"
  )

  bad_list <- list(diag(2), matrix(1, nrow = 1, ncol = 1))
  expect_error(
    with_mocked_bindings(
      lcc:::extract_random_effects_cov("model"),
      getVarCov = function(model, type) bad_list,
      .package = "nlme"
    ),
    class = "lcc_error_internal"
  )

  expect_error(
    with_mocked_bindings(
      lcc:::extract_random_effects_cov("model"),
      getVarCov = function(model, type) NULL,
      .package = "nlme"
    ),
    class = "lcc_error_internal"
  )
})

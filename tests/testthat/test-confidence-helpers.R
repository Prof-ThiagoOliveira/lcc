test_that("build_ci_metric returns empty structure when bootstrap missing", {
  ci <- lcc:::build_ci_metric(
    boot_list = NULL,
    alpha = 0.1,
    transform = lcc:::ZFisher,
    inv_transform = lcc:::ZFisher_inv,
    percentile = FALSE,
    bounds = c(-1, 1)
  )
  expect_identical(dim(ci), c(2L, 0L))
  expect_identical(attr(ci, "ci_level"), 0.9)
  expect_identical(attr(ci, "ci_method"), "normal")

  empty_list <- list()
  ci2 <- lcc:::build_ci_metric(
    boot_list = empty_list,
    alpha = 0.2,
    transform = lcc:::ZFisher,
    inv_transform = lcc:::ZFisher_inv,
    percentile = TRUE,
    bounds = c(-1, 1)
  )
  expect_identical(dim(ci2), c(2L, 0L))
  expect_identical(attr(ci2, "ci_level"), 0.8)
  expect_identical(attr(ci2, "ci_method"), "percentile")
})

test_that("build_ci_metric falls back to NA matrix on computation error", {
  boot_list <- list(rnorm(5))
  res <- with_mocked_bindings(
    lcc:::build_ci_metric(
      boot_list = boot_list,
      alpha = 0.05,
      transform = lcc:::ZFisher,
      inv_transform = lcc:::ZFisher_inv,
      percentile = FALSE,
      bounds = c(-1, 1)
    ),
    .build_ci_from_boot = function(...) stop("boom"),
    .package = "lcc"
  )
  expect_identical(dim(res), c(2L, 5L))
  expect_true(all(is.na(res)))
})

test_that("build_ci_metric handles percentile and normal intervals", {
  boot_list_normal <- replicate(5, seq(-0.8, 0.8, length.out = 4), simplify = FALSE)
  ci_normal <- lcc:::build_ci_metric(
    boot_list = boot_list_normal,
    alpha = 0.1,
    transform = lcc:::ZFisher,
    inv_transform = lcc:::ZFisher_inv,
    percentile = FALSE,
    bounds = c(-1, 1)
  )
  expect_identical(dim(ci_normal), c(2L, length(boot_list_normal[[1]])))
  expect_true(all(ci_normal[1, ] <= ci_normal[2, ]))

  boot_list <- list(
    c(-0.8, -0.4, 0, 0.4),
    c(-0.6, -0.2, 0.2, 0.6),
    c(-0.7, -0.3, 0.3, 0.7)
  )
  ci_percentile <- lcc:::build_ci_metric(
    boot_list = boot_list,
    alpha = 0.2,
    transform = lcc:::ZFisher,
    inv_transform = lcc:::ZFisher_inv,
    percentile = TRUE,
    bounds = c(-1, 1)
  )
  expect_identical(dim(ci_percentile), c(2L, length(boot_list[[1]])))
  expect_true(all(ci_percentile[1, ] <= ci_percentile[2, ]))
})

test_that(".build_ci_from_boot clamps, warns and transforms", {
  mat <- matrix(c(-2, 0, 2, 0.5, 0.5, 0.5, Inf, NaN), nrow = 4)
  old <- getOption("lcc.show.warnings")
  on.exit(options(lcc.show.warnings = old), add = TRUE)
  options(lcc.show.warnings = TRUE)

  expect_warning(
    ci <- lcc:::.build_ci_from_boot(
      boot_list = mat,
      alpha = 0.1,
      transform = lcc:::ZFisher,
      inv_transform = lcc:::ZFisher_inv,
      percentile = FALSE,
      bounds = c(-1, 1),
      warn_label = "metric"
    ),
    class = "lcc_warning"
  )
  expect_identical(dim(ci), c(2L, 4L))
  expect_true(all(ci >= -1 & ci <= 1, na.rm = TRUE))
  expect_true(anyNA(ci))
})

test_that("CI alignment helpers reshape matrices and lists", {
  base_mat <- matrix(runif(4), nrow = 2)
  expect_error(lcc:::.align_ci_to_grid(1:5, len_time = 3), class = "lcc_error_internal")

  mat_transposed <- matrix(1:4, nrow = 4, ncol = 2)
  aligned <- lcc:::.align_ci_to_grid(mat_transposed, len_time = 3)
  expect_identical(dim(aligned), c(2L, 3L))

  long_mat <- matrix(runif(8), nrow = 2)
  trimmed <- lcc:::.align_ci_to_grid(long_mat, len_time = 3)
  expect_identical(dim(trimmed), c(2L, 3L))

  short_mat <- matrix(runif(2), nrow = 2)
  padded <- lcc:::.align_ci_to_grid(short_mat, len_time = 4)
  expect_identical(dim(padded), c(2L, 4L))

  ci_list <- list(long_mat, short_mat)
  aligned_list <- lcc:::.align_ci_list_to_grid(ci_list, len_time = 3)
  expect_identical(length(aligned_list), 2L)
  expect_identical(dim(aligned_list[[1]]), c(2L, 3L))
  expect_identical(dim(aligned_list[[2]]), c(2L, 3L))

  expect_error(lcc:::.align_ci_list_to_grid(1, len_time = 2), class = "lcc_error_internal")
})

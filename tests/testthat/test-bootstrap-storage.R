library(testthat)
library(lcc)

context("bootstrap storage and reuse")

data(hue)
set.seed(123)

test_that("bootstrap outputs use matrix storage for single comparison", {
  set.seed(123)
  fit <- lcc(
    data = hue, subject = "Fruit", resp = "H_mean",
    method = "Method", time = "Time",
    qf = 1, qr = 0, ci = FALSE, components = TRUE
  )
  
  tk_grid <- sort(unique(fit$model$data$time))
  
  bs <- lcc:::bootstrapSamples(
    nboot         = 5,
    model         = fit$model,
    q_f           = 1,
    q_r           = 0,
    interaction   = TRUE,
    covar         = NULL,
    var.class     = NULL,
    pdmat         = pdSymm,
    weights.form  = NULL,
    show.warnings = FALSE,
    tk            = tk_grid,
    diffbeta      = list(0),
    ldb           = 1,
    components    = TRUE,
    lme.control   = NULL,
    method.init   = "REML",
    numCore       = 1
  )
  
  expect_true(is.matrix(bs$LCC_Boot))
  expect_equal(dim(bs$LCC_Boot), c(length(tk_grid), 5))
  expect_true(is.matrix(bs$LPC_Boot))
  expect_equal(dim(bs$LPC_Boot), c(length(tk_grid), 5))
  expect_true(is.matrix(bs$Cb_Boot))
  expect_equal(dim(bs$Cb_Boot), c(length(tk_grid), 5))
})

test_that("cached basis can be reused for repeated calls", {
  set.seed(456)
  fit <- lcc(
    data = hue, subject = "Fruit", resp = "H_mean",
    method = "Method", time = "Time",
    qf = 1, qr = 0, ci = FALSE, components = FALSE
  )
  
  tk_grid <- sort(unique(fit$model$data$time))
  
  # First run (builds basis)
  bs1 <- lcc:::bootstrapSamples(
    nboot         = 2,
    model         = fit$model,
    q_f           = 1,
    q_r           = 0,
    interaction   = TRUE,
    covar         = NULL,
    var.class     = NULL,
    pdmat         = pdSymm,
    weights.form  = NULL,
    show.warnings = FALSE,
    tk            = tk_grid,
    diffbeta      = list(0),
    ldb           = 1,
    components    = FALSE,
    lme.control   = NULL,
    method.init   = "REML",
    numCore       = 1
  )
  
  # Second run should hit the cached basis path without error
  bs2 <- lcc:::bootstrapSamples(
    nboot         = 2,
    model         = fit$model,
    q_f           = 1,
    q_r           = 0,
    interaction   = TRUE,
    covar         = NULL,
    var.class     = NULL,
    pdmat         = pdSymm,
    weights.form  = NULL,
    show.warnings = FALSE,
    tk            = tk_grid,
    diffbeta      = list(0),
    ldb           = 1,
    components    = FALSE,
    lme.control   = NULL,
    method.init   = "REML",
    numCore       = 1
  )
  
  expect_true(is.matrix(bs1$LCC_Boot))
  expect_true(is.matrix(bs2$LCC_Boot))
})

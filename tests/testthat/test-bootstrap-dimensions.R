library(testthat)
library(lcc)

test_that("bootstrap outputs have expected dimensions", {
  data(hue)
  fit <- lcc(
    data = hue, subject = "Fruit", resp = "H_mean",
    method = "Method", time = "Time",
    qf = 1, qr = 0, ci = FALSE, components = FALSE
  )
  tk_grid <- sort(unique(fit$model$data$time))
  
  bs <- lcc:::bootstrapSamples(
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
  
  expect_true(is.matrix(bs$LCC_Boot))
  expect_equal(dim(bs$LCC_Boot), c(length(tk_grid), 2))
})

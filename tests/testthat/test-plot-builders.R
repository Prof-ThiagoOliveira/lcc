test_that("plotBuilder_lpc produces ggplot output with confidence intervals", {
  arg <- plotControl(plot = FALSE)
  model <- list(data = data.frame(method = factor(c("A", "B"), levels = c("A", "B"))))
  tk_plot <- c(0, 1)
  tk_plot2 <- c(0, 1)
  lpc_vals <- c(0.4, 0.5)
  env_lpc <- matrix(c(0.3, 0.35, 0.45, 0.55), nrow = 2, byrow = TRUE)
  pearson_vals <- list(data.frame(V1 = c(0.42, 0.48)))

  plt <- lcc:::plotBuilder_lpc(
    LPC         = lpc_vals,
    ENV.LPC     = env_lpc,
    tk.plot     = tk_plot,
    Pearson_vals = pearson_vals,
    tk.plot2    = tk_plot2,
    ldb         = 1,
    model       = model,
    ci          = TRUE,
    arg         = arg
  )

  expect_s3_class(plt, "ggplot")
  expect_equal(plt$labels$y, arg$ylab)
})

test_that("plotBuilder_la handles multiple comparisons without confidence intervals", {
  arg <- plotControl(plot = FALSE)
  model <- list(data = data.frame(method = factor(c("A", "B", "C"), levels = c("A", "B", "C"))))
  tk_plot  <- c(0, 1)
  tk_plot2 <- c(0, 1)
  CCC_vals <- list(
    data.frame(V1 = c(0.8, 0.75)),
    data.frame(V1 = c(0.7, 0.65))
  )
  Pearson_vals <- list(
    data.frame(V1 = c(0.85, 0.8)),
    data.frame(V1 = c(0.78, 0.74))
  )
  Cb <- matrix(c(0.6, 0.7, 0.55, 0.65), nrow = 2, byrow = FALSE)

  plt <- lcc:::plotBuilder_la(
    CCC_vals     = CCC_vals,
    Pearson_vals = Pearson_vals,
    Cb           = Cb,
    ENV.Cb       = NULL,
    tk.plot      = tk_plot,
    tk.plot2     = tk_plot2,
    ldb          = 2,
    model        = model,
    ci           = FALSE,
    arg          = arg
  )

  expect_s3_class(plt, "ggplot")
})

test_that("plotBuilder_lcc returns ggplot objects for basic inputs", {
  arg <- plotControl(plot = FALSE)
  model <- list(data = data.frame(method = factor(c("A", "B"), levels = c("A", "B"))))
  tk_plot  <- c(0, 1)
  tk_plot2 <- c(0, 1)
  rho <- c(0.5, 0.6)
  CCC_vals <- list(data.frame(V1 = c(0.48, 0.58)))

  plt <- lcc:::plotBuilder_lcc(
    rho      = rho,
    ENV.LCC  = NULL,
    tk.plot  = tk_plot,
    CCC      = CCC_vals,
    tk.plot2 = tk_plot2,
    ldb      = 1,
    model    = model,
    ci       = FALSE,
    arg      = arg
  )

  expect_s3_class(plt, "ggplot")
  expect_equal(plt$labels$title, "B vs. A")
})

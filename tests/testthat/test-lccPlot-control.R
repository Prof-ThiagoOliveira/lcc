test_that("plotControl axis overrides are honoured", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  data(hue)

  fm <- lcc(
    data    = hue,
    subject = "Fruit",
    resp    = "H_mean",
    method  = "Method",
    time    = "Time",
    qf      = 2, qr = 2,
    components = TRUE
  )

  ctrl <- plotControl(
    scale_y_continuous = list(
      limits = c(0.7, 1),
      breaks = seq(0.7, 1, by = 0.1)
    ),
    expand_y = c(0, 0)
  )

  p <- suppressWarnings(lccPlot(fm, control = ctrl))

  b <- ggplot2::ggplot_build(p)
  params <- b$layout$panel_params[[1]]

  # Extract range/breaks in a way that is robust across ggplot2 versions
  yrange <- params$y.range %||% params$y_range %||% params$y$range
  ybreaks <- params$y.breaks %||% params$y$breaks %||% params$y$get_breaks()

  expect_equal(yrange, c(0.7, 1))
  expect_equal(ybreaks, seq(0.7, 1, by = 0.1))
})

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

  first_non_null <- function(...) {
    candidates <- list(...)
    for (x in candidates) {
      if (!is.null(x)) {
        return(x)
      }
    }
    NULL
  }

  safe_axis_component <- function(axis, attr_name, accessor = NULL) {
    if (is.null(axis)) {
      return(NULL)
    }
    val <- axis[[attr_name]]
    if (!is.null(val)) {
      return(val)
    }
    if (!is.null(accessor) && is.function(axis[[accessor]])) {
      return(axis[[accessor]]())
    }
    NULL
  }

  # Extract range/breaks in a way that is robust across ggplot2 versions
  yrange <- first_non_null(
    params$y.range,
    params$y_range,
    safe_axis_component(params$y, "range", "range")
  )
  ybreaks <- first_non_null(
    params$y.breaks,
    safe_axis_component(params$y, "breaks", "get_breaks")
  )

  expect_equal(yrange, c(0.7, 1))
  expect_equal(ybreaks, seq(0.7, 1, by = 0.1))
})

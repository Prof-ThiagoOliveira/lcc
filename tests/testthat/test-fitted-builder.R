test_that("fittedBuilder builds data.frame for single comparison", {
  obj <- list(
    Summary.lcc = list(
      comp = "B vs A",
      fitted = data.frame(
        Time = c(0, 1),
        LCC  = c(0.4, 0.5),
        LPC  = c(0.55, 0.6),
        LA   = c(0.7, 0.72)
      )
    )
  )

  res_lcc <- lcc:::fittedBuilder(obj, type = "lcc")
  expect_s3_class(res_lcc, "data.frame")
  expect_equal(res_lcc$Methods, rep("B vs A", 2))
  expect_equal(res_lcc$`fitted.LCC`, c(0.4, 0.5))
})

test_that("fittedBuilder unwraps component lists for single comparison", {
  fitted_lists <- list(
    data.frame(Time = c(0, 1), LCC = c(0.4, 0.5), LPC = c(0.55, 0.6), LA = c(0.7, 0.72))
  )
  obj <- list(
    Summary.lcc = list(
      comp = "B vs A",
      fitted = list(
        fitted_lists[[1]],
        fitted_lists[[1]],
        fitted_lists[[1]]
      )
    )
  )

  res_lpc <- lcc:::fittedBuilder(obj, type = "lpc")
  expect_equal(res_lpc$`fitted.LPC`, c(0.55, 0.6))

  res_la <- lcc:::fittedBuilder(obj, type = "la")
  expect_equal(res_la$`fitted.LA`, c(0.7, 0.72))
})

test_that("fittedBuilder stacks multi-comparison outputs", {
  comp_labels <- list("B vs A", "C vs A")
  lcc_tables <- list(
    data.frame(Time = c(0, 1), LCC = c(0.4, 0.45)),
    data.frame(Time = c(0, 1), LCC = c(0.5, 0.55))
  )
  lpc_tables <- list(
    data.frame(Time = c(0, 1), LPC = c(0.6, 0.65)),
    data.frame(Time = c(0, 1), LPC = c(0.7, 0.75))
  )
  la_tables <- list(
    data.frame(Time = c(0, 1), LA = c(0.8, 0.82)),
    data.frame(Time = c(0, 1), LA = c(0.85, 0.87))
  )
  obj <- list(
    Summary.lcc = list(
      comp = comp_labels,
      fitted = list(
        LCC = lcc_tables,
        LPC = lpc_tables,
        LA  = la_tables
      )
    )
  )

  res_multi <- lcc:::fittedBuilder(obj, type = "lpc")
  expect_equal(res_multi$Methods, rep(unlist(comp_labels), each = 2))
  expect_equal(res_multi$Time, rep(c(0, 1), times = 2))
  expect_equal(res_multi$`fitted.LPC`, unlist(lapply(lpc_tables, `[[`, "LPC")))
})

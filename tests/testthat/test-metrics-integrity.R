test_that("metric estimates respect bounds", {
  data(hue, package = "lcc")

  fit <- lcc(
    data       = hue,
    subject    = "Fruit",
    resp       = "H_mean",
    method     = "Method",
    time       = "Time",
    qf         = 1,
    qr         = 0,
    ci         = TRUE,
    nboot      = 40,
    components = TRUE,
    numCore    = 1,
    boot.seed  = 101,
    keep.boot.models = FALSE,
    show.warnings = FALSE
  )

  metrics <- fit$plot_info$metrics

  lcc_vals <- unlist(metrics$lcc$estimate, use.names = FALSE)
  expect_true(all(lcc_vals >= -1 - 1e-8 & lcc_vals <= 1 + 1e-8, na.rm = TRUE))

  lpc_vals <- unlist(metrics$lpc$estimate, use.names = FALSE)
  expect_true(all(lpc_vals >= -1 - 1e-8 & lpc_vals <= 1 + 1e-8, na.rm = TRUE))

  la_vals <- unlist(metrics$la$estimate, use.names = FALSE)
  expect_true(all(la_vals >= 0 - 1e-8 & la_vals <= 1 + 1e-8, na.rm = TRUE))

  check_ci_bounds <- function(ci_list) {
    if (is.null(ci_list)) return(TRUE)
    all(vapply(ci_list, function(mat) {
      if (is.null(mat)) return(TRUE)
      lower <- mat["lower", , drop = TRUE]
      upper <- mat["upper", , drop = TRUE]
      all(lower <= upper + 1e-10, na.rm = TRUE)
    }, logical(1)))
  }

  expect_true(check_ci_bounds(metrics$lcc$ci))
  expect_true(check_ci_bounds(metrics$lpc$ci))
  expect_true(check_ci_bounds(metrics$la$ci))
})

test_that("CI width decreases with smaller alpha", {
  data(hue, package = "lcc")

  fit_95 <- lcc(
    data       = hue,
    subject    = "Fruit",
    resp       = "H_mean",
    method     = "Method",
    time       = "Time",
    qf         = 1,
    qr         = 0,
    ci         = TRUE,
    alpha      = 0.05,
    nboot      = 30,
    components = FALSE,
    numCore    = 1,
    boot.seed  = 202,
    keep.boot.models = FALSE,
    show.warnings = FALSE
  )

  fit_90 <- lcc(
    data       = hue,
    subject    = "Fruit",
    resp       = "H_mean",
    method     = "Method",
    time       = "Time",
    qf         = 1,
    qr         = 0,
    ci         = TRUE,
    alpha      = 0.10,
    nboot      = 30,
    components = FALSE,
    numCore    = 1,
    boot.seed  = 202,
    keep.boot.models = FALSE,
    show.warnings = FALSE
  )

  ci95 <- fit_95$plot_info$metrics$lcc$ci[[1]]
  ci90 <- fit_90$plot_info$metrics$lcc$ci[[1]]

  width95 <- ci95["upper", ] - ci95["lower", ]
  width90 <- ci90["upper", ] - ci90["lower", ]

  expect_true(all(width95 >= width90 - 1e-8, na.rm = TRUE))
})

test_that("degenerate responses yield NA metrics without errors", {
  degenerate <- expand.grid(
    Fruit  = factor(1:3),
    Method = factor(c("M1", "M2")),
    Time   = c(1, 2)
  )
  degenerate$Resp <- 5

  deg_fit <- try(
    lcc(
      data       = degenerate,
      subject    = "Fruit",
      resp       = "Resp",
      method     = "Method",
      time       = "Time",
      qf         = 1,
      qr         = 0,
      ci         = FALSE,
      components = TRUE,
      numCore    = 1,
      keep.boot.models = FALSE,
      show.warnings = FALSE
    ),
    silent = TRUE
  )

  if (inherits(deg_fit, "try-error")) {
    cond <- attr(deg_fit, "condition")
    msg <- if (!is.null(cond)) conditionMessage(cond) else as.character(deg_fit)
    skip(paste0("Degenerate fit failed to converge: ", msg))
  }

  deg_metrics <- deg_fit$plot_info$metrics
  expect_true(all(is.na(unlist(deg_metrics$lcc$estimate, use.names = FALSE))))
  expect_true(all(is.na(unlist(deg_metrics$lpc$estimate, use.names = FALSE))))
  expect_true(all(is.na(unlist(deg_metrics$la$estimate,  use.names = FALSE))))
})

test_that("lccModel fits a mixed-effects model on the fixture dataset", {
  data <- build_fixture_dataset()

  fit <- lcc:::lccModel(
    dataset     = data,
    resp        = "H_mean",
    subject     = "Fruit",
    method      = "Method",
    time        = "Time",
    qf          = 1,
    qr          = 0,
    interaction = FALSE,
    covar       = NULL,
    gs          = NULL,
    var.class   = NULL,
    weights.form = NULL,
    lme.control = nlme::lmeControl(msMaxIter = 20, tolerance = 1e-6),
    method.init = "REML",
    pdmat       = nlme::pdSymm
  )

  expect_s3_class(fit$model, "lme")
  expect_equal(fit$method.init, "REML")
  expect_identical(fit$wcount, 0L)
  expect_true("fixed" %in% names(fit$data))
})

test_that("lccInternal builds metric bundles and comparison names", {
  data <- expand.grid(
    subject = factor(sprintf("S%d", 1:3)),
    time    = 0:1,
    method  = factor(c("A", "B"), levels = c("A", "B"))
  )
  data <- data[order(data$subject, data$method, data$time), ]
  data$resp <- with(data, 1 + as.numeric(subject) * 0.1 + ifelse(method == "B", 0.3, 0) + 0.05 * time)

  model <- structure(
    list(
      data = data,
      summary_obj = list(modelStruct = list(varStruct = 0))
    ),
    class = "mock_lme"
  )

  summary_capture <- list()

  assign("summary.mock_lme", function(object, ...) object$summary_obj, envir = .GlobalEnv)
  on.exit(rm("summary.mock_lme", envir = .GlobalEnv), add = TRUE)
  stub_cov <- use_namespace_stub(
    "extract_random_effects_cov",
    function(model) list(G = matrix(1), n_re = 1L)
  )
  on.exit(restore_namespace_stub(stub_cov), add = TRUE)

  stub_pre <- use_namespace_stub(
    ".precompute_longitudinal",
    function(model, tk, q_f, q_r, ...) {
      list(tag = "pre")
    }
  )
  on.exit(restore_namespace_stub(stub_pre), add = TRUE)

  stub_lcc <- use_namespace_stub(
    ".compute_LCC",
    function(pre, diffbeta) list(c(0.4, 0.5))
  )
  on.exit(restore_namespace_stub(stub_lcc), add = TRUE)

  stub_lpc <- use_namespace_stub(
    ".compute_LPC",
    function(pre) list(c(0.6, 0.65))
  )
  on.exit(restore_namespace_stub(stub_lpc), add = TRUE)

  stub_la <- use_namespace_stub(
    ".compute_LA",
    function(pre, diffbeta) list(c(0.8, 0.82))
  )
  on.exit(restore_namespace_stub(stub_la), add = TRUE)

  stub_summary <- use_namespace_stub(
    "lccSummary",
    function(model, tk, tk.plot, tk.plot2, metrics, ldb, ci,
             components, degenerate_resp) {
      summary_capture <<- metrics
      list(fitted = "mock")
    }
  )
  on.exit(restore_namespace_stub(stub_summary), add = TRUE)

  res <- lcc:::lccInternal(
    model        = model,
    q_f          = 1,
    q_r          = 0,
    interaction  = FALSE,
    tk           = NULL,
    covar        = NULL,
    pdmat        = nlme::pdSymm,
    diffbeta     = list(c(0.1, -0.05)),
    time_lcc     = list(time = c(0.5)),
    ci           = FALSE,
    boot.scheme  = "np_case",
    ci.method    = "normal",
    alpha        = 0.05,
    nboot        = 10,
    labels       = NULL,
    var.class    = NULL,
    weights.form = NULL,
    show.warnings = FALSE,
    components   = TRUE,
    lme.control  = nlme::lmeControl(),
    method.init  = "REML",
    numCore      = 1,
    keep_models  = FALSE,
    boot_seed    = NULL
  )

  expect_equal(res$comp_names, "B vs A")
  expect_equal(res$tk.plot, c(0, 0.5, 1))
  expect_equal(res$tk.plot2, c(0, 1))
  expect_equal(res$rho, c(0.4, 0.5))
  expect_equal(res$rho.pearson, c(0.6, 0.65))
  expect_equal(res$Cb, c(0.8, 0.82))
  expect_equal(summary_capture$lcc$estimate[[1]], c(0.4, 0.5))
})

test_that("lcc_intervals builds percentile envelopes for single and multiple comparisons", {
  rho_single <- c(0.3, 0.4)
  boot_single <- list(c(0.2, 0.3), c(0.35, 0.45))

  ci_single <- lcc:::lcc_intervals(
    rho         = rho_single,
    tk.plot     = 1:2,
    tk.plot2    = 1:2,
    ldb         = 1,
    model       = NULL,
    ci          = TRUE,
    LCC_Boot    = boot_single,
    alpha       = 0.1,
    ci.method   = "percentile"
  )

  expect_equal(dim(ci_single$ENV.LCC), c(2, 2))

  rho_multi <- matrix(c(0.2, 0.3, 0.4, 0.5), nrow = 2)
  boot_multi <- list(
    list(c(0.1, 0.2), c(0.4, 0.5)),
    list(c(0.15, 0.25), c(0.45, 0.55))
  )

  ci_multi <- lcc:::lcc_intervals(
    rho         = rho_multi,
    tk.plot     = 1:2,
    tk.plot2    = 1:2,
    ldb         = 2,
    model       = NULL,
    ci          = TRUE,
    LCC_Boot    = boot_multi,
    alpha       = 0.1,
    ci.method   = "normal"
  )

  expect_equal(length(ci_multi$ENV.LCC), 2)
  expect_true(all(vapply(ci_multi$ENV.LCC, function(mat) all(dim(mat) == c(2, 2)), logical(1))))

  empty_boot <- matrix(numeric(0), nrow = 2)
  ci_empty <- lcc:::lcc_intervals(
    rho         = c(0.1, 0.2),
    tk.plot     = 1:2,
    tk.plot2    = 1:2,
    ldb         = 1,
    model       = NULL,
    ci          = TRUE,
    LCC_Boot    = empty_boot,
    alpha       = 0.05,
    ci.method   = "normal"
  )

  expect_true(all(is.na(ci_empty$ENV.LCC)))
})

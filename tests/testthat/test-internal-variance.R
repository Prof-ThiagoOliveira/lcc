test_that("getDelta defaults to constant scaling when variance structure is absent", {
  summary_stub <- list(modelStruct = list(varStruct = NULL))

  res <- lcc:::getDelta(model = NULL, summary_obj = summary_stub)

  expect_equal(res$delta, 0)
  expect_equal(res$deltal, 0)
  expect_equal(res$g(5), 1)
})

test_that("getDelta handles varIdent variance structures", {
  lambda <- log(1.3)
  var_struct <- structure(lambda, class = "varIdent")
  summary_stub <- list(modelStruct = list(varStruct = var_struct))

  res <- lcc:::getDelta(model = NULL, summary_obj = summary_stub)

  expect_equal(res$delta, 1)
  expect_equal(res$deltal, rep(exp(lambda)^2, length(lambda)))
  expect_equal(res$g(2), 2)
})

test_that("getDelta recognises varExp formulations", {
  var_exp_time <- structure(0.2, class = "varExp")
  attr(var_exp_time, "formula") <- ~time
  summary_time <- list(modelStruct = list(varStruct = var_exp_time))

  res_time <- lcc:::getDelta(model = NULL, summary_obj = summary_time)
  expect_s3_class(res_time$delta, "varExp")
  expect_equal(as.numeric(res_time$delta), 0.2)
  expect_equal(as.numeric(res_time$deltal), 0.2)
  expect_equal(res_time$g(0.2, tk = c(0, 1)), exp(c(0, 0.4)))

  var_exp_method <- structure(c(0.25, 0.1, -0.05), class = "varExp")
  attr(var_exp_method, "formula") <- ~time | method
  summary_method <- list(modelStruct = list(varStruct = var_exp_method))

  res_method <- lcc:::getDelta(model = NULL, summary_obj = summary_method)
  expect_equal(res_method$delta, 0.25)
  expect_equal(res_method$deltal, c(0.1, -0.05))
  expect_equal(res_method$g(0.25, tk = 0:1), exp(0.5 * 0:1))
})

test_that("getDelta errors on unsupported variance structures", {
  bogus <- structure(1, class = "varPower")
  summary_stub <- list(modelStruct = list(varStruct = bogus))

  expect_error(
    lcc:::getDelta(model = NULL, summary_obj = summary_stub),
    "Method only implemented for classes"
  )
})

test_that("init returns expected defaults for simple inputs", {
  data <- build_fixture_dataset()

  res <- lcc:::init(
    var.class   = NULL,
    weights.form = NULL,
    REML        = TRUE,
    qf          = 1,
    qr          = 0,
    pdmat       = "pdSymm()",
    dataset     = data,
    resp        = "H_mean",
    subject     = "Fruit",
    method      = "Method",
    time        = "Time",
    gs          = NULL,
    numCore     = 1
  )

  expect_equal(res$MethodREML, "REML")
  expect_identical(res$pdmat, nlme::pdSymm)
  expect_null(res$var.class)
})

test_that("init validates variance function pairing", {
  data <- build_fixture_dataset()

  expect_error(
    lcc:::init(
      var.class   = "varIdent()",
      weights.form = "time",
      REML        = FALSE,
      qf          = 1,
      qr          = 0,
      pdmat       = nlme::pdSymm,
      dataset     = data,
      resp        = "H_mean",
      subject     = "Fruit",
      method      = "Method",
      time        = "Time",
      gs          = NULL,
      numCore     = 1
    ),
    "Please specify 'weights.form' correctly for varExp class"
  )

  expect_error(
    lcc:::init(
      var.class   = "varExp()",
      weights.form = NULL,
      REML        = FALSE,
      qf          = 1,
      qr          = 0,
      pdmat       = nlme::pdSymm,
      dataset     = data,
      resp        = "H_mean",
      subject     = "Fruit",
      method      = "Method",
      time        = "Time",
      gs          = NULL,
      numCore     = 1
    ),
    "Please specify the 'weights.form' argument"
  )
})

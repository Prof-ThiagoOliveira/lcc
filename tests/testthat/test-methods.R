stopifnot(require("testthat"), require("lcc"))

context("Testing methods")

data(hue)
#-----------------------------------------------------------------------
# summary,  anova
#-----------------------------------------------------------------------
test_that("Test if interaction works",{
  expect_that(fm1 <- lcc(data = hue, subject = "Fruit",
                         resp = "H_mean", method = "Method",
                         time = "Time", qf = 2, qr = 2) ,is_a("lcc"))
  expect_that(summary(fm1), is_a("summary.lcc"))
  expect_that(summary(fm1,  type = "lcc"), is_a("summary.lcc"))
  expect_that(anova(fm1), is_a("anova.lcc"))
  expect_that(print(fm1), is_a("lcc"))
})
#=======================================================================
#-----------------------------------------------------------------------
# AIC,  BIC
#-----------------------------------------------------------------------
test_that("Testing methods",{
  expect_that(fmeint2<-lcc(data = hue, subject = "Fruit",
                         resp = "H_mean", method = "Method",
                         time = "Time", qf = 2, qr = 2),is_a("lcc"))
  expect_equal(AIC(fmeint2), AIC(fmeint2$model))
  expect_equal(BIC(fmeint2), BIC(fmeint2$model))
  expect_equal(residuals(fmeint2), residuals(fmeint2$model))
})
#=======================================================================
# Test class of methods
#=======================================================================
test_that("Test if interaction works",{
  expect_that(fm1 <- lcc(data = hue, subject = "Fruit",
                         resp = "H_mean", method = "Method",
                         time = "Time", qf = 2, qr = 2) ,is_a("lcc"))
  expect_equal(class(coef(fm1)),  c("coef.lcc", "ranef.lcc",
                                    "data.frame"))
  expect_equal(class(ranef(fm1)),  c("ranef.lcc", "data.frame"))
  expect_equal(class(getVarCov(fm1)),  c("random.effects", "VarCov"))
  expect_equal(class(residuals(fm1)),  c("numeric"))
  expect_equal(class(AIC(fm1)),  c("numeric"))
  expect_equal(class(BIC(fm1)),  c("numeric"))
  expect_equal(class(logLik(fm1)),  c("logLik"))
  expect_equal(class(anova(fm1)),  c("anova.lcc", "data.frame"))
})
#=======================================================================
# Test var-cov when changing time_lcc
#=======================================================================
test_that("Test if interaction works",{
  expect_that(fm1 <- lcc(data = hue, subject = "Fruit",
                         resp = "H_mean", method = "Method",
                         time = "Time", qf = 2, qr = 2) ,is_a("lcc"))
    expect_that(fm2 <- lcc(data = hue, subject = "Fruit",
                         resp = "H_mean", method = "Method",
                         time = "Time", qf = 2, qr = 2,
                         time_lcc = list(from = 0, to = 10, n = 30)),
                         is_a("lcc"))
  expect_equal(getVarCov(fm1), getVarCov(fm2))
  expect_equal(anova(fm1), anova(fm2))
  expect_length(fm1$Summary.lcc$fitted[, 1],
                n = length(unique(hue$Time)))
  expect_length(fm2$Summary.lcc$fitted[, 1],
                n = length(unique(c(hue$Time,
                                    seq(0, 10, length.out = 30)))))
})
#=======================================================================
test_that("fixture-based lcc object exercises S3 helpers", {
  fit <- build_test_lcc(components = TRUE)

  expect_true(is.lcc(fit))
  expect_false(is.lcc(list()))

  fitted_lcc <- testthat::capture_output(res_lcc <- fitted(fit, type = "lcc"))
  expect_s3_class(res_lcc, "data.frame")
  expect_true("fitted.LCC" %in% names(res_lcc))

  fitted_lpc <- testthat::capture_output(res_lpc <- fitted(fit, type = "lpc"))
  expect_s3_class(res_lpc, "data.frame")
  expect_true("fitted.LPC" %in% names(res_lpc))

  fitted_la <- testthat::capture_output(res_la <- fitted(fit, type = "la"))
  expect_s3_class(res_la, "data.frame")
  expect_true("fitted.LA" %in% names(res_la))

  expect_error(fitted(fit, type = "invalid"), "Available 'type'")

  printed <- testthat::capture_output(print_return <- print(fit))
  expect_true(nzchar(printed))
  expect_identical(print_return, fit)

  summary_model <- testthat::capture_output(sum_model <- summary(fit, type = "model"))
  expect_s3_class(sum_model, "summary.lcc")

  summary_metrics <- testthat::capture_output(sum_metrics <- summary(fit, type = "lcc"))
  expect_s3_class(sum_metrics, "summary.lcc")
  expect_true("fitted" %in% names(sum_metrics))

  expect_error(summary(fit, type = "unknown"), "Available 'type' are 'lcc' or 'model'")
})
#=======================================================================

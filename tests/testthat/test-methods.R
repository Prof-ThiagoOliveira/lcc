stopifnot(require("testthat"), require("lcc"))

context("Testing methods")

data(hue)
#-----------------------------------------------------------------------
# summary,  anova
#-----------------------------------------------------------------------
test_that("Test if interaction works",{
  expect_that(fm1 <- lcc(dataset = hue, subject = "Fruit",
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
  expect_that(fmeint2<-lcc(dataset = hue, subject = "Fruit",
                         resp = "H_mean", method = "Method",
                         time = "Time", qf = 2, qr = 2),is_a("lcc"))
  expect_equal(AIC(fmeint2), AIC(fmeint2$model))
  expect_equal(BIC(fmeint2), BIC(fmeint2$model))
  expect_equal(residuals(fmeint2), residuals(fmeint2$model))
})
#=======================================================================

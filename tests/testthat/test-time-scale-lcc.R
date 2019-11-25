stopifnot(require("testthat"), require("lcc"))

context("scale handling in time variable")

test_that("changes in the time scale does not affect lcc estimates", {
  data(hue)
  hue$Time2<-hue$Time+10
  hue$Time3<-scale(hue$Time)
  #---------------------------------------------------------------------
  expect_that(fm1<-lcc(dataset = hue, subject = "Fruit",
                       resp = "H_mean", method = "Method",
                       time = "Time", qf = 2, qr = 2,
                       components = TRUE), is_a("lcc"))
  #---------------------------------------------------------------------
  expect_that(fm2<-lcc(dataset = hue, subject = "Fruit",
                       resp = "H_mean", method = "Method",
                       time = "Time2", qf = 2, qr = 2,
                       components = TRUE), is_a("lcc"))
  #---------------------------------------------------------------------
  expect_that(fm3<-lcc(dataset = hue, subject = "Fruit",
                       resp = "H_mean", method = "Method",
                       time = "Time3", qf = 2, qr = 2,
                       components = TRUE), is_a("lcc"))
  #---------------------------------------------------------------------
  expect_equal(fm1$Summary.lcc$fitted[,-1],
               fm2$Summary.lcc$fitted[,-1],tolerance = 1e-4)
  expect_equal(fm1$Summary.lcc$fitted[,-1],
               fm3$Summary.lcc$fitted[,-1],tolerance = 1e-4)
})

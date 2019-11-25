stopifnot(require("testthat"), require("lcc"))

context("running regular sequence")

test_that("runing regular sequence with experimental time from 0 to 10",{
  time<-seq(0,10)
  expected<-unique(sort(c(time,seq(0,10,length.out = 30))))
  expect_equal(expected, time_lcc(time = time, from = 0, to = 10, n = 30))
})

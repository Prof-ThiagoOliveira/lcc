stopifnot(require("testthat"), require("lcc"))
context("running sampled CCC")

test_that("runing sampled CCC",{
  require(MASS)
  S<-matrix(c(2.0, 0.3,
              0.3, 1.8), byrow = TRUE, ncol=2)
  data<-mvrnorm(n=100, mu=c(10,9.5), Sigma=S)
  ccc_lin<-2*cov(data[,1],data[,2])/(var(data[,1])+var(data[,2])+(mean(data[,1])-mean(data[,2]))^2)
  ccc_lcc <- lcc:::CCC(data[,1],data[,2])
  expect_equal(ccc_lin, ccc_lcc)
})

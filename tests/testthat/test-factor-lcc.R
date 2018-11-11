library("testthat")
library("lcc")

context("factor handling in grouping variables")

test_that("factors", {
data(hue)
hue2<-transform(hue, "Method"=as.factor(as.numeric(hue$Method)),
  "Fruit"=as.factor(as.numeric(hue$Fruit)))
    
expect_that(fm1<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, components = TRUE), is_a("lcc"))

expect_that(fm2<-lcc(dataset = hue2, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, components = TRUE), is_a("lcc"))
expect_equivalent(fm1$Summary.lcc,fm2$Summary.lcc)
})

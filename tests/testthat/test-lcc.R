stopifnot(require("testthat"), require("lcc"))

context("fitting lcc models")

data(hue)

test_that("Name of variables", {
  expect_error(lcc(dataset = ddd, subject = "Fruit", resp = "H_mean",
    method = "Method", time = "Time", qf = 2, qr = 2),"object 'ddd' not found")
  expect_error(lcc(dataset = hue, subject = "AAAAA", resp = "H_mean",
    method = "Method", time = "Time", qf = 2, qr = 2),"Please, verify the name of 'resp', 'subject', 'method', and 'time' variables")
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "AAAAA",
    method = "Method", time = "Time", qf = 2, qr = 2),"Please, verify the name of 'resp', 'subject', 'method', and 'time' variables")
  expect_error(lcc(dataset = hue, subject = "Fruits", resp = "H_mean",
    method = "AAAAA", time = "Time", qf = 2, qr = 2),"Please, verify the name of 'resp', 'subject', 'method', and 'time' variables")
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean",
    method = "Method", time = "AAAAA", qf = 2, qr = 2),"Please, verify the name of 'resp', 'subject', 'method', and 'time' variables")
  })

test_that("If qr>qf should result in error",{
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 3),"'qr' should be less or equal 'qf'")
 })

test_that("pdmat" ,{
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, pdmat = pdIdent()),"Do not include brackets after the pdmat function, e.g. pdSymm()")
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, pdmat = AAAA),"object 'AAAA' not found")

# compatibility names in pdmat
  expect_that(fm1<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, pdmat = "pdIdent()"), is_a("lcc"))
  
  expect_that(fm2<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, pdmat = pdIdent), is_a("lcc"))
  
  expect_that(fm3<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, pdmat = "pdIdent"), is_a("lcc"))
  
  expect_equivalent(fm1,fm2)
  expect_equivalent(fm1,fm3)
  expect_equivalent(fm2,fm3)
})



test_that("var.class and weights.form" ,{
  # without declare weights.forms
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, var.class = varIdent),"Please specify the 'weights.form' argument.")
  
  # wrong name for var.class
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, var.class = varIdent(), weights.form = "method"),"Do not include brackets after the var.class function, e.g. varExp()")
  
  # wrong name for weights.form
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, var.class = varIdent, weights.form = "AAAA"),"The weights.form argument are \"time\", \"method\", \"time.ident\", or \"both\".")
  
  # wrong name for weights.form
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, var.class = varIdent, weights.form = "both"),"Please specify the 'weight.form' correctly for varIdent class")
  })


context("running longitudinal concordance correlation")
context("lcc")

test_that("Acessing an undefined variable should result in error",{
  data(hue)
  expect_error(
    {
      fm1<-lcc(dataset = hue, subject = "Fruits", resp = "H_mean",
        method = "Method", time = "Time", qf = 2, qr = 2)   
    },
    "Please, verify the name of 'resp', 'subject', 'method', and 'time' variables"
  )
  expect_error(
    {
      fm1<-lcc(dataset = hue, subject = "Fruit", resp = "H_means",
        method = "Method", time = "Time", qf = 2, qr = 2)   
    },
    "Please, verify the name of 'resp', 'subject', 'method', and 'time' variables"
  )
  expect_error(
    {
      fm1<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean",
        method = "Methods", time = "Time", qf = 2, qr = 2)   
    },
    "Please, verify the name of 'resp', 'subject', 'method', and 'time' variables"
  )
  expect_error(
    {
      fm1<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean",
        method = "Method", time = "Times", qf = 2, qr = 2)   
    },
    "Please, verify the name of 'resp', 'subject', 'method', and 'time' variables"
  )
  })

test_that("If qr>qf should result in error",{
  data(hue)
  expect_error(
    {
      fm1<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean",
        method = "Method", time = "Time", qf = 2, qr = 3)   
    },
    "'qr' should be less or equal 'qf'"
  )
})

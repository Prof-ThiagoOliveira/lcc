context("running longitudinal concordance correlation")

test_that("Acessing an undefined variable should result in error",{
  data(hue)
  expect_error(
    {
      fm1<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean",
        method = "Methods", time = "Time", qf = 2, qr = 2)   
    },
    "Please, verify the name of 'resp', 'subject', 'method', and 'time' variables"
  )
  })

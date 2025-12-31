test_that("time_lcc respects input bounds", {
  time <- c(1, 3, 5)
  seq_out <- time_lcc(time = time, from = 0, to = 6, n = 5)
  expect_equal(seq_out[1], min(time))
  expect_equal(seq_out[length(seq_out)], max(time))
  expect_true(all(time %in% seq_out))
})

test_that("time_lcc returns sorted unique grid", {
  time <- c(2, 4, 4, 7)
  seq_out <- time_lcc(time = time, from = 2, to = 7, n = 3)
  expect_false(is.unsorted(seq_out))
  expect_equal(seq_out, sort(unique(c(seq(2, 7, length.out = 3), time))))
})

test_that("time_lcc keeps regular sequence when already aligned", {
  time <- seq(0, 1, length.out = 5)
  seq_out <- time_lcc(time = time, from = 0, to = 1, n = 5)
  expect_equal(seq_out, time)
})

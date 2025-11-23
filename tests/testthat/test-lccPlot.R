# Testing lcc plot
stopifnot(require("testthat"), require("lcc"))

## Simulate dataset
Data <- function(N, time) {
  Beta01 <- 120 # Method A and Region A
  Beta02 <- 130 # Method A and Region B
  Beta03 <- 120 # Method B and Region A
  Beta04 <- 105 # Method B and Region B
  Beta11 <- -3   # Method A and Region A
  Beta12 <- -3.5 # Method A and Region B
  Beta13 <- -3.01# Method B and Region A
  Beta14 <- -2.5 # Method B and Region B
  sigmab0 <- 10
  sigmab1 <- 0.5
  corb01  <- 0.8
  covb0b1 <- corb01 * sqrt(sigmab0 * sigmab1)
  error   <- c(0.5, 0.5)
  
  par <- data.frame(
    Par       = c(Beta01, Beta02, Beta03, Beta04,
                  Beta11, Beta12, Beta13, Beta14,
                  sigmab0, sigmab1, corb01, error[1], error[2]),
    Parameters = c("Beta01", "Beta02", "Beta03", "Beta04",
                   "Beta11", "Beta12", "Beta13", "Beta14",
                   "sigmab0", "sigmab1", "corb01", "error[1]", "error[2]")
  )
  
  VAR1 <- function(time) {
    sigmab0 + sigmab1 * time^2 + 2 * covb0b1 * time + error[1]
  }
  
  VAR2 <- function(time) {
    sigmab0 + sigmab1 * time^2 + 2 * covb0b1 * time + error[2]
  }
  
  COV <- function(time) {
    sigmab0 + sigmab1 * time^2 + 2 * covb0b1 * time
  }
  
  E1 <- function(time) sum(c(1, time) * c(Beta01, Beta11))
  E2 <- function(time) sum(c(1, time) * c(Beta02, Beta12))
  E3 <- function(time) sum(c(1, time) * c(Beta03, Beta13))
  E4 <- function(time) sum(c(1, time) * c(Beta04, Beta14))
  
  mu1T <- Vectorize(function(t) E1(t), "t")
  mu2T <- Vectorize(function(t) E2(t), "t")
  mu3T <- Vectorize(function(t) E3(t), "t")
  mu4T <- Vectorize(function(t) E4(t), "t")
  
  S12T <- Vectorize(function(t) E1(t) - E2(t), "t")
  S13T <- Vectorize(function(t) E1(t) - E3(t), "t")
  S14T <- Vectorize(function(t) E1(t) - E4(t), "t")
  
  LCC12T <- Vectorize(function(t) 2 * COV(t) /
                        (VAR1(t) + VAR1(t) + (E1(t) - E2(t))^2), "t")
  LCC13T <- Vectorize(function(t) 2 * COV(t) /
                        (VAR2(t) + VAR2(t) + (E1(t) - E3(t))^2), "t")
  LCC14T <- Vectorize(function(t) 2 * COV(t) /
                        (VAR1(t) + VAR2(t) + (E1(t) - E4(t))^2), "t")
  
  Fruit  <- gl(N, length(time), 4 * N * length(time))
  Method <- gl(4, k = length(time) * N)
  Time   <- rep(time, N * 4)
  
  E1y <- Vectorize(function(t) E1(t), "t")
  E2y <- Vectorize(function(t) E2(t), "t")
  E3y <- Vectorize(function(t) E3(t), "t")
  E4y <- Vectorize(function(t) E4(t), "t")
  
  Evalue <- c(rep(E1y(time), N),
              rep(E2y(time), N),
              rep(E3y(time), N),
              rep(E4y(time), N))
  
  require(MASS)
  bi <- MASS::mvrnorm(
    N,
    mu    = c(0, 0),
    Sigma = matrix(c(sigmab0, covb0b1,
                     covb0b1, sigmab1),
                   byrow = TRUE, ncol = 2),
    empirical = FALSE
  )
  
  Zi  <- model.matrix(~ time)
  Zb. <- vector("list", N)
  for (i in seq_len(N)) {
    Zb.[[i]] <- Zi %*% bi[i, ]
  }
  Zb <- rep(unlist(Zb.), 4)
  
  residual <- c(
    rnorm(length(Zb) / 2, mean = 0, sd = sqrt(error[1])),
    rnorm(length(Zb) / 2, mean = 0, sd = sqrt(error[2]))
  )
  
  Response <- Evalue + Zb + residual
  
  dataset <- data.frame(Fruit, Method, Response, Time)
  
  Time2 <- time
  list(
    data   = dataset,
    par    = par,
    LCC12T = LCC12T(Time2),
    LCC13T = LCC13T(Time2),
    LCC14T = LCC14T(Time2),
    mu1T   = mu1T(Time2),
    mu2T   = mu2T(Time2),
    mu3T   = mu3T(Time2),
    mu4T   = mu4T(Time2),
    S12T   = S12T(Time2),
    S13T   = S13T(Time2),
    S14T   = S14T(Time2)
  )
}

set.seed(5925670)
dataset <- Data(N = 30, time = seq(0, 15, 1))

## ----------------------------------------------------------------------
## Basic error check
## ----------------------------------------------------------------------

test_that("Object does not inherit from class lcc", {
  aaa <- lm(rnorm(100, 10, 1) ~ rnorm(100, 50, 3))
  expect_error(lccPlot(aaa), "Object must inherit from class \"lcc\"")
})

data(hue)

## ----------------------------------------------------------------------
## LCC, LPC and LA plots (no CI)
## ----------------------------------------------------------------------

test_that("LCC, LPC and LA plot test", {
  # Two methods
  expect_that(
    fme1 <- lcc(
      data = hue, subject = "Fruit",
      resp = "H_mean", method = "Method",
      time = "Time", qf = 1, qr = 1
    ),
    is_a("lcc")
  )
  expect_silent(p1 <- suppressWarnings(lccPlot(fme1)))
  expect_s3_class(p1, "ggplot")
  
  # Components TRUE
  expect_that(
    fme2 <- lcc(
      data = hue, subject = "Fruit",
      resp = "H_mean", method = "Method",
      time = "Time", qf = 1, qr = 1,
      components = TRUE
    ),
    is_a("lcc")
  )
  expect_silent(p2 <- suppressWarnings(lccPlot(fme2)))
  expect_s3_class(p2, "ggplot")
  
  # More than two methods
  expect_that(
    fme3 <- lcc(
      data = dataset$data, subject = "Fruit",
      resp = "Response", method = "Method",
      time = "Time", qf = 1, qr = 1
    ),
    is_a("lcc")
  )
  expect_silent(p3 <- suppressWarnings(lccPlot(fme3)))
  expect_s3_class(p3, "ggplot")
  
  # Components TRUE
  expect_that(
    fme4 <- lcc(
      data = dataset$data, subject = "Fruit",
      resp = "Response", method = "Method",
      time = "Time", qf = 1, qr = 1,
      components = TRUE
    ),
    is_a("lcc")
  )
  expect_silent(p4 <- suppressWarnings(lccPlot(fme4)))
  expect_s3_class(p4, "ggplot")
})

## ----------------------------------------------------------------------
## Confidence interval plots
## ----------------------------------------------------------------------

test_that("Confidence intervals plot", {
  # Two methods
  expect_that(
    fme5 <- lcc(
      data = hue, subject = "Fruit",
      resp = "H_mean", method = "Method",
      time = "Time", qf = 1, qr = 1,
      ci = TRUE, nboot = 100
    ),
    is_a("lcc")
  )
  expect_silent(p5 <- suppressWarnings(lccPlot(fme5)))
  expect_s3_class(p5, "ggplot")
  
  # Components TRUE
  expect_that(
    fme6 <- lcc(
      data = hue, subject = "Fruit",
      resp = "H_mean", method = "Method",
      time = "Time", qf = 1, qr = 1,
      components = TRUE, ci = TRUE, nboot = 100
    ),
    is_a("lcc")
  )
  expect_silent(p6 <- suppressWarnings(lccPlot(fme6)))
  expect_s3_class(p6, "ggplot")
  
  # More than two methods
  expect_that(
    fme7 <- lcc(
      data = dataset$data, subject = "Fruit",
      resp = "Response", method = "Method",
      time = "Time", qf = 1, qr = 1,
      ci = TRUE, nboot = 100
    ),
    is_a("lcc")
  )
  expect_silent(p7 <- suppressWarnings(lccPlot(fme7)))
  expect_s3_class(p7, "ggplot")
  
  # Components TRUE
  expect_that(
    fme8 <- lcc(
      data = dataset$data, subject = "Fruit",
      resp = "Response", method = "Method",
      time = "Time", qf = 1, qr = 1,
      components = TRUE, ci = TRUE, nboot = 100
    ),
    is_a("lcc")
  )
  expect_silent(p8 <- suppressWarnings(lccPlot(fme8)))
  expect_s3_class(p8, "ggplot")
})

## ----------------------------------------------------------------------
## Labels, shape and colour
## ----------------------------------------------------------------------

test_that("labels, shape and colour", {
  expect_that(
    fm <- lcc(
      data = dataset$data, subject = "Fruit",
      resp = "Response", method = "Method",
      time = "Time", qf = 1, qr = 1,
      components = TRUE
    ),
    is_a("lcc")
  )
  
  ## LCC
  expect_silent(p_lcc <- suppressWarnings(lccPlot(
    fm, type = "lcc",
    control = list(
      shape = 2, colour = "red", size = 2,
      xlab = "Time (hours)", ylab = "Longitudinal CC"
    )
  )))
  expect_s3_class(p_lcc, "ggplot")
  expect_identical(p_lcc$labels$x, "Time (hours)")
  expect_identical(p_lcc$labels$y, "Longitudinal CC")
  expect_identical(p_lcc$layers[[1]]$aes_params$colour, "red")
  expect_identical(p_lcc$layers[[1]]$aes_params$linewidth, 2)
  expect_identical(p_lcc$layers[[2]]$aes_params$shape, 2)
  
  ## LPC
  expect_silent(p_lpc <- suppressWarnings(lccPlot(
    fm, type = "lpc",
    control = list(
      shape = 2, colour = "red", size = 2,
      xlab = "Time (hours)", ylab = "Longitudinal PC"
    )
  )))
  expect_s3_class(p_lpc, "ggplot")
  expect_identical(p_lpc$labels$x, "Time (hours)")
  expect_identical(p_lpc$labels$y, "Longitudinal PC")
  
  ## LA
  expect_silent(p_la <- suppressWarnings(lccPlot(
    fm, type = "la",
    control = list(
      shape = 2, colour = "red", size = 2,
      xlab = "Time (hours)", ylab = "Longitudinal Accuracy"
    )
  )))
  expect_s3_class(p_la, "ggplot")
  expect_identical(p_la$labels$x, "Time (hours)")
  expect_identical(p_la$labels$y, "Longitudinal Accuracy")
})

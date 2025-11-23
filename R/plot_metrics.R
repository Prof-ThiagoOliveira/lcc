##' @title Internal Function to Estimate the Sampled Concordance
##'   Correlation Coefficient.
##'
##' @description This function is internally called to estimate
##'   the sampled concordance correlation coefficient.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @importFrom stats cor cov
##'
##' @keywords internal
CCC_lin <- function(dataset, resp, subject, method, time) {
  # We assume here that `dataset` has already been prepared by dataBuilder
  selectedData <- subset(dataset, select = c(resp, method, time, subject))
  dataByMethod <- split(selectedData, selectedData[[method]])
  methodLevels <- levels(selectedData[[method]])
  
  calculateCCC_fast <- function(Y1, Y2, time) {
    n          <- length(time)
    time_fac   <- as.factor(time)
    idx_by_t   <- split(seq_len(n), time_fac)
    
    Y1_full    <- rep(Y1, length.out = n)
    Y2_full    <- rep(Y2, length.out = n)
    
    ccc_vec <- vapply(
      idx_by_t,
      function(idx) {
        y1 <- Y1_full[idx]
        y2 <- Y2_full[idx]
        m1 <- mean(y1)
        m2 <- mean(y2)
        s1 <- var(y1)
        s2 <- var(y2)
        s12 <- cov(y1, y2)
        2 * s12 / (s1 + s2 + (m1 - m2)^2)
      },
      numeric(1L)
    )
    data.frame(V1 = unname(ccc_vec))
  }
  
  CCC.Lin <- lapply(
    seq(2L, length(methodLevels)),
    function(i) {
      calculateCCC_fast(
        Y1   = dataByMethod[[1L]][[resp]],
        Y2   = dataByMethod[[i]][[resp]],
        time = selectedData[[time]]
      )
    }
  )
  
  CCC.Lin
}


#' @title Estimate Sampled Pearson Correlation
#'
#' @description Internally called function to estimate the sampled Pearson correlation.
#'
#' @usage NULL
#'
#' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
#'
#' @importFrom stats cor
#'
#' @keywords internal
Pearson <- function(dataset, resp, subject, method, time) {
  selectedData <- subset(dataset, select = c(resp, method, time, subject))
  dataByMethod <- split(selectedData, selectedData[[method]])
  methodLevels <- levels(selectedData[[method]])
  
  calculateCorrelation_fast <- function(Y1, Y2, time) {
    n        <- length(time)
    time_fac <- as.factor(time)
    idx_by_t <- split(seq_len(n), time_fac)
    
    Y1_full  <- rep(Y1, length.out = n)
    Y2_full  <- rep(Y2, length.out = n)
    
    cor_vec <- vapply(
      idx_by_t,
      function(idx) cor(Y1_full[idx], Y2_full[idx]),
      numeric(1L)
    )
    
    data.frame(V1 = unname(cor_vec))
  }
  
  pearsonResults <- vector("list", length(methodLevels) - 1L)
  for (i in seq(2L, length(methodLevels))) {
    pearsonResults[[i - 1L]] <- calculateCorrelation_fast(
      Y1   = dataByMethod[[1L]][[resp]],
      Y2   = dataByMethod[[i]][[resp]],
      time = selectedData[[time]]
    )
  }
  
  pearsonResults
}
#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: plot_metrics.R                                                #
# Contains: CCC_lin, Pearson                                          #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 23/11/2025                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

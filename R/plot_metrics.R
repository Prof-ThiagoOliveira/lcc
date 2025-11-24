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
  df <- dataset[, c(resp, method, time, subject)]
  names(df) <- c("resp", "method", "time", "subject")
  method_levels <- levels(df$method)
  ref_level     <- method_levels[1L]
  
  calculateCCC_fast <- function(df_time, target_level) {
    wide <- stats::reshape(
      df_time[, c("subject", "method", "resp")],
      idvar   = "subject",
      timevar = "method",
      direction = "wide"
    )
    col_ref <- paste0("resp.", ref_level)
    col_tgt <- paste0("resp.", target_level)
    if (!all(c(col_ref, col_tgt) %in% names(wide))) return(NA_real_)
    y1 <- wide[[col_ref]]
    y2 <- wide[[col_tgt]]
    keep <- stats::complete.cases(y1, y2)
    if (sum(keep) < 2L) return(NA_real_)
    m1  <- mean(y1[keep])
    m2  <- mean(y2[keep])
    s1  <- stats::var(y1[keep])
    s2  <- stats::var(y2[keep])
    s12 <- stats::cov(y1[keep], y2[keep])
    2 * s12 / (s1 + s2 + (m1 - m2)^2)
  }
  
  split_time <- split(df, df$time)
  
  lapply(
    seq(2L, length(method_levels)),
    function(i) {
      target <- method_levels[i]
      vals <- vapply(split_time, calculateCCC_fast, numeric(1L), target_level = target)
      data.frame(V1 = unname(vals))
    }
  )
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
  df <- dataset[, c(resp, method, time, subject)]
  names(df) <- c("resp", "method", "time", "subject")
  method_levels <- levels(df$method)
  ref_level     <- method_levels[1L]
  split_time    <- split(df, df$time)
  
  calculateCorrelation_fast <- function(df_time, target_level) {
    wide <- stats::reshape(
      df_time[, c("subject", "method", "resp")],
      idvar   = "subject",
      timevar = "method",
      direction = "wide"
    )
    col_ref <- paste0("resp.", ref_level)
    col_tgt <- paste0("resp.", target_level)
    if (!all(c(col_ref, col_tgt) %in% names(wide))) return(NA_real_)
    y1 <- wide[[col_ref]]
    y2 <- wide[[col_tgt]]
    keep <- stats::complete.cases(y1, y2)
    if (sum(keep) < 2L) return(NA_real_)
    stats::cor(y1[keep], y2[keep])
  }
  
  lapply(
    seq(2L, length(method_levels)),
    function(i) {
      target <- method_levels[i]
      vals <- vapply(split_time, calculateCorrelation_fast, numeric(1L), target_level = target)
      data.frame(V1 = unname(vals))
    }
  )
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

#' @importFrom stats reshape

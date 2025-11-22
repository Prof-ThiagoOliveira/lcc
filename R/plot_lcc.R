#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: plot_lcc.R                                                    #
# Contains: CCC_lin, plot_lcc                                         #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

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


##' @title Internal Function to Prepare the \code{plotBuilder_lcc} Function
##'
##' @description This function is internally called to prepare
##'   the \code{\link[lcc]{plotBuilder_lcc}} function.
##'
##' @usage NULL
##'##' @title Internal Function to Prepare the \code{plotBuilder_la} Function
##'
##' @description This function is internally called to prepare
##'   the \code{\link[lcc]{plotBuilder_la}} function.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @keywords internal
plot_la <- function(Cb, ENV.Cb, tk.plot, tk.plot2, ldb, model, ci, arg, ...) {
  CCC <- CCC_lin(dataset = model$data, resp = "resp", subject = "subject", method = "method", time = "time")
  Pearson <- Pearson(dataset = model$data, resp = "resp", subject = "subject", method = "method", time = "time")
  
  if (ci) {
    plotBuilder_la(CCC = CCC, Pearson = Pearson, ENV.Cb = ENV.Cb, 
                   tk.plot = tk.plot, tk.plot2 = tk.plot2, ldb = ldb, Cb = Cb, 
                   model = model, ci = TRUE, arg = arg, ...)
  } else {
    plotBuilder_la(CCC = CCC, Pearson = Pearson, tk.plot = tk.plot, 
                   tk.plot2 = tk.plot2, ldb = ldb, Cb = Cb, model = model, 
                   ci = FALSE, arg = arg, ...)
  }
}

##'
##' @keywords internal
plot_lcc <- function(rho, ENV.LCC, tk.plot, tk.plot2, ldb, model, ci, arg, ...) {
  CCC <- CCC_lin(dataset = model$data, resp = "resp", subject = "subject", 
                 method = "method", time = "time")
  
  if (ci) {
    plotBuilder_lcc(rho = rho, ENV.LCC = ENV.LCC, tk.plot = tk.plot, 
                    tk.plot2 = tk.plot2, ldb = ldb, CCC = CCC, model = model,
                    ci = TRUE, arg = arg, ...)
  } else {
    plotBuilder_lcc(rho = rho, tk.plot = tk.plot, tk.plot2 = tk.plot2, ldb = ldb,
                    CCC = CCC, model = model, ci = FALSE, arg = arg, ...)
  }
}


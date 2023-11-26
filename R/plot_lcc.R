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
  selectedData <- subset(dataset, select = c(resp, method, time, subject))
  dataByMethod <- split(selectedData, selectedData$method)
  
  calculateCCC <- function(Y1, Y2, time) {
    dataFrame <- data.frame(Y1, Y2, time)
    means <- with(dataFrame, sapply(list(Y1, Y2), function(Y) tapply(Y, time, mean)))
    variances <- with(dataFrame, sapply(list(Y1, Y2), function(Y) tapply(Y, time, var)))
    covariances <- by(dataFrame[, 1:2], dataFrame$time, function(x) cov(x$Y1, x$Y2))
    
    dataLin <- data.frame(
      time = unique(time),
      M1 = means[[1]],
      M2 = means[[2]],
      S1 = variances[[1]],
      S2 = variances[[2]],
      S12 = sapply(covariances, '[', 1)
    )
    
    CCC.results <- by(dataLin, dataLin$time, function(x) {
      2 * x$S12 / (x$S1 + x$S2 + (x$M1 - x$M2)^2)
    })
    
    CCC.results <- sapply(CCC.results, function(x) unlist(x))
    return(as.data.frame(as.matrix(CCC.results)))
  }
  
  CCC.Lin <- lapply(seq(2, length(levels(selectedData$method))), function(i) {
    calculateCCC(Y1 = dataByMethod[[1]]$resp, Y2 = dataByMethod[[i]]$resp, time = selectedData$time)
  })
  
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


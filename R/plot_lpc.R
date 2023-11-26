#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: plot_lpc.R                                                    #
# Contains: Pearson, plot_lpc                                         #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

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
  dataByMethod <- split(selectedData, selectedData$method)
  
  calculateCorrelation <- function(Y1, Y2, time) {
    dataFrame <- data.frame(Y1, Y2, time)
    correlations <- by(dataFrame[, 1:2], dataFrame$time, 
                       function(x) cor(x$Y1, x$Y2))
    as.data.frame(as.matrix(correlations))
  }
  
  pearsonResults <- list()
  methodLevels <- levels(selectedData$method)
  
  for (i in 2:length(methodLevels)) {
    pearsonResults[[i - 1]] <- calculateCorrelation(
      Y1 = dataByMethod[[1]]$resp, 
      Y2 = dataByMethod[[i]]$resp, 
      time = selectedData$time
    )
  }
  
  pearsonResults
}

#' @title Prepare Plot for LPC Function
#'
#' @description Internally called function to prepare data for the `plotBuilder_lpc` function.
#'
#' @usage NULL
#'
#' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
#'
#' @keywords internal
plot_lpc <- function(LPC, ENV.LPC, tk.plot, tk.plot2, ldb, model,
                              ci, arg, ...) {
  Pearson <- Pearson(
    dataset = model$data, 
    resp = "resp", 
    subject = "subject", 
    method = "method", 
    time = "time"
  )
  
  if (ci) {
    plotBuilder_lpc(
      LPC = LPC, ENV.LPC = ENV.LPC, tk.plot = tk.plot,
      tk.plot2 = tk.plot2, ldb = ldb, Pearson = Pearson,
      model = model, ci = TRUE, arg = arg, ...
    )
  } else {
    plotBuilder_lpc(
      LPC = LPC, tk.plot = tk.plot, tk.plot2 = tk.plot2, ldb = ldb,
      Pearson = Pearson, model = model, ci = FALSE, arg = arg, ...
    )
  }
}

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

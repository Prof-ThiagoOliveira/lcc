#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: lpcWrapper.R                                                  #
# Contains: lpcWrapper function                                       #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal Function to Prepare the \code{lpcBuilder} Function
##'
##' @description This function is internally called to prepare
##'   the \code{\link[lcc]{lpcBuilder}} function.
##'
##' @usage NULL
##'
##' @details Returns a vector or list containing the longitudinal
##'   Pearson correlation estimates.
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @keywords internal
lpcWrapper <- function(model, q_f, tk, n.delta) {
  covarianceMatrix <- getVarCov(model)
  q_r <- dim(covarianceMatrix)[1] - 1
  deltaValues <- getDelta(model = model)
  delta <- deltaValues$delta
  deltal <- deltaValues$deltal
  g <- deltaValues$g
  varianceEpsilon <- model$sigma^2
  
  rhoPearson <- lpcBuilder(
    G = covarianceMatrix, tk = tk, q_r = q_r, q_f = q_f, g = g, 
    sig2_epsilon = varianceEpsilon, delta = delta, deltal = deltal, model = model
  )
  
  rhoPearson[[n.delta]]
}


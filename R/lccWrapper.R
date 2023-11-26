#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: lccWrapper.R                                                  #
# Contains: lccWrapper function                                       #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal Function to Prepare the \code{lccBuilder} Function
##'
##' @description Internally called function to prepare the \code{\link[lcc]{lccBuilder}} function.
##'
##' @details Returns a vector or list containing the longitudinal concordance correlation estimates.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @keywords internal
lccWrapper <- function(model, q_f, tk, diffbeta, n.delta) {
  G <- getVarCov(model)
  q_r <- dim(G)[1] - 1
  deltas <- getDelta(model = model)
  sig2_epsilon <- model$sigma^2
  
  rho <- lccBuilder(
    G = G, diffbeta = diffbeta, tk = tk, q_r = q_r,
    q_f = q_f, g = deltas$g, sig2_epsilon = sig2_epsilon,
    delta = deltas$delta, deltal = deltas$deltal, model = model
  )
  
  if (length(rho) == 1 || sum(is.na(rho[[2]])) != 0) {
    rho[[1]]
  } else {
    rho[[n.delta]]
  }
}


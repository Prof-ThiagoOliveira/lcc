#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: lpcBuilder.R                                                  #
# Contains: lpcBuilder function                                       #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal Function to Estimate the Longitudinal Pearson Correlation
##'
##' @description Internally used function to estimate the longitudinal Pearson correlation (LPC).
##'
##' @details Returns a vector or list containing the longitudinal Pearson correlation estimates.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @keywords internal
lpcBuilder <- function(G, tk, q_r, q_f, g, sig2_epsilon, delta, deltal, model) {
  Tk_r <- sapply(0:q_r, function(x) tk^x)
  tGt <- diag(Tk_r %*% G %*% t(Tk_r))
  varcomp <- summary(model)
  var.f <- class(varcomp$modelStruct$varStruct)[1]
  rho.pearson <- list()
  
  calculateRhoPearson <- function(sig2, gd, gdl) {
    as.numeric(tGt / sqrt((tGt + sig2 * gd) * (tGt + sig2 * gdl)))
  }
  
  switch(var.f,
         "varIdent" = {
           gd <- g(delta)
           gdl <- g(deltal)
           rho.pearson <- lapply(gdl, function(gdli) calculateRhoPearson(sig2_epsilon, gd, gdli))
         },
         "varExp" = {
           if(attr(varcomp$modelStruct$varStruct, "formula") == ~time) {
             gd <- g(delta, tk)
             gdl <- g(deltal, tk)
             rho.pearson <- list(tGt / sqrt((tGt + sig2_epsilon * gd) * (tGt + sig2_epsilon * gdl)), NA)
           } else if(attr(varcomp$modelStruct$varStruct, "formula") == ~time | method) {
             gd <- g(delta, tk)
             gdl <- lapply(deltal, function(d) g(d, tk))
             rho.pearson <- lapply(gdl, function(gdli) calculateRhoPearson(sig2_epsilon, gd, gdli))
           } else {
             print("Method not implemented yet")
           }
         },
         "NULL" = {
           rho.pearson <- list(tGt / (tGt + sig2_epsilon), NA)
         }
  )
  
  rho.pearson
}


#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: lccBuilder.R                                                 #
# Contains: lccBuilder function                                       #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 18/06/2018                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal function to estimate the longitudinal concordance correlation.
##'
##' @description This is an internally called function used to estimate the longitudinal concordance correlation (LCC).
##'
##' @usage NULL
##'
##' @details returns a vector or list containing the longitudinal concordance correlation estimates.
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@usp.br}
##'
##' @keywords internal
lccBuilder <- function(G, diffbeta, tk, q_r, q_f, g, sig2_epsilon, delta, deltal, model) {
  Tk_r <- sapply(0:q_r, function(x) tk^x)
  Tk_f <- sapply(0:q_f, function(x) tk^x)
  S2 <- try((as.matrix(Tk_f[,c(1:(length(diffbeta)))]) %*% (diffbeta))^2, silent = TRUE)
  if(class(S2) == "try-error") {
    stop(print("Please, change the name of factor levels associated with method argument. Example:  levels M1, M2, M3, ..., Mn, with n=0, 1, 2, ..., N."), call. = FALSE)
  }
  tGt <- diag(Tk_r %*% G %*% t(Tk_r))
  varcomp <- summary(model)
  var.f <- class(varcomp$modelStruct$varStruct)[1]
  rho <- list()
  if(var.f == "varIdent") {
    gd <- g(delta)
    gdl <- g(deltal)
    ldb2 <- length(gdl)
    for(i in 1:ldb2){
      rho[[i]] <- as.numeric(
        tGt / (tGt + 0.5 * (sig2_epsilon * (gd + gdl[i]) + S2))
      )
    }
  }
  if(var.f == "varExp") {
    if(attr(varcomp$modelStruct$varStruct, "formula")==~time){
      gd <- g(delta,tk)
      gdl <- g(deltal,tk)
      rho <-list(tGt / (tGt + 0.5 * (sig2_epsilon * (gd + gdl) + S2)),NA)
    } else {if(attr(varcomp$modelStruct$varStruct, "formula")==~time | FacA){
      gd <- g(delta,tk)
      ldb2 <- length(deltal)
      gdl<-list()
      for(i in 1:ldb2){
        gdl[[i]] <- g(deltal[i],tk)
        rho[[i]] <- as.numeric(unlist(tGt / (tGt + 0.5 * (sig2_epsilon * (gd + gdl[[i]]) + S2))))
      }
    }else{print("method not implemented yet")}
    }
  }
  if(var.f == "NULL"){
    rho<-list(tGt / (tGt + sig2_epsilon + 0.5*S2),NA)
  }
  return(rho)
}

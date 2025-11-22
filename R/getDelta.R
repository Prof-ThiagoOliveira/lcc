#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: getDelta.R                                                    #
# Contains: getDelta function                                         #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################
##' @title Internal Function to Extract Variance Components Estimates.
##'
##' @usage NULL
##'
##' @description This is an internally called function used to extract
##'   variance components estimate of \eqn{\Sigma} matrix based on
##'   specified structure.
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@alumni.usp.br} and Rafael de Andrade Moral,
##'   \email{rafael_moral@@yahoo.com.br}
##'
##' @keywords internal
getDelta <- function(model) {
  varcomp   <- summary(model)
  varStruct <- varcomp$modelStruct$varStruct
  
  # default placeholders
  delta  <- 0
  deltal <- 0
  g_mode <- "const"   # which flavour of g() to use: "const", "ident", "exp"
  
  #-------------------------------------------------------------------
  # No variance structure
  #-------------------------------------------------------------------
  if (is.null(varStruct)) {
    # delta = deltal = 0, g(x) = 1
    g_mode <- "const"
    
  } else {
    # There *is* a variance structure: varIdent or varExp
    var.f <- class(varStruct)[1L]
    
    if (var.f == "varIdent") {
      # components2 = (1, lambda_2, ..., lambda_k) with lambda_j = exp(theta_j)^2
      components2 <- c(1, exp(as.numeric(varStruct))^2)
      delta       <- 1
      deltal      <- components2[-1L]
      g_mode      <- "ident"
      
    } else if (var.f == "varExp") {
      # distinguish ~time vs ~time | method
      form     <- attr(varStruct, "formula")
      form_chr <- gsub(" ", "", paste(deparse(form), collapse = ""))
      
      if (form_chr == "~time") {
        # single time-dependent parameter
        delta  <- varStruct
        deltal <- varStruct
        g_mode <- "exp"
        
      } else if (form_chr == "~time|method") {
        # first coefficient: time; others: method-specific weights
        components2 <- as.numeric(varStruct)
        delta       <- components2[1L]
        deltal      <- components2[-1L]
        g_mode      <- "exp"
        
      } else {
        stop("Method not implemented for the specified varStruct formula.",
             call. = FALSE)
      }
      
    } else {
      stop("Method only implemented for classes varIdent and varExp",
           call. = FALSE)
    }
  }
  
  #-------------------------------------------------------------------
  # Define g() once, according to mode
  #-------------------------------------------------------------------
  g <- switch(
    g_mode,
    "const" = function(x, tk = NULL) 1,
    "ident" = function(x, tk = NULL) x,
    "exp"   = function(x, tk) exp(2 * x * tk),
    # fallback (should not be used)
    function(x, tk = NULL) 1
  )
  
  list(delta = delta, deltal = deltal, g = g)
}

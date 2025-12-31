#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: variance_components.R                                         #
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
extract_random_effects_cov <- function(model, tolerance = 1e-8) {
  G_obj <- nlme::getVarCov(model, type = "random.effects")

  coerce_matrix <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    if (is.matrix(x)) {
      return(x)
    }
    if (inherits(x, "pdMat") || inherits(x, "VarCov")) {
      return(as.matrix(x))
    }
    if (is.numeric(x) && length(x) == 1L) {
      val <- as.numeric(x)
      return(matrix(val, nrow = 1L, ncol = 1L))
    }
    NULL
  }

  mats <- NULL
  if (is.list(G_obj)) {
    mats <- lapply(G_obj, coerce_matrix)
    mats <- Filter(Negate(is.null), mats)
  } else {
    mats <- list(coerce_matrix(G_obj))
  }

  if (!length(mats) || any(vapply(mats, is.null, logical(1L)))) {
    abort_internal("Failed to extract random-effects covariance matrix.")
  }

  G <- mats[[1L]]
  for (mat in mats[-1L]) {
    if (!all(dim(mat) == dim(G))) {
      abort_internal(
        "Random-effects covariance matrix dimensions differ across groups."
      )
    }
    if (max(abs(mat - G)) > tolerance) {
      warn_general(
        "Random-effects covariance differs across groups beyond tolerance; using first instance."
      )
      break
    }
  }

  if (!is.matrix(G) || !nrow(G) || !ncol(G)) {
    abort_internal("Failed to extract random-effects covariance matrix.")
  }
  if (nrow(G) != ncol(G)) {
    abort_internal(
      "Random-effects covariance matrix must be square; got {.val {nrow(G)}} x {.val {ncol(G)}}."
    )
  }

  list(G = G, n_re = nrow(G))
}

##' @keywords internal
getDelta <- function(model, summary_obj = NULL) {
  if (is.null(summary_obj)) {
    summary_obj <- summary(model)
  }
  varStruct <- summary_obj$modelStruct$varStruct
  
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
        abort_input("Method not implemented for the specified varStruct formula.")
      }
      
    } else {
      abort_input("Method only implemented for classes varIdent and varExp")
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

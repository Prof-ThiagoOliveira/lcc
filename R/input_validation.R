#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: input_validation.R                                            #
# Contains: init function                                             #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################


##' @title Internal Function to Prepare \code{\link[lcc]{lccModel}}
##'   Function
##'
##' @description This is an internally called function used to verify
##'   the specification of variance-covariance matrices and likelihood
##'   based method.
##'
##' @usage NULL
##' @return A named list with resolved arguments used internally by lccModel
##' @importFrom nlme lme
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @keywords internal
init <- function(var.class, weights.form, REML, qf, qr, pdmat, dataset,
                 resp, subject, method, time, gs, numCore) {
  resp    <- resp
  subject <- subject
  method  <- method
  time    <- time
  
  #-------------------------------------------------------------------
  # Rename columns to standard internal names: y, ind, FacA, time
  #-------------------------------------------------------------------
  Data <- data.frame(dataset)
  Data <- try(
    rename.vars(
      Data,
      from = c(resp, subject, method, time),
      to   = c("y", "ind", "FacA", "time"),
      info = FALSE
    ),
    TRUE
  )
  
  if (inherits(Data, "try-error")) {
    abort_input("Please, verify the name of 'resp', 'subject', 'method', and 'time' variables")
  }
  
  check_is_factor_col(Data, "ind", "subject")
  check_is_factor_col(Data, "FacA", "method")
  check_is_numeric_col(Data, "time", "time")
  check_is_numeric_col(Data, "y", "resp")
  
  #-------------------------------------------------------------------
  # Test / normalize pdmat
  #-------------------------------------------------------------------
  if (!is.function(pdmat)) {
    if (is.character(pdmat)) {
      if (substr(pdmat, nchar(pdmat) - 1L, nchar(pdmat)) == "()") {
        pdmat <- substr(pdmat, 1L, nchar(pdmat) - 2L)
      }
      pdmat <- get(pdmat)
    } else {
      abort_input("Do not include brackets after the pdmat function, e.g. pdSymm()")
    }
  }
  
  #-------------------------------------------------------------------
  # Test for var.class
  #-------------------------------------------------------------------
  if (is.null(var.class) == FALSE) {
    if (!is.function(var.class)) {
      if (is.character(var.class)) {
        if (substr(var.class, nchar(var.class) - 1L, nchar(var.class)) == "()") {
          var.class <- substr(var.class, 1L, nchar(var.class) - 2L)
        }
        var.class <- get(var.class)
      } else {
        abort_input("Do not include brackets after the var.class function, e.g. varExp()")
      }
    }
  }
  
  #-------------------------------------------------------------------
  # Test for var.class and weights.form
  # (NB: this intentionally calls var.class() as in the original code)
  #-------------------------------------------------------------------
  if (is.null(var.class) == FALSE) {
    vc <- class(var.class())[1L]
    
    if (is.null(weights.form)) {
      abort_input("Please specify the 'weights.form' argument.")
    }
    
    if (!weights.form %in% c("time", "method", "time.ident", "both")) {
      abort_input(
        "The weights.form argument are \"time\", \"method\", \"time.ident\", or \"both\"."
      )
    }
    
    if (weights.form == "time.ident" || weights.form == "method") {
      vc <- class(var.class())[1L]
      if (vc != "varIdent") {
        abort_input("Please specify 'weights.form' correctly for varIdent class")
      }
    }
    
    if (weights.form == "time" || weights.form == "both") {
      vc <- class(var.class())[1L]
      if (vc != "varExp") {
        abort_input("Please specify 'weights.form' correctly for varExp class")
      }
    }
    
    if (vc != "varIdent" && vc != "varExp") {
      abort_input(
        "Method only implemented for classes {.cls varIdent} and {.cls varExp}."
      )
    }
  }
  
  #-------------------------------------------------------------------
  # Test for REML
  #-------------------------------------------------------------------
  check_REML_flag(REML)
  MethodREML <- if (REML) "REML" else "ML"
  
  #-------------------------------------------------------------------
  # Test for qr and qf
  #-------------------------------------------------------------------
  check_polynomial_degrees(qf, qr)
  
  initList <- list(
    "MethodREML" = MethodREML,
    "pdmat"      = pdmat,
    "var.class"  = var.class
  )
  
  #-------------------------------------------------------------------
  # Test for gs
  #-------------------------------------------------------------------
  if (is.null(gs) == FALSE) {
    check_gs(gs, dataset = dataset, method = method)
  }
  
  #-------------------------------------------------------------------
  # Number of cores
  #-------------------------------------------------------------------
  check_num_core(numCore)
  
  initList
}

#-----------------------------------------------------------------------
# Rename vars
#-----------------------------------------------------------------------
##' @title Internal Function to Prepare \code{\link[lcc]{lccModel}}
##'   Function
##'
##' @description This is an internally called function used to verify
##'   the specification of variance-covariance matrices and likelihood
##'   based method.
##'
##' @usage NULL
##'
##' @author Code by Don MacQueen
##'
##' @keywords internal
rename.vars <- function(data, from = "", to = "", info = TRUE) {
  
  dsn <- deparse(substitute(data))
  dfn <- names(data)
  
  if (length(from) != length(to)) {
    abort_internal("from and to not same length")
  }
  
  if (length(dfn) < length(to)) {
    abort_internal("too many new names")
  }
  
  chng <- match(from, dfn)
  
  frm.in <- from %in% dfn
  if (!all(frm.in)) {
    abort_internal("Some of the 'from' names were not found in {.arg {dsn}}.")
  }
  
  if (length(to) != length(unique(to))) {
    abort_internal("New names not unique")
  }
  
  dfn.new <- dfn
  dfn.new[chng] <- to
  if (info) {
    inform_general("Changing names in {.arg {dsn}}.")
  }
  tmp <- rbind(from, to)
  dimnames(tmp)[[1]] <- c("From:", "To:")
  dimnames(tmp)[[2]] <- rep("", length(from))
  if (info) {
    print(tmp, quote = FALSE)
  }
  names(data) <- dfn.new
  data
}

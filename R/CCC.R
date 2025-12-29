#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: CCC.R                                                         #
# Contains: CCC function                                              #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 25/11/2023                                            #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal Function to Compute the Sampled Concordance
##'   Correlation Values.
##'
##' @description This function computes the sampled concordance correlation
##'   coefficient (CCC) between two numeric vectors. It is used internally
##'   for statistical analysis.
##'
##' @param Y1 A numeric vector.
##' @param Y2 A numeric vector, typically paired with Y1.
##'
##' @return The concordance correlation coefficient between Y1 and Y2.
##'
##' @importFrom stats var cov
##'
##' @keywords internal
##' @examples
##' # Example usage:
##' # CCC(c(1, 2, 3), c(3, 2, 1))
CCC <- function(Y1, Y2) {
  validate_numeric_no_na(Y1, "Y1")
  validate_numeric_no_na(Y2, "Y2")
  validate_equal_length(Y1, Y2, "Y1", "Y2")
  ok1 <- validate_non_degenerate_var(Y1, "Y1")
  ok2 <- validate_non_degenerate_var(Y2, "Y2")
  if (!ok1 || !ok2) {
    return(NA_real_)
  }

  s12 <- stats::cov(Y1, Y2)
  if (!is.finite(s12)) {
    warn_general("Covariance between {.arg Y1} and {.arg Y2} is not finite; returning {.val NA_real_}.")
    return(NA_real_)
  }

  m1  <- mean(Y1)
  m2  <- mean(Y2)
  v1  <- stats::var(Y1)
  v2  <- stats::var(Y2)

  2 * s12 / (v1 + v2 + (m1 - m2)^2)
}

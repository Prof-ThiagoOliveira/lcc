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
  m1  <- mean(Y1)
  m2  <- mean(Y2)
  v1  <- var(Y1)
  v2  <- var(Y2)
  s12 <- cov(Y1, Y2)
  
  2 * s12 / (v1 + v2 + (m1 - m2)^2)
}

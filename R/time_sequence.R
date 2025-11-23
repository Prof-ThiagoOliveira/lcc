#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: time_sequence.R                                               #
# Contains: time_lcc function                                         #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Regular Sequence Generator for Time Variable
##'
##' @description Generates a regular sequence for the time variable, including
##'   the unique values from the input time vector. This function is used
##'   internally to construct LCC, LPC, and LA curves and their simultaneous
##'   confidence intervals.
##'
##' @param time A numeric vector of unique time values.
##' @param from The starting (minimum) value for the time sequence.
##' @param to The ending (maximum) value for the time sequence.
##' @param n Desired length of the sequence (integer). Typically, a value
##'   between 30 and 50 is adequate.
##'
##' @return A numeric vector containing a regular sequence of time values,
##'   including the unique values from the input time vector.
##'
##' @keywords internal
##'
##' @examples
##' data(hue)
##' attach(hue)
##' time_lcc(time = Time, from = min(Time), to = max(Time), n = 30)
##' detach(hue)
##'
##' @export

time_lcc <- function(time, from, to, n) {
  if (from < min(time)) from <- min(time)
  if (to > max(time)) to <- max(time)
  
  seq_time <- seq.int(from, to, length.out = n)
  
  if (all(seq_time %in% time)) {
    tk.new <- time
  } else {
    tk.new <- unique(c(seq_time, time))
  }
  
  return(sort(tk.new))
}

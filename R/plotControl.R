#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: plotControl.R                                                 #
# Contains: plotControl function                                      #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 03/06/2020                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################
#' Control Settings for \code{lcc} Plots
#'
#' This function customizes the graphical control settings for plots of the 
#' \code{lcc} class. It allows for the adjustment of various aspects such as 
#' shape, color, and size of plot elements, as well as axis labels. The function 
#' returns a list containing all these settings, which can be used in plotting 
#' functions for \code{lcc} objects.
#'
#' @param plot Logical flag to include an initial plot. If set to \code{TRUE}, 
#'   a \code{\link[ggplot2]{ggplot}} object with an initial plot for \code{lcc} 
#'   class is returned. Defaults to \code{TRUE}.
#' @param shape Numeric value specifying the shape of points in the plot. 
#'   Acceptable values are from 0 to 25, and 32 to 127. See 
#'   \code{\link[ggplot2]{aes}} for details on setting shape. Default is \code{1}.
#' @param colour String specifying the color of lines in the plot. 
#'   Default color is \code{"black"}.
#' @param size Numeric value specifying the size of lines in the plot, given in 
#'   millimeters. See \code{\link[ggplot2]{aes}} for details on setting size. 
#'   Default is \code{0.5}.
#' @param xlab Title for the x-axis, defaulting to \code{"Time"}.
#' @param ylab Title for the y-axis, defaulting to \code{"LCC"}.
#'
#' @return A list with the specified graphical parameters.
#'
#' @author Thiago de Paula Oliveira,
#'   \email{thiago.paula.oliveira@@alumni.usp.br}
#'
#' @importFrom ggplot2 ggplot aes
#'
#' @keywords internal
plotControl <- function(plot = TRUE, shape = 1, colour = "black", size = 0.5,
                        xlab = "Time", ylab = "LCC") {
  list(plot = plot, shape = shape, colour = colour, size = size,
       xlab = xlab, ylab = ylab)
}

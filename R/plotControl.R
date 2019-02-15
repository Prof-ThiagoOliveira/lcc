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
# Last update: 18/06/2018                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################
##' @title Specifying graphical control values for \code{lcc} class
##'
##' @description The values supplied in the \code{plotControl()} call
##'   replace the defaults, and a \code{\link{list}} with all settings
##'   is returned.
##'
##' @return a list with components for each of the possible arguments.
##'
##' @param plot an optional to include an initial plot. If \code{TRUE},
##'   the default, returns a \code{\link[ggplot2]{ggplot}} object with a
##'   initial plot for \code{lcc} class. If \code{FALSE} never includes.
##'
##' @param shape Draw points considering a shape parameter. Legal shape
##'   values are the numbers 0 to 25, and 32 to 127; see
##'   \code{\link[ggplot2]{aes_linetype_size_shape}}. Default is
##'   \code{1}.
##'
##' @param colour an specification for lines color. Default is
##'   \code{"black"}.
##'
##' @param size an specification for lines size. Should be specified
##'   with a numerical value (in millimetres); see
##'   \code{\link[ggplot2]{aes_linetype_size_shape}}. Default is
##'   \code{0.5}.
##'
##' @param xlab a title for the \code{x} axis.  Default is
##'   \code{"Time"}.
##'
##' @param LCC_ylab title for the \code{y} axis related with
##'   longitudinal concordance correlation (LCC). Default is
##'   \code{"LCC"}.
##'
##' @param LPC_ylab title for the \code{y} axis related with
##'   longitudinal Pearson correlation (LPC). Default is \code{"LPC"}.
##'
##' @param LA_ylab title for the \code{y} axis related with longitudinal
##'   accuracy (LA). Default is \code{"LA"}.
##'
##' @param LCC_scale_y_continuous a numeric vector of length two
##'   providing limits of the scale related with longitudinal
##'   concordance correlation (LCC). Default is \code{c(0,1)}.
##'
##' @param LPC_scale_y_continuous a numeric vector of length two
##'   providing limits of the scale related with longitudinal
##'   concordance correlation (LCC). Default is \code{c(0,1)}.
##'
##' @param LA_scale_y_continuous a numeric vector of length two
##'   providing limits of the scale related with longitudinal
##'   concordance correlation (LCC). Default is \code{c(0,1)}.
##'
##' @param all.plot viewports functions for \code{lcc} class. If
##'   \code{TRUE}, the default, returns a object created by
##'   \code{\link[grid]{viewport}} function with multiple plots in a
##'   single page. If \code{FALSE} returns a single
##'   \code{\link[ggplot2]{ggplot}} object by page.
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@usp.br}
##'
##' @import ggplot2
##'
##' @importFrom grid viewport
##'
##' @keywords internal
plotControl<-function(plot= TRUE, shape=1, colour="black", size=0.5,
                       xlab = "Time", LCC_ylab = "LCC",
                       LPC_ylab = "LPC", LA_ylab = "LA",
                       LCC_scale_y_continuous=c(0,1),
                       LPC_scale_y_continuous=c(0,1),
                       LA_scale_y_continuous=c(0,1),
                       all.plot=TRUE)
{
       list(plot = plot, shape=shape, colour=colour, size=size,
            xlab = xlab, LCC_ylab = LCC_ylab, LPC_ylab = LPC_ylab,
            LA_ylab = LA_ylab,
            LCC_scale_y_continuous = LCC_scale_y_continuous,
            LPC_scale_y_continuous = LPC_scale_y_continuous,
            LA_scale_y_continuous = LA_scale_y_continuous,
            all.plot = all.plot)
}

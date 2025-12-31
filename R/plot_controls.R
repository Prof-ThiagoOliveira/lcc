#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: plot_controls.R                                               #
# Contains: plotControl helpers                                       #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 23/11/2025                                             #
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
#' @param ci_fill Fill colour used for confidence interval ribbons. If
#'   \code{NULL}, defaults to the main line colour.
#' @param ci_alpha Alpha level for the confidence interval ribbons. If
#'   \code{NULL}, a sensible default is applied.
#' @param point_alpha Alpha level for points (e.g., CCC markers). If
#'   \code{NULL}, a sensible default is applied.
#' @param xlab Title for the x-axis, defaulting to \code{"Time"}.
#' @param ylab Title for the y-axis, defaulting to \code{"LCC"}.
#' @param scale_y_continuous Optional list of arguments passed to
#'   \code{ggplot2::scale_y_continuous()} for full control of the y-axis.
#' @param expand_y Optional expansion for the y-axis (passed to the
#'   \code{expand} argument of \code{scale_y_continuous}). If \code{NULL},
#'   a light default expansion is applied for readability.
#'
#' @return A list with the specified graphical parameters.
#'
#' @author Thiago de Paula Oliveira,
#'   \email{thiago.paula.oliveira@@alumni.usp.br}
#'
#' @importFrom ggplot2 ggplot aes
#'
#' @examples
#' ctrl <- plotControl(
#'   plot        = FALSE,
#'   colour      = "steelblue",
#'   point_alpha = 0.4,
#'   ci_alpha    = 0.2,
#'   xlab        = "Measurement time",
#'   ylab        = "Concordance"
#' )
#' str(ctrl)
#'
#' @export
plotControl <- function(plot        = TRUE,
                        shape       = 16,
                        colour      = "#1B4F72",
                        size        = 0.7,
                        ci_fill     = NULL,
                        ci_alpha    = NULL,
                        point_alpha = NULL,
                        xlab        = "Time",
                        ylab        = "LCC",
                        scale_y_continuous = NULL,
                        expand_y    = NULL) {
  list(
    plot        = plot,
    shape       = shape,
    colour      = colour,
    size        = size,
    ci_fill     = ci_fill,
    ci_alpha    = ci_alpha,
    point_alpha = point_alpha,
    xlab        = xlab,
    ylab        = ylab,
    scale_y_continuous = scale_y_continuous,
    expand_y    = expand_y
  )
}


#' Internal base theme used by all ggplot-based summaries
#' @keywords Internal
.lcc_default_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position    = "bottom",
      plot.title         = ggplot2::element_text(hjust = 0.5, face = "bold"),
      strip.background   = ggplot2::element_blank(),
      strip.text         = ggplot2::element_text(face = "bold"),
      plot.margin        = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
    )
}


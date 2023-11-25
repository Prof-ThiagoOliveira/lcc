#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: plotmethod_lcc.R                                              #
# Contains: lccPlot function                                          #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################
#' Plot Fitted Curves from an \code{lcc} Object
#'
#' This function generates a plot of predictions versus the time covariate for 
#' an \code{lcc} object. Predicted values are connected by lines, while actual 
#' observations are denoted by circles. If \code{components=TRUE} was used in the 
#' \code{lcc} object, individual plots for each statistic (LCC, LPC, and LA) are 
#' produced on separate pages.
#'
#' @param obj An object inheriting from class "lcc", representing a fitted lcc model.
#' @param type Character string specifying the type of plot to generate. 
#'   \itemize{
#'     \item \code{"lcc"}: Produces the LCC plot.
#'     \item \code{"lpc"}: Produces the LPC plot. Available only if \code{components = TRUE}.
#'     \item \code{"la"}: Produces the LA plot. Available only if \code{components = TRUE}.
#'   }
#' @param control A list of graphical control values or character strings returned 
#'   by the \code{\link{plotControl}} function. Defaults to an empty list. 
#'   The list can contain components like \code{shape}, \code{colour}, \code{size},
#'   \code{xlab}, \code{ylab}, \code{scale_y_continuous}, and \code{all.plot}.
#' @param ... Additional arguments passed to the 
#'   \code{\link[ggplot2]{facet_wrap}} function.
#'
#' @return An object of class \code{ggplot} or \code{viewport}, depending on the 
#'   \code{all.plot} setting in \code{control}.
#'
#' @examples
#' data(hue)
#' # Second degree polynomial model with random intercept, slope and quadratic term
#' fm1 <- lcc(data = hue, subject = "Fruit", resp = "H_mean",
#'           method = "Method", time = "Time", qf = 2, qr = 2, components = TRUE)
#' lccPlot(fm1, type = "lcc")
#' lccPlot(fm1, type = "lpc")
#' lccPlot(fm1, type = "la")
#'
#' # Using ggplot2 themes
#' lccPlot(fm1, type = "lpc") + theme_bw() + labs(x = "Time (Days)", y = "LPC Value")
#'
#' # Generating and saving plots
#' \dontrun{
#'   ggsave("lccPlot.pdf", lccPlot(fm1, type = "lcc"))
#' }
#'
#' @seealso \code{\link[lcc]{lcc}}, \code{\link{plotControl}}
#' @importFrom ggplot2 ggplot facet_wrap
#' @importFrom grid viewport
#' @author Thiago de Paula Oliveira,
#'   \email{thiago.paula.oliveira@@alumni.usp.br}
#' @export
lccPlot<-function(obj, type = "lcc", control = list(), ...){
  if (!inherits(obj, "lcc")) stop("Object must inherit from class \"lcc\"",
                                  call.=FALSE)
  # Arguments for the plot
  plot.cons<-plotControl(shape=1, colour="black",
                         size=0.5, xlab = "Time",
                         ylab = "LCC")
  if (type == "lpc") plot.cons$ylab = "LPC"
  if (type == "la") plot.cons$ylab = "LA"
  if(length(control)){
    nms <- names(control)
    if (!is.list(control) || is.null(nms))
      stop("'control' argument must be a named list")
    pos <- pmatch(nms, names(plot.cons))
    if (any(nap <- is.na(pos))) {
      warning(sprintf(ngettext(length(nap), "unrecognized plot element named %s ignored",
                               "unrecognized plot elements named %s ignored"),
                      paste(sQuote(nms[nap]), collapse = ", ")), domain = NA)
      pos <- pos[!nap]
      control <- control[!nap]
    }
    for(i in seq_len(length(pos))){
      plot.cons[[pos[i]]]<-control[[i]]
    }
  }
  #---------------------------------------------------------------------
  #Standard arguments
  #---------------------------------------------------------------------
  nd<-obj$plot_info$nd
  model<-obj$model
  tk.plot<-obj$plot_info$tk.plot
  tk.plot2<-obj$plot_info$tk.plot2
  ldb<-obj$plot_info$ldb
  ci<-obj$plot_info$ci
  components<-obj$plot_info$components
  if (components == FALSE &  type != "lcc") {
    stop("'lpc' and 'la' plots are only available if 'components = TRUE' in the 'lcc' call",
         call.= FALSE)
  }
  #---------------------------------------------------------------------
  if(ci==FALSE) {
    if(ldb == 1) {
      if (type == "lcc") {
        lccplot <- plot_lcc(rho=obj$plot_info$rho, tk.plot= tk.plot,
                            tk.plot2=tk.plot2, ldb=ldb,
                            model=model, ci = ci,
                            arg=plot.cons, ...)
      }
      if(components==TRUE){
        if (type == "lpc") {
          lccplot <- plot_lpc(LPC=obj$plot_info$rho.pearson,
                              tk.plot= tk.plot,
                              tk.plot2=tk.plot2, ldb=ldb,
                              model=model, ci = ci, arg = plot.cons,
                              ...)
        }
        if (type == "la") {
          lccplot <- plot_la(Cb=obj$plot_info$Cb, tk.plot= tk.plot,
                             tk.plot2=tk.plot2, ldb=ldb,
                             model=model, ci = ci, arg = plot.cons,
                             ...)
        }
      }
    } else {
      if (type == "lcc") {
        lccplot <- plot_lcc(rho=obj$plot_info$rho, tk.plot= tk.plot,
                            tk.plot2=tk.plot2, ldb=ldb, model=model,
                            ci = ci, arg = plot.cons, ...)
      }
      if(components==TRUE){
        if (type == "lpc") {
          lccplot <- plot_lpc(LPC=obj$plot_info$rho.pearson,
                              tk.plot= tk.plot,
                              tk.plot2=tk.plot2, ldb=ldb, model=model,
                              ci = ci, arg = plot.cons, ...)
        }
        if (type == "la") {
          lccplot <- plot_la(Cb=obj$plot_info$Cb, tk.plot= tk.plot,
                             tk.plot2=tk.plot2, ldb=ldb, model=model,
                             ci = ci, arg = plot.cons, ...)
        }
      }
    }
  }else{
    ENV.LCC<-obj$plot_info$ENV.LCC
    if (type == "lcc") {
      lccplot <- plot_lcc(rho=obj$plot_info$rho, ENV.LCC=ENV.LCC,
                          tk.plot= tk.plot, tk.plot2=tk.plot2, ldb=ldb,
                          model=model, ci = ci, arg = plot.cons, ...)
    }
    if(components==TRUE){
      if (type == "lpc") {
        ENV.LPC<-obj$plot_info$ENV.LPC
        lccplot <- plot_lpc(LPC=obj$plot_info$rho.pearson,ENV.LPC=ENV.LPC,
                            tk.plot= tk.plot, tk.plot2=tk.plot2, ldb=ldb,
                            model=model, ci = ci, arg = plot.cons, ...)
      }
      if (type == "la") {
        ENV.Cb<-obj$plot_info$ENV.LA
        lccplot <- plot_la(Cb=obj$plot_info$Cb,ENV.Cb = ENV.Cb, tk.plot= tk.plot,
                           tk.plot2=tk.plot2, ldb=ldb, model=model, ci = ci,
                           arg = plot.cons, ...)
      }
    }
  }
  return(invisible(lccplot))
}
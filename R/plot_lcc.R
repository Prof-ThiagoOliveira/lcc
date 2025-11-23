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
lccPlot <- function(obj, type = "lcc", control = list(), .) {
  if (!inherits(obj, "lcc"))
    stop("Object must inherit from class \"lcc\"", call. = FALSE)
  
  ## Base control defaults
  plot.cons <- plotControl(
    shape  = 1,
    colour = "black",
    size   = 0.5,
    xlab   = "Time",
    ylab   = "LCC"
  )
  if (type == "lpc") plot.cons$ylab <- "LPC"
  if (type == "la")  plot.cons$ylab <- "LA"
  
  ## User overrides
  if (length(control)) {
    nms <- names(control)
    if (!is.list(control) || is.null(nms)) {
      stop("'control' argument must be a named list", call. = FALSE)
    }
    pos <- pmatch(nms, names(plot.cons))
    if (any(nap <- is.na(pos))) {
      warning(sprintf(
        ngettext(
          length(nap),
          "unrecognized plot element named %s ignored",
          "unrecognized plot elements named %s ignored"
        ),
        paste(sQuote(nms[nap]), collapse = ", ")
      ), domain = NA)
      pos     <- pos[!nap]
      control <- control[!nap]
    }
    for (i in seq_along(pos)) {
      plot.cons[[pos[i]]] <- control[[i]]
    }
  }
  
  ## Standard arguments from fitted lcc object
  nd         <- obj$plot_info$nd
  model      <- obj$model
  tk.plot    <- obj$plot_info$tk.plot
  tk.plot2   <- obj$plot_info$tk.plot2
  ldb        <- obj$plot_info$ldb
  ci         <- obj$plot_info$ci
  components <- obj$plot_info$components
  
  if (!components && type != "lcc") {
    stop(
      "'lpc' and 'la' plots are only available if 'components = TRUE' ",
      "in the 'lcc' call",
      call. = FALSE
    )
  }
  
  ## Helper to call the correct internal plot_* wrapper
  call_plotter <- function(type, ci) {
    if (type == "lcc") {
      if (ci) {
        ENV.LCC <- obj$plot_info$ENV.LCC
        plot_lcc(
          rho     = obj$plot_info$rho,
          ENV.LCC = ENV.LCC,
          tk.plot = tk.plot,
          tk.plot2 = tk.plot2,
          ldb     = ldb,
          model   = model,
          ci      = TRUE,
          arg     = plot.cons,
          .
        )
      } else {
        plot_lcc(
          rho     = obj$plot_info$rho,
          ENV.LCC = NULL,
          tk.plot = tk.plot,
          tk.plot2 = tk.plot2,
          ldb     = ldb,
          model   = model,
          ci      = FALSE,
          arg     = plot.cons,
          .
        )
      }
    } else if (type == "lpc") {
      if (ci) {
        ENV.LPC <- obj$plot_info$ENV.LPC
        plot_lpc(
          LPC     = obj$plot_info$rho.pearson,
          ENV.LPC = ENV.LPC,
          tk.plot = tk.plot,
          tk.plot2 = tk.plot2,
          ldb     = ldb,
          model   = model,
          ci      = TRUE,
          arg     = plot.cons,
          .
        )
      } else {
        plot_lpc(
          LPC     = obj$plot_info$rho.pearson,
          ENV.LPC = NULL,
          tk.plot = tk.plot,
          tk.plot2 = tk.plot2,
          ldb     = ldb,
          model   = model,
          ci      = FALSE,
          arg     = plot.cons,
          .
        )
      }
    } else if (type == "la") {
      if (ci) {
        ENV.Cb <- obj$plot_info$ENV.LA
        plot_la(
          Cb      = obj$plot_info$Cb,
          ENV.Cb  = ENV.Cb,
          tk.plot = tk.plot,
          tk.plot2 = tk.plot2,
          ldb     = ldb,
          model   = model,
          ci      = TRUE,
          arg     = plot.cons,
          .
        )
      } else {
        plot_la(
          Cb      = obj$plot_info$Cb,
          ENV.Cb  = NULL,
          tk.plot = tk.plot,
          tk.plot2 = tk.plot2,
          ldb     = ldb,
          model   = model,
          ci      = FALSE,
          arg     = plot.cons,
          .
        )
      }
    } else {
      stop("Unknown 'type' in lccPlot: ", type, call. = FALSE)
    }
  }
  
  ## ci / no-ci handled only through ENV.* arguments
  lccplot <- call_plotter(type = type, ci = ci)
  
  invisible(lccplot)
}


##' @keywords internal
plot_lcc <- function(rho, ENV.LCC, tk.plot, tk.plot2, ldb, model, ci, arg, .) {
  CCC <- CCC_lin(
    dataset = model$data,
    resp    = "resp",
    subject = "subject",
    method  = "method",
    time    = "time"
  )
  
  plotBuilder_lcc(
    rho     = rho,
    ENV.LCC = ENV.LCC,
    tk.plot = tk.plot,
    CCC     = CCC,
    tk.plot2 = tk.plot2,
    ldb     = ldb,
    model   = model,
    ci      = ci,
    arg     = arg,
    .
  )
}
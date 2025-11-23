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
lccPlot <- function(obj, type = "lcc", control = list(), ...) {
  if (!inherits(obj, "lcc"))
    stop("Object must inherit from class \"lcc\"", call. = FALSE)
  
  ## Base control defaults (modernised)
  plot.cons <- plotControl(
    shape       = 16,        # solid points
    colour      = "#1B4F72", # muted blue
    size        = 0.7,       # line width (mm)
    ci_fill     = NULL,      # defaults set just below
    ci_alpha    = NULL,
    point_alpha = NULL,
    xlab        = "Time",
    ylab        = "LCC"
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
  
  ## Fill in modern defaults for new elements if user did not specify
  if (is.null(plot.cons$ci_fill))     plot.cons$ci_fill     <- plot.cons$colour
  if (is.null(plot.cons$ci_alpha))    plot.cons$ci_alpha    <- 0.15
  if (is.null(plot.cons$point_alpha)) plot.cons$point_alpha <- 0.8
  
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
      plot_lcc(
        rho     = obj$plot_info$rho,
        ENV.LCC = if (ci) obj$plot_info$ENV.LCC else NULL,
        tk.plot = tk.plot,
        tk.plot2 = tk.plot2,
        ldb     = ldb,
        model   = model,
        ci      = ci,
        arg     = plot.cons,
        ...
      )
    } else if (type == "lpc") {
      plot_lpc(
        LPC     = obj$plot_info$rho.pearson,
        ENV.LPC = if (ci) obj$plot_info$ENV.LPC else NULL,
        tk.plot = tk.plot,
        tk.plot2 = tk.plot2,
        ldb     = ldb,
        model   = model,
        ci      = ci,
        arg     = plot.cons,
        ...
      )
    } else if (type == "la") {
      plot_la(
        Cb      = obj$plot_info$Cb,
        ENV.Cb  = if (ci) obj$plot_info$ENV.LA else NULL,
        tk.plot = tk.plot,
        tk.plot2 = tk.plot2,
        ldb     = ldb,
        model   = model,
        ci      = ci,
        arg     = plot.cons,
        ...
      )
    } else {
      stop("Unknown 'type' in lccPlot: ", type, call. = FALSE)
    }
  }
  
  lccplot <- call_plotter(type = type, ci = ci)
  invisible(lccplot)
}




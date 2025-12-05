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
  if (!inherits(obj, "lcc")) {
    abort_input("Object must inherit from class \"lcc\"")
  }
  
  ## Base defaults
  plot.cons <- plotControl(
    shape       = 16,
    colour      = "#1B4F72",
    size        = 0.7,
    ci_fill     = NULL,
    ci_alpha    = NULL,
    point_alpha = NULL,
    xlab        = "Time",
    ylab        = "LCC"
  )
  if (type == "lpc") plot.cons$ylab <- "LPC"
  if (type == "la")  plot.cons$ylab <- "LA"
  
  ## Apply user overrides from 'control'
  if (length(control)) {
    nms <- names(control)
    if (!is.list(control) || is.null(nms)) {
      abort_input("'control' argument must be a named list")
    }
    pos <- pmatch(nms, names(plot.cons))
    if (any(nap <- is.na(pos))) {
      warn_general(sprintf(
        ngettext(
          length(nap),
          "unrecognized plot element named %s ignored",
          "unrecognized plot elements named %s ignored"
        ),
        paste(sQuote(nms[nap]), collapse = ", ")
      ))
      pos     <- pos[!nap]
      control <- control[!nap]
    }
    for (i in seq_along(pos)) {
      plot.cons[[pos[i]]] <- control[[i]]
    }
  }
  
  ## Fill in defaults for CI aesthetics if still NULL
  if (is.null(plot.cons$ci_fill))     plot.cons$ci_fill     <- plot.cons$colour
  if (is.null(plot.cons$ci_alpha))    plot.cons$ci_alpha    <- 0.15
  if (is.null(plot.cons$point_alpha)) plot.cons$point_alpha <- 0.8
  
  ## Extract plotting info from fitted lcc object
  model      <- obj$model
  tk.plot    <- obj$plot_info$tk.plot
  tk.plot2   <- obj$plot_info$tk.plot2
  ldb        <- obj$plot_info$ldb
  ci         <- obj$plot_info$ci
  components <- obj$plot_info$components
  
  if (!components && type != "lcc") {
    abort_input("'lpc' and 'la' plots are only available if 'components = TRUE' in the 'lcc' call")
  }
  
  # Precompute CCC/Pearson once per call to avoid re-splitting data
  CCC_vals <- CCC_lin(
    dataset = model$data,
    resp    = "resp",
    subject = "subject",
    method  = "method",
    time    = "time"
  )
  Pearson_vals <- NULL
  if (type != "lcc") {
    Pearson_vals <- Pearson(
      dataset = model$data,
      resp    = "resp",
      subject = "subject",
      method  = "method",
      time    = "time"
    )
  }
  
  ## Dispatch to internal plotters, all of which RETURN a ggplot object
  res <- switch(
    type,
    "lcc" = plot_lcc(
      rho      = obj$plot_info$rho,
      ENV.LCC  = if (ci) obj$plot_info$ENV.LCC else NULL,
      tk.plot  = tk.plot,
      tk.plot2 = tk.plot2,
      ldb      = ldb,
      model    = model,
      ci       = ci,
      arg      = plot.cons,
      CCC_vals = CCC_vals,
      ...
    ),
    "lpc" = plot_lpc(
      LPC      = obj$plot_info$rho.pearson,
      ENV.LPC  = if (ci) obj$plot_info$ENV.LPC else NULL,
      tk.plot  = tk.plot,
      tk.plot2 = tk.plot2,
      ldb      = ldb,
      model    = model,
      ci       = ci,
      arg      = plot.cons,
      Pearson_vals = Pearson_vals,
      ...
    ),
    "la"  = plot_la(
      Cb       = obj$plot_info$Cb,
      ENV.Cb   = if (ci) obj$plot_info$ENV.LA else NULL,
      tk.plot  = tk.plot,
      tk.plot2 = tk.plot2,
      ldb      = ldb,
      model    = model,
      ci       = ci,
      arg      = plot.cons,
      CCC_vals = CCC_vals,
      Pearson_vals = Pearson_vals,
      ...
    ),
    abort_input("Unknown 'type' in lccPlot: {type}")
  )
  
  ## Optionally print for interactive use
  if (isTRUE(plot.cons$plot)) {
    print(res)
  }
  
  invisible(res)
}

#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: plotBuilder_lcc.R                                             #
# Contains: plotBuilder_lcc function                                  #
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
plotControl <- function(plot        = TRUE,
                        shape       = 16,
                        colour      = "#1B4F72",
                        size        = 0.7,
                        ci_fill     = NULL,
                        ci_alpha    = NULL,
                        point_alpha = NULL,
                        xlab        = "Time",
                        ylab        = "LCC") {
  list(
    plot        = plot,
    shape       = shape,
    colour      = colour,
    size        = size,
    ci_fill     = ci_fill,
    ci_alpha    = ci_alpha,
    point_alpha = point_alpha,
    xlab        = xlab,
    ylab        = ylab
  )
}


#' Internal base theme used by all ggplot-based summaries
#' @keywords Internal
.lcc_default_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position    = "bottom",
      plot.title         = ggplot2::element_text(hjust = 0.5, face = "bold"),
      strip.background   = ggplot2::element_blank(),
      strip.text         = ggplot2::element_text(face = "bold")
    )
}


##' @title Internal Function to Estimate the Sampled Concordance
##'   Correlation Coefficient.
##'
##' @description This function is internally called to estimate
##'   the sampled concordance correlation coefficient.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @importFrom stats cor cov
##'
##' @keywords internal
CCC_lin <- function(dataset, resp, subject, method, time) {
  # We assume here that `dataset` has already been prepared by dataBuilder
  selectedData <- subset(dataset, select = c(resp, method, time, subject))
  dataByMethod <- split(selectedData, selectedData[[method]])
  methodLevels <- levels(selectedData[[method]])
  
  calculateCCC_fast <- function(Y1, Y2, time) {
    n          <- length(time)
    time_fac   <- as.factor(time)
    idx_by_t   <- split(seq_len(n), time_fac)
    
    Y1_full    <- rep(Y1, length.out = n)
    Y2_full    <- rep(Y2, length.out = n)
    
    ccc_vec <- vapply(
      idx_by_t,
      function(idx) {
        y1 <- Y1_full[idx]
        y2 <- Y2_full[idx]
        m1 <- mean(y1)
        m2 <- mean(y2)
        s1 <- var(y1)
        s2 <- var(y2)
        s12 <- cov(y1, y2)
        2 * s12 / (s1 + s2 + (m1 - m2)^2)
      },
      numeric(1L)
    )
    data.frame(V1 = unname(ccc_vec))
  }
  
  CCC.Lin <- lapply(
    seq(2L, length(methodLevels)),
    function(i) {
      calculateCCC_fast(
        Y1   = dataByMethod[[1L]][[resp]],
        Y2   = dataByMethod[[i]][[resp]],
        time = selectedData[[time]]
      )
    }
  )
  
  CCC.Lin
}


#' @title Estimate Sampled Pearson Correlation
#'
#' @description Internally called function to estimate the sampled Pearson correlation.
#'
#' @usage NULL
#'
#' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
#'
#' @importFrom stats cor
#'
#' @keywords internal
Pearson <- function(dataset, resp, subject, method, time) {
  selectedData <- subset(dataset, select = c(resp, method, time, subject))
  dataByMethod <- split(selectedData, selectedData[[method]])
  methodLevels <- levels(selectedData[[method]])
  
  calculateCorrelation_fast <- function(Y1, Y2, time) {
    n        <- length(time)
    time_fac <- as.factor(time)
    idx_by_t <- split(seq_len(n), time_fac)
    
    Y1_full  <- rep(Y1, length.out = n)
    Y2_full  <- rep(Y2, length.out = n)
    
    cor_vec <- vapply(
      idx_by_t,
      function(idx) cor(Y1_full[idx], Y2_full[idx]),
      numeric(1L)
    )
    
    data.frame(V1 = unname(cor_vec))
  }
  
  pearsonResults <- vector("list", length(methodLevels) - 1L)
  for (i in seq(2L, length(methodLevels))) {
    pearsonResults[[i - 1L]] <- calculateCorrelation_fast(
      Y1   = dataByMethod[[1L]][[resp]],
      Y2   = dataByMethod[[i]][[resp]],
      time = selectedData[[time]]
    )
  }
  
  pearsonResults
}

#######################################################################
# Internal wrappers: plot_lcc / plot_lpc / plot_la
#######################################################################

##' @keywords internal
plot_lcc <- function(rho, ENV.LCC, tk.plot, tk.plot2, ldb, model, ci, arg, ...) {
  CCC <- CCC_lin(
    dataset = model$data,
    resp    = "resp",
    subject = "subject",
    method  = "method",
    time    = "time"
  )
  
  plotBuilder_lcc(
    rho      = rho,
    ENV.LCC  = ENV.LCC,
    tk.plot  = tk.plot,
    CCC      = CCC,
    tk.plot2 = tk.plot2,
    ldb      = ldb,
    model    = model,
    ci       = ci,
    arg      = arg,
    ...
  )
}

#' @title Prepare Plot for LPC Function
#'
#' @description Internally called function to prepare data for the `plotBuilder_lpc` function.
#'
#' @usage NULL
#'
#' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
#'
#' @keywords internal
plot_lpc <- function(LPC, ENV.LPC, tk.plot, tk.plot2, ldb, model, ci, arg, ...) {
  Pearson <- Pearson(
    dataset = model$data,
    resp    = "resp",
    subject = "subject",
    method  = "method",
    time    = "time"
  )
  
  plotBuilder_lpc(
    LPC     = LPC,
    ENV.LPC = ENV.LPC,
    tk.plot = tk.plot,
    Pearson = Pearson,
    tk.plot2 = tk.plot2,
    ldb     = ldb,
    model   = model,
    ci      = ci,
    arg     = arg,
    ...
  )
}

##' @title Internal Function to Prepare the \code{plotBuilder_la} Function
##'
##' @description This function is internally called to prepare
##'   the \code{\link[lcc]{plotBuilder_la}} function.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @keywords internal
plot_la <- function(Cb, ENV.Cb, tk.plot, tk.plot2, ldb, model, ci, arg, ...) {
  CCC <- CCC_lin(
    dataset = model$data,
    resp    = "resp",
    subject = "subject",
    method  = "method",
    time    = "time"
  )
  Pearson <- Pearson(
    dataset = model$data,
    resp    = "resp",
    subject = "subject",
    method  = "method",
    time    = "time"
  )
  
  plotBuilder_la(
    CCC     = CCC,
    Pearson = Pearson,
    Cb      = Cb,
    ENV.Cb  = ENV.Cb,
    tk.plot = tk.plot,
    tk.plot2 = tk.plot2,
    ldb     = ldb,
    model   = model,
    ci      = ci,
    arg     = arg,
    ...
  )
}

#######################################################################
# Low-level builders: LCC / LPC / LA
#######################################################################

#' Generate a Longitudinal Concordance Correlation Plot
#'
#' Produces a longitudinal concordance correlation plot from fitted 
#' and sampled values with optional non-parametric confidence intervals.
#'
#' @param rho Vector of LCC values.
#' @param ENV.LCC Environment matrix for LCC values, used for confidence intervals.
#' @param tk.plot Time points for LCC values.
#' @param CCC Matrix or list of CCC values.
#' @param tk.plot2 Time points for CCC values.
#' @param ldb Number of levels in the data.
#' @param model The model object from which data was extracted.
#' @param ci Logical, indicating if confidence intervals should be included.
#' @param arg List of graphical arguments.
#' @param ... Additional arguments for ggplot.
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon labs theme 
#' element_text geom_path
#' @keywords internal
plotBuilder_lcc <- function(rho, ENV.LCC, tk.plot, CCC,
                            tk.plot2, ldb, model, ci, arg, ...) {
  
  method_levels <- levels(model$data$method)
  level_label <- function(i) {
    paste(method_levels[i + 1L], "vs.", method_levels[1L])
  }
  
  ci_fill     <- if (!is.null(arg$ci_fill)) arg$ci_fill else arg$colour
  ci_alpha    <- if (!is.null(arg$ci_alpha)) arg$ci_alpha else 0.15
  point_alpha <- if (!is.null(arg$point_alpha)) arg$point_alpha else 0.8
  
  if (!ci) {
    ## No CI
    if (ldb == 1L) {
      data_plot  <- data.frame(LCC = rho,          Time = tk.plot)
      data_plot2 <- data.frame(CCC = CCC[[1L]]$V1, Time = tk.plot2)
      
      Plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = Time, y = LCC)) +
        ggplot2::geom_path(colour = arg$colour, linewidth = arg$size) +
        ggplot2::geom_point(
          data  = data_plot2,
          ggplot2::aes(x = Time, y = CCC),
          shape = arg$shape,
          alpha = point_alpha
        ) +
        ggplot2::labs(
          x     = arg$xlab,
          y     = arg$ylab,
          title = level_label(1L)
        ) +
        .lcc_default_theme()
      
    } else {
      data_main <- vector("list", ldb)
      data_ccc  <- vector("list", ldb)
      for (i in seq_len(ldb)) {
        lab_i <- level_label(i)
        data_main[[i]] <- data.frame(
          LCC   = rho[, i],
          Time  = tk.plot,
          Level = lab_i
        )
        data_ccc[[i]] <- data.frame(
          CCC   = CCC[[i]]$V1,
          Time  = tk.plot2,
          Level = lab_i
        )
      }
      data_plot  <- do.call(rbind.data.frame, data_main)
      data_plot2 <- do.call(rbind.data.frame, data_ccc)
      
      Plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = Time, y = LCC)) +
        ggplot2::geom_line(colour = arg$colour, linewidth = arg$size) +
        ggplot2::geom_point(
          data  = data_plot2,
          ggplot2::aes(x = Time, y = CCC),
          shape = arg$shape,
          alpha = point_alpha
        ) +
        ggplot2::facet_wrap(~Level) +
        ggplot2::labs(x = arg$xlab, y = arg$ylab) +
        .lcc_default_theme()
    }
    
  } else {
    ## With CI
    if (ldb == 1L) {
      data_plot <- data.frame(
        LCC       = rho,
        Time      = tk.plot,
        lower_rho = t(ENV.LCC)[, 1L],
        upper_rho = t(ENV.LCC)[, 2L]
      )
      data_plot2 <- data.frame(
        CCC  = CCC[[1L]]$V1,
        Time = tk.plot2
      )
      
      Plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = Time, y = LCC)) +
        ggplot2::geom_line(colour = arg$colour, linewidth = arg$size) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower_rho, ymax = upper_rho),
          fill        = ci_fill,
          alpha       = ci_alpha,
          colour      = NA,
          show.legend = FALSE
        ) +
        ggplot2::geom_point(
          data  = data_plot2,
          ggplot2::aes(x = Time, y = CCC),
          shape = arg$shape,
          alpha = point_alpha
        ) +
        ggplot2::labs(
          x     = arg$xlab,
          y     = arg$ylab,
          title = level_label(1L)
        ) +
        .lcc_default_theme()
      
    } else {
      data_main <- vector("list", ldb)
      data_ccc  <- vector("list", ldb)
      for (i in seq_len(ldb)) {
        lab_i <- level_label(i)
        env_i <- ENV.LCC[[i]]
        data_main[[i]] <- data.frame(
          LCC       = rho[, i],
          Time      = tk.plot,
          lower_rho = t(env_i)[, 1L],
          upper_rho = t(env_i)[, 2L],
          Level     = lab_i
        )
        data_ccc[[i]] <- data.frame(
          CCC   = CCC[[i]]$V1,
          Time  = tk.plot2,
          Level = lab_i
        )
      }
      data_plot  <- do.call(rbind.data.frame, data_main)
      data_plot2 <- do.call(rbind.data.frame, data_ccc)
      
      Plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = Time, y = LCC)) +
        ggplot2::geom_line(colour = arg$colour, linewidth = arg$size) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower_rho, ymax = upper_rho),
          fill        = ci_fill,
          alpha       = ci_alpha,
          colour      = NA,
          show.legend = FALSE
        ) +
        ggplot2::geom_point(
          data  = data_plot2,
          ggplot2::aes(x = Time, y = CCC),
          shape = arg$shape,
          alpha = point_alpha
        ) +
        ggplot2::facet_wrap(~Level) +
        ggplot2::labs(x = arg$xlab, y = arg$ylab) +
        .lcc_default_theme()
    }
  }
  
  Plot
}


##' @title Internal Function to Produce a Longitudinal Pearson Correlation Plot
##'
##' @description Produces a longitudinal Pearson correlation plot from fitted 
##'   and sampled values, with optional non-parametric confidence intervals.
##'
##' @details Returns an initial plot for the longitudinal Pearson correlation.
##'
##' @usage NULL
##'
##' @importFrom ggplot2 ggplot geom_line geom_point geom_ribbon labs theme element_text facet_wrap
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##' @keywords internal
plotBuilder_lpc <- function(LPC, ENV.LPC, tk.plot, Pearson,
                            tk.plot2, ldb, model, ci, arg, ...) {
  
  method_levels <- levels(model$data$method)
  level_label <- function(i) {
    paste(method_levels[i + 1L], "vs.", method_levels[1L])
  }
  
  ci_fill     <- if (!is.null(arg$ci_fill)) arg$ci_fill else arg$colour
  ci_alpha    <- if (!is.null(arg$ci_alpha)) arg$ci_alpha else 0.15
  point_alpha <- if (!is.null(arg$point_alpha)) arg$point_alpha else 0.8
  
  if (!ci) {
    if (ldb == 1L) {
      data_plot  <- data.frame(LPC = LPC,                 Time = tk.plot)
      data_plot2 <- data.frame(Pearson = Pearson[[1L]]$V1, Time = tk.plot2)
      
      Plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = Time, y = LPC)) +
        ggplot2::geom_line(colour = arg$colour, linewidth = arg$size) +
        ggplot2::geom_point(
          data  = data_plot2,
          ggplot2::aes(x = Time, y = Pearson),
          shape = arg$shape,
          alpha = point_alpha
        ) +
        ggplot2::labs(
          x     = arg$xlab,
          y     = arg$ylab,
          title = level_label(1L)
        ) +
        .lcc_default_theme()
      
    } else {
      data_main <- vector("list", ldb)
      data_pear <- vector("list", ldb)
      for (i in seq_len(ldb)) {
        lab_i <- level_label(i)
        data_main[[i]] <- data.frame(
          LPC   = LPC[, i],
          Time  = tk.plot,
          Level = lab_i
        )
        data_pear[[i]] <- data.frame(
          Pearson = Pearson[[i]]$V1,
          Time    = tk.plot2,
          Level   = lab_i
        )
      }
      data_plot  <- do.call(rbind.data.frame, data_main)
      data_plot2 <- do.call(rbind.data.frame, data_pear)
      
      Plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = Time, y = LPC)) +
        ggplot2::geom_line(colour = arg$colour, linewidth = arg$size) +
        ggplot2::geom_point(
          data  = data_plot2,
          ggplot2::aes(x = Time, y = Pearson),
          shape = arg$shape,
          alpha = point_alpha
        ) +
        ggplot2::facet_wrap(~Level) +
        ggplot2::labs(x = arg$xlab, y = arg$ylab) +
        .lcc_default_theme()
    }
    
  } else {
    if (ldb == 1L) {
      data_plot <- data.frame(
        LPC        = LPC,
        Time       = tk.plot,
        lower_LPC  = t(ENV.LPC)[, 1L],
        upper_LPC  = t(ENV.LPC)[, 2L]
      )
      data_plot2 <- data.frame(
        Pearson = Pearson[[1L]]$V1,
        Time    = tk.plot2
      )
      
      Plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = Time, y = LPC)) +
        ggplot2::geom_line(colour = arg$colour, linewidth = arg$size) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower_LPC, ymax = upper_LPC),
          fill        = ci_fill,
          alpha       = ci_alpha,
          colour      = NA,
          show.legend = FALSE
        ) +
        ggplot2::geom_point(
          data  = data_plot2,
          ggplot2::aes(x = Time, y = Pearson),
          shape = arg$shape,
          alpha = point_alpha
        ) +
        ggplot2::labs(
          x     = arg$xlab,
          y     = arg$ylab,
          title = level_label(1L)
        ) +
        .lcc_default_theme()
      
    } else {
      data_main <- vector("list", ldb)
      data_pear <- vector("list", ldb)
      for (i in seq_len(ldb)) {
        lab_i <- level_label(i)
        env_i <- ENV.LPC[[i]]
        data_main[[i]] <- data.frame(
          LPC       = LPC[, i],
          Time      = tk.plot,
          lower_LPC = t(env_i)[, 1L],
          upper_LPC = t(env_i)[, 2L],
          Level     = lab_i
        )
        data_pear[[i]] <- data.frame(
          Pearson = Pearson[[i]]$V1,
          Time    = tk.plot2,
          Level   = lab_i
        )
      }
      data_plot  <- do.call(rbind.data.frame, data_main)
      data_plot2 <- do.call(rbind.data.frame, data_pear)
      
      Plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = Time, y = LPC)) +
        ggplot2::geom_line(colour = arg$colour, linewidth = arg$size) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower_LPC, ymax = upper_LPC),
          fill        = ci_fill,
          alpha       = ci_alpha,
          colour      = NA,
          show.legend = FALSE
        ) +
        ggplot2::geom_point(
          data  = data_plot2,
          ggplot2::aes(x = Time, y = Pearson),
          shape = arg$shape,
          alpha = point_alpha
        ) +
        ggplot2::facet_wrap(~Level) +
        ggplot2::labs(x = arg$xlab, y = arg$ylab) +
        .lcc_default_theme()
    }
  }
  
  Plot
}


##' @title Internal Function to Produces a Longitudinal Accuracy Plot.
##'
##' @description Produces a longitudinal accuracy plot from fitted and sampled values 
##'   with optional non-parametric confidence intervals.
##'
##' @details Returns an initial plot for the longitudinal accuracy correlation.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @importFrom ggplot2 ggplot geom_line geom_point geom_ribbon labs theme element_text ggtitle
##' @keywords internal
plotBuilder_la <- function(CCC, Pearson, Cb, ENV.Cb,
                           tk.plot, tk.plot2, ldb, model, ci, arg, ...) {
  
  method_levels <- levels(model$data$method)
  level_label <- function(i) {
    paste(method_levels[i + 1L], "vs.", method_levels[1L])
  }
  
  LA_fun <- function(i) CCC[[i]]$V1 / Pearson[[i]]$V1
  
  ci_fill     <- if (!is.null(arg$ci_fill)) arg$ci_fill else arg$colour
  ci_alpha    <- if (!is.null(arg$ci_alpha)) arg$ci_alpha else 0.15
  point_alpha <- if (!is.null(arg$point_alpha)) arg$point_alpha else 0.8
  
  if (!ci) {
    if (ldb == 1L) {
      data_plot  <- data.frame(LA = Cb,         Time = tk.plot)
      data_plot2 <- data.frame(Cb = LA_fun(1L), Time = tk.plot2)
      
      Plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = Time, y = LA)) +
        ggplot2::geom_line(colour = arg$colour, linewidth = arg$size) +
        ggplot2::geom_point(
          data  = data_plot2,
          ggplot2::aes(x = Time, y = Cb),
          shape = arg$shape,
          alpha = point_alpha
        ) +
        ggplot2::labs(
          x     = arg$xlab,
          y     = arg$ylab,
          title = level_label(1L)
        ) +
        .lcc_default_theme()
      
    } else {
      data_main <- vector("list", ldb)
      data_la   <- vector("list", ldb)
      for (i in seq_len(ldb)) {
        lab_i <- level_label(i)
        data_main[[i]] <- data.frame(
          LA    = Cb[, i],
          Time  = tk.plot,
          Level = lab_i
        )
        data_la[[i]] <- data.frame(
          Cb    = LA_fun(i),
          Time  = tk.plot2,
          Level = lab_i
        )
      }
      data_plot  <- do.call(rbind.data.frame, data_main)
      data_plot2 <- do.call(rbind.data.frame, data_la)
      
      Plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = Time, y = LA)) +
        ggplot2::geom_line(colour = arg$colour, linewidth = arg$size) +
        ggplot2::geom_point(
          data  = data_plot2,
          ggplot2::aes(x = Time, y = Cb),
          shape = arg$shape,
          alpha = point_alpha
        ) +
        ggplot2::facet_wrap(~Level) +
        ggplot2::labs(x = arg$xlab, y = arg$ylab) +
        .lcc_default_theme()
    }
    
  } else {
    if (ldb == 1L) {
      data_plot <- data.frame(
        LA       = Cb,
        Time     = tk.plot,
        lower_LA = t(ENV.Cb)[, 1L],
        upper_LA = t(ENV.Cb)[, 2L]
      )
      data_plot2 <- data.frame(
        Cb   = LA_fun(1L),
        Time = tk.plot2
      )
      
      Plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = Time, y = LA)) +
        ggplot2::geom_line(colour = arg$colour, linewidth = arg$size) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower_LA, ymax = upper_LA),
          fill        = ci_fill,
          alpha       = ci_alpha,
          colour      = NA,
          show.legend = FALSE
        ) +
        ggplot2::geom_point(
          data  = data_plot2,
          ggplot2::aes(x = Time, y = Cb),
          shape = arg$shape,
          alpha = point_alpha
        ) +
        ggplot2::labs(
          x     = arg$xlab,
          y     = arg$ylab,
          title = level_label(1L)
        ) +
        .lcc_default_theme()
      
    } else {
      data_main <- vector("list", ldb)
      data_la   <- vector("list", ldb)
      for (i in seq_len(ldb)) {
        lab_i <- level_label(i)
        env_i <- ENV.Cb[[i]]
        data_main[[i]] <- data.frame(
          LA       = Cb[, i],
          Time     = tk.plot,
          lower_LA = t(env_i)[, 1L],
          upper_LA = t(env_i)[, 2L],
          Level    = lab_i
        )
        data_la[[i]] <- data.frame(
          Cb    = LA_fun(i),
          Time  = tk.plot2,
          Level = lab_i
        )
      }
      data_plot  <- do.call(rbind.data.frame, data_main)
      data_plot2 <- do.call(rbind.data.frame, data_la)
      
      Plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = Time, y = LA)) +
        ggplot2::geom_line(colour = arg$colour, linewidth = arg$size) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower_LA, ymax = upper_LA),
          fill        = ci_fill,
          alpha       = ci_alpha,
          colour      = NA,
          show.legend = FALSE
        ) +
        ggplot2::geom_point(
          data  = data_plot2,
          ggplot2::aes(x = Time, y = Cb),
          shape = arg$shape,
          alpha = point_alpha
        ) +
        ggplot2::facet_wrap(~Level) +
        ggplot2::labs(x = arg$xlab, y = arg$ylab) +
        .lcc_default_theme()
    }
  }
  
  Plot
}
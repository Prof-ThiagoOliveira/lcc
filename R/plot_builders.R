#######################################################################
# Internal wrappers: plot_lcc / plot_lpc / plot_la
#######################################################################

clamp01 <- function(x) pmin(pmax(x, 0), 1)

##' @keywords internal
plot_lcc <- function(rho, ENV.LCC, tk.plot, tk.plot2, ldb, model, ci, arg,
                     CCC_vals, ...) {
  plotBuilder_lcc(
    rho      = rho,
    ENV.LCC  = ENV.LCC,
    tk.plot  = tk.plot,
    CCC      = CCC_vals,
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
plot_lpc <- function(LPC, ENV.LPC, tk.plot, tk.plot2, ldb, model, ci, arg,
                     Pearson_vals, ...) {
  plotBuilder_lpc(
    LPC          = LPC,
    ENV.LPC      = ENV.LPC,
    tk.plot      = tk.plot,
    Pearson_vals = Pearson_vals,
    tk.plot2     = tk.plot2,
    ldb          = ldb,
    model        = model,
    ci           = ci,
    arg          = arg,
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
plot_la <- function(Cb, ENV.Cb, tk.plot, tk.plot2, ldb, model, ci, arg,
                    CCC_vals, Pearson_vals, ...) {
  plotBuilder_la(
    CCC_vals     = CCC_vals,
    Pearson_vals = Pearson_vals,
    Cb           = Cb,
    ENV.Cb       = ENV.Cb,
    tk.plot      = tk.plot,
    tk.plot2     = tk.plot2,
    ldb          = ldb,
    model        = model,
    ci           = ci,
    arg          = arg,
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
      data_plot  <- data.frame(LCC = clamp01(rho),          Time = tk.plot)
      data_plot2 <- data.frame(CCC = clamp01(CCC[[1L]]$V1), Time = tk.plot2)
      
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
          LCC   = clamp01(rho[, i]),
          Time  = tk.plot,
          Level = lab_i
        )
        data_ccc[[i]] <- data.frame(
          CCC   = clamp01(CCC[[i]]$V1),
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
      if (ncol(ENV.LCC) == 0L) {
        ENV.LCC <- matrix(NA_real_, nrow = 2L, ncol = length(rho))
      }
      data_plot <- data.frame(
        LCC       = clamp01(rho),
        Time      = tk.plot,
        lower_rho = clamp01(t(ENV.LCC)[, 1L]),
        upper_rho = clamp01(t(ENV.LCC)[, 2L])
      )
      data_plot2 <- data.frame(
        CCC  = clamp01(CCC[[1L]]$V1),
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
        if (ncol(env_i) == 0L) {
          env_i <- matrix(NA_real_, nrow = 2L, ncol = nrow(rho))
        }
        data_main[[i]] <- data.frame(
          LCC       = clamp01(rho[, i]),
          Time      = tk.plot,
          lower_rho = clamp01(t(env_i)[, 1L]),
          upper_rho = clamp01(t(env_i)[, 2L]),
          Level     = lab_i
        )
        data_ccc[[i]] <- data.frame(
          CCC   = clamp01(CCC[[i]]$V1),
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
plotBuilder_lpc <- function(LPC, ENV.LPC, tk.plot, Pearson_vals,
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
      data_plot  <- data.frame(LPC = LPC,                      Time = tk.plot)
      data_plot2 <- data.frame(Pearson = Pearson_vals[[1L]]$V1, Time = tk.plot2)
      
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
          Pearson = Pearson_vals[[i]]$V1,
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
        LPC        = clamp01(LPC),
        Time       = tk.plot,
        lower_LPC  = clamp01(t(ENV.LPC)[, 1L]),
        upper_LPC  = clamp01(t(ENV.LPC)[, 2L])
      )
      data_plot2 <- data.frame(
        Pearson = Pearson_vals[[1L]]$V1,
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
          LPC       = clamp01(LPC[, i]),
          Time      = tk.plot,
          lower_LPC = clamp01(t(env_i)[, 1L]),
          upper_LPC = clamp01(t(env_i)[, 2L]),
          Level     = lab_i
        )
        data_pear[[i]] <- data.frame(
          Pearson = Pearson_vals[[i]]$V1,
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
plotBuilder_la <- function(CCC_vals, Pearson_vals, Cb, ENV.Cb,
                           tk.plot, tk.plot2, ldb, model, ci, arg, ...) {
  
  method_levels <- levels(model$data$method)
  level_label <- function(i) {
    paste(method_levels[i + 1L], "vs.", method_levels[1L])
  }
  
  LA_fun <- function(i) CCC_vals[[i]]$V1 / Pearson_vals[[i]]$V1
  
  ci_fill     <- if (!is.null(arg$ci_fill)) arg$ci_fill else arg$colour
  ci_alpha    <- if (!is.null(arg$ci_alpha)) arg$ci_alpha else 0.15
  point_alpha <- if (!is.null(arg$point_alpha)) arg$point_alpha else 0.8
  
  if (!ci) {
    if (ldb == 1L) {
      data_plot  <- data.frame(LA = clamp01(Cb),         Time = tk.plot)
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
            LA    = clamp01(Cb[, i]),
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
#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: plot_builders.R                                               #
# Contains: plot_lcc, plot_lpc, plot_la and plot builders             #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 23/11/2025                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

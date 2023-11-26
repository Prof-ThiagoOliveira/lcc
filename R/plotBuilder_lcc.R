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
# Last update: 03/06/2020                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
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
                            tk.plot2, ldb, model, ci, arg, ...){
  if(ci==FALSE){
    #-------------------------------------------------------------------
    if(ldb == 1) {
      #-----------------------------------------------------------------
      data_plot<-data.frame("LCC"=rho,
                            "Time"=tk.plot)
      data_plot2<-data.frame("CCC"=CCC[[1]]$V1,
                             "Time"=tk.plot2)
      LCC<-data_plot$LCC
      Time<-data_plot$Time
      Plot<-ggplot(data_plot, aes(y=LCC, x=Time))+
        geom_path(data=data_plot, colour=arg$colour, size=arg$size)+
        geom_point(data=data_plot2, aes(y=CCC, x=Time), shape=arg$shape)+
        ggtitle(paste(levels(model$data$method)[2], "vs.",
                      levels(model$data$method)[1]))+
        labs(x = paste0(arg$xlab))+
        labs(y = paste0(arg$ylab))+
        theme(plot.title = element_text(hjust = 0.5))
      print(Plot)
    } else{
      data_plot<-list(NA)
      data_plot2<-list(NA)
      for(i in 1:ldb){
        data_plot[[i]] <- data.frame("LCC"=rho[,i],
                                     "Time"=tk.plot)
        data_plot[[i]]$Level <-
          paste(levels(model$data$method)[i+1], "vs.",
                levels(model$data$method)[1])
        data_plot2[[i]] <- data.frame("CCC"=CCC[[i]]$V1,
                                      "Time"=tk.plot2)
        data_plot2[[i]]$Level <-
          paste(levels(model$data$method)[i+1], "vs.",
                levels(model$data$method)[1])
      }
      data_plot_ <- do.call(rbind.data.frame, data_plot)
      data_plot2_ <- do.call(rbind.data.frame, data_plot2)
      Plot <- ggplot(data_plot_, aes(y=LCC, x=Time))+
        geom_line(data=data_plot_, colour=arg$colour, size=arg$size)+
        geom_point(data=data_plot2_, aes(y=CCC, x=Time), shape=arg$shape)+
        facet_wrap(~ Level, ...) +
        labs(x = paste0(arg$xlab)) +
        labs(y = paste0(arg$ylab)) +
        theme(plot.title = element_text(hjust = 0.5))
      print(Plot)
    }
  }else{
    #-------------------------------------------------------------------
    if(ldb == 1) {
      data_plot <- data.frame("LCC"=rho,
                              "Time"=tk.plot,
                              "lower_rho"=t(ENV.LCC)[,1],
                              "upper_rho"=t(ENV.LCC)[,2])
      data_plot2 <- data.frame("CCC"=CCC[[1]]$V1,
                               "Time"=tk.plot2)
      LCC<-data_plot$LCC
      Time<-data_plot$Time
      lower_rho<-data_plot$lower_rho
      upper_rho<-data_plot$upper_rho
      Plot<-ggplot(data_plot, aes(y=LCC, x=Time))+
        geom_line(data=data_plot, colour=arg$colour, size=arg$size)+
        geom_point(data=data_plot2, aes(y=CCC, x=Time), shape=arg$shape)+
        geom_ribbon(data=data_plot,aes(ymin=lower_rho,ymax=upper_rho),
                    fill="grey70", alpha=0.3,show.legend = TRUE)+
        ggtitle(paste(levels(model$data$method)[2], "vs.",
                      levels(model$data$method)[1]))+
        labs(x = paste0(arg$xlab))+
        labs(y = paste0(arg$ylab))+
        theme(plot.title = element_text(hjust = 0.5))
      print(Plot)
    } else{
      data_plot<-list(NA)
      data_plot2<-list(NA)
      for(i in 1:ldb){
        data_plot[[i]] <- data.frame("LCC"=rho[,i],
                                     "Time"=tk.plot,
                                     "lower_rho"=t(ENV.LCC[[i]])[,1],
                                     "upper_rho"=t(ENV.LCC[[i]])[,2])
        data_plot[[i]]$Level <-
          paste(levels(model$data$method)[i+1], "vs.",
                levels(model$data$method)[1])
        data_plot2[[i]] <- data.frame("CCC"=CCC[[i]]$V1,
                                      "Time"=tk.plot2)
        data_plot2[[i]]$Level <-
          paste(levels(model$data$method)[i+1], "vs.",
                levels(model$data$method)[1])
      }
      data_plot_ <- do.call(rbind.data.frame, data_plot)
      data_plot2_ <- do.call(rbind.data.frame, data_plot2)
      Plot <- ggplot(data_plot_, aes(y=LCC, x=Time))+
        geom_line(data=data_plot_, colour=arg$colour, size=arg$size)+
        geom_point(data=data_plot2_, aes(y=CCC, x=Time), shape=arg$shape)+
        geom_ribbon(data=data_plot_, aes(ymin = lower_rho,
                                         ymax = upper_rho),
                    fill="grey70", alpha = 0.3, show.legend = TRUE) +
        facet_wrap(~Level, ...) +
        labs(x = paste0(arg$xlab)) +
        labs(y = paste0(arg$ylab)) +
        theme(plot.title = element_text(hjust = 0.5))
      print(Plot)
    }
  }
}
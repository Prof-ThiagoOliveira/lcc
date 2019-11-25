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
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal Function to Produces a Longitudinal Concordance
##'   Correlation Plot.
##'
##' @description This is an internally called function used to produces
##'   a longitudinal concordance correlation plot from fitted ans
##'   sampled values with or not non-parametric confidence intervals.
##'
##' @details returns a inital plot for the longitudinal concordance
##'   correlation.
##'
##' @usage NULL
##'
##' @importFrom grid grid.newpage
##'
##' @importFrom gridExtra marrangeGrob
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@usp.br}
##'
##' @keywords internal
plotBuilder_lcc<-function(rho, ENV.LCC, tk.plot, CCC,
                           tk.plot2, ldb, model, ci, arg){
if(ci==FALSE){
  if(ldb == 1) {
    data_plot<-data.frame("LCC"=rho,
                          "Time"=tk.plot)
    data_plot2<-data.frame("CCC"=CCC[[1]]$V1,
                           "Time"=tk.plot2)
    LCC<-data_plot$LCC
    Time<-data_plot$Time
    Plot<-ggplot(data_plot, aes(y=LCC, x=Time))+
      geom_path(data=data_plot, colour=arg$colour, size=arg$size)+
      geom_point(data=data_plot2, aes(y=CCC, x=Time), shape=arg$shape)+
      scale_y_continuous(limits = arg$scale_y_continuous)+
      ggtitle(paste(levels(model$data$method)[2], "vs.",
                    levels(model$data$method)[1]))+
      labs(x = paste0(arg$xlab))+
      labs(y = paste0(arg$ylab))+
      theme(plot.title = element_text(hjust = 0.5))
    if(arg$scale_y_continuous[2]==1){
      Plot<-Plot+geom_hline(yintercept = 1, linetype="dashed")
    }
    print(Plot)
  } else{
    data_plot<-list(NA)
    data_plot2<-list(NA)
    Plot<-list(NA)
    for(i in 1:ldb){
      data_plot[[i]]<-data.frame("LCC"=rho[,i],
                                 "Time"=tk.plot)
      data_plot2[[i]]<-data.frame("CCC"=CCC[[i]]$V1,
                                  "Time"=tk.plot2)
      LCC<-data_plot[[i]]$LCC
      Time<-data_plot[[i]]$Time
      Plot[[i]]<-ggplot(data_plot[[i]], aes(y=LCC, x=Time))+
        geom_line(data=data_plot[[i]], colour=arg$colour, size=arg$size)+
        geom_point(data=data_plot2[[i]], aes(y=CCC, x=Time), shape=arg$shape)+
        scale_y_continuous(limits = arg$scale_y_continuous)+
        ggtitle(paste(levels(model$data$method)[i+1], "vs.",
                      levels(model$data$method)[1]))+
        labs(x = paste0(arg$xlab))+
        labs(y = paste0(arg$ylab))+
        theme(plot.title = element_text(hjust = 0.5))
        if(arg$scale_y_continuous[2]==1){
          Plot[[i]]<-Plot[[i]]+geom_hline(yintercept = 1, linetype="dashed")
        }
    }
    numPlots = length(Plot)
    if(arg$all.plot){
    grid.newpage()
    cols<-ceiling(sqrt(ldb))
    rows<-signif(sqrt(ldb),1)
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = rows)
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(Plot[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                     layout.pos.col = matchidx$col))
    }
    }else{
      all_plots <- lapply(1:numPlots, function(x) Plot[[x]])
      ml <- gridExtra::marrangeGrob(all_plots, nrow = 1, ncol = 1,
                                    top = " ")
      print(ml)
    }
  }
}else{
if(ldb == 1) {
  data_plot<-data.frame("LCC"=rho,
                        "Time"=tk.plot,
                        "lower_rho"=t(ENV.LCC)[,1],
                        "upper_rho"=t(ENV.LCC)[,2])
  data_plot2<-data.frame("CCC"=CCC[[1]]$V1,
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
    scale_y_continuous(limits = arg$scale_y_continuous)+
    ggtitle(paste(levels(model$data$method)[2], "vs.",
                  levels(model$data$method)[1]))+
    labs(x = paste0(arg$xlab))+
    labs(y = paste0(arg$ylab))+
    theme(plot.title = element_text(hjust = 0.5))
    if(arg$scale_y_continuous[2]==1){
      Plot<-Plot+geom_hline(yintercept = 1, linetype="dashed")
    }
    print(Plot)
} else{
  data_plot<-list(NA)
  data_plot2<-list(NA)
  Plot<-list(NA)
  for(i in 1:ldb){
    data_plot[[i]]<-data.frame("LCC"=rho[,i],
                               "Time"=tk.plot,
                               "lower_rho"=t(ENV.LCC[[i]])[,1],
                               "upper_rho"=t(ENV.LCC[[i]])[,2])
    data_plot2[[i]]<-data.frame("CCC"=CCC[[i]]$V1,
                                "Time"=tk.plot2)
    LCC<-data_plot[[i]]$LCC
    Time<-data_plot[[i]]$Time
    lower_rho<-data_plot[[i]]$lower_rho
    upper_rho<-data_plot[[i]]$upper_rho
    Plot[[i]]<-ggplot(data_plot[[i]], aes(y=LCC, x=Time))+
      geom_line(data=data_plot[[i]], colour=arg$colour, size=arg$size)+
      geom_point(data=data_plot2[[i]], aes(y=CCC, x=Time), shape=arg$shape)+
      geom_ribbon(data=data_plot[[i]],aes(ymin=lower_rho,ymax=upper_rho),
                  fill="grey70", alpha=0.3,show.legend = TRUE)+
      scale_y_continuous(limits = arg$scale_y_continuous)+
      ggtitle(paste(levels(model$data$method)[i+1], "vs.",
                    levels(model$data$method)[1]))+
      labs(x = paste0(arg$xlab))+
      labs(y = paste0(arg$ylab))+
      theme(plot.title = element_text(hjust = 0.5))
      if(arg$scale_y_continuous[2]==1){
        Plot[[i]]<-Plot[[i]]+geom_hline(yintercept = 1, linetype="dashed")
      }
  }
  numPlots = length(Plot)
  if(arg$all.plot){
  grid.newpage()
  cols<-ceiling(sqrt(ldb))
  rows<-signif(sqrt(ldb),1)
  layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                   ncol = cols, nrow = rows)
  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
  for (i in 1:numPlots) {
    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
    print(Plot[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                   layout.pos.col = matchidx$col))
        }
  }else{
     all_plots <- lapply(1:numPlots, function(x) Plot[[x]])
     ml <- gridExtra::marrangeGrob(all_plots, nrow = 1, ncol = 1,
                                   top = " ")
      print(ml)
      }
    }
  }
}

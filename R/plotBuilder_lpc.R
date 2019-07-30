#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: plotBuilder_lpc.R                                             #
# Contains: plotBuilder_lpc function                                  #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal function to produces a longitudinal Perason
##'   correlation plot.
##'
##' @description This is an internally called function used to produces
##'   a longitudinal Perason correlation plot from fitted ans sampled
##'   values with or not non-parametric confidence intervals.
##'
##' @details returns a inital plot for the longitudinal Pearson
##'   correlation.
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@usp.br}
##'
##' @usage NULL
##'
##' @importFrom grid pushViewport grid.layout
##'
##' @keywords internal
plotBuilder_lpc<-function(LPC, ENV.LPC, tk.plot, Pearson,
                           tk.plot2, ldb, model, ci, arg){
  Time<-tk.plot
if(ci==FALSE){
  if(ldb == 1) {
    data_plot<-data.frame("LPC"=LPC,
                          "Time"=tk.plot)
    data_plot2<-data.frame("Pearson"=Pearson[[1]]$V1,
                           "Time"=tk.plot2)
    LPC<-data_plot$LPC
    Time<-data_plot$Time
    Plot<-ggplot(data_plot, aes(y=LPC, x=Time))+
      geom_line(data=data_plot, colour=arg$colour, size=arg$size)+
      geom_point(data=data_plot2, aes(y=Pearson, x=Time),
                 shape=arg$shape)+
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
      data_plot[[i]]<-data.frame("LPC"=LPC[,i],
                                 "Time"=tk.plot)

      data_plot2[[i]]<-data.frame("Pearson"=Pearson[[i]]$V1,
                                  "Time"=tk.plot2)
      LPC2<-data_plot[[i]]$LPC
      Time<-data_plot[[i]]$Time
      Plot[[i]]<-ggplot(data_plot[[i]], aes(y=LPC2, x=Time))+
        geom_line(data=data_plot[[i]], colour=arg$colour, size=arg$size)+
        geom_point(data=data_plot2[[i]], aes(y=Pearson, x=Time),
                   shape=arg$shape)+
        scale_y_continuous(limits = arg$scale_y_continuous)+
        ggtitle(paste(levels(model$data$method)[i+1], "vs.",
                      levels(model$data$method)[1]))+
        labs(x = paste0(arg$xlab))+
        labs(y = paste0(arg$ylab))+
        theme(plot.title = element_text(hjust = 0.5))
      if(arg$scale_y_continuous[2]==1){
        Plot[[i]]<-Plot[[i]]+geom_hline(yintercept = 1,
                                        linetype="dashed")
      }
    }
    numPlots = length(Plot)
    if(arg$all.plot){
    grid.newpage()
    cols<-ceiling(sqrt(ldb))
    rows<-signif(sqrt(ldb),1)
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = rows)
    pushViewport(viewport(layout = grid.layout(nrow(layout),
                                               ncol(layout))))
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
  data_plot<-data.frame("LPC"=LPC,
                        "Time"=tk.plot,
                        "lower_LPC"=t(ENV.LPC)[,1],
                        "upper_LPC"=t(ENV.LPC)[,2]
  )
  data_plot2<-data.frame("Pearson"=Pearson[[1]]$V1,
                         "Time"=tk.plot2)
  LPC<-data_plot$LPC
  Time<-data_plot$Time
  lower_LPC<-data_plot$lower_LPC
  upper_LPC<-data_plot$upper_LPC
  Plot<-ggplot(data_plot, aes(y=LPC, x=Time))+
    geom_line(data=data_plot, colour=arg$colour, size=arg$size)+
    geom_point(data=data_plot2, aes(y=Pearson, x=Time),
               shape=arg$shape)+
    geom_ribbon(data=data_plot,aes(ymin=lower_LPC,ymax=upper_LPC),
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
    data_plot[[i]]<-data.frame("LPC"=LPC[,i],
                               "Time"=tk.plot,
                               "lower_LPC"=t(ENV.LPC[[i]])[,1],
                               "upper_LPC"=t(ENV.LPC[[i]])[,2])

    data_plot2[[i]]<-data.frame("Pearson"=Pearson[[i]]$V1,
                                "Time"=tk.plot2)
    LPC2<-data_plot[[i]]$LPC
    Time<-data_plot[[i]]$Time
    lower_LPC<-data_plot[[i]]$lower_LPC
    upper_LPC<-data_plot[[i]]$upper_LPC
    Plot[[i]]<-ggplot(data_plot[[i]], aes(y=LPC2, x=Time))+
      geom_line(data=data_plot[[i]], colour=arg$colour, size=arg$size)+
      geom_point(data=data_plot2[[i]], aes(y=Pearson, x=Time),
                 shape=arg$shape)+
      geom_ribbon(data=data_plot[[i]],aes(ymin=lower_LPC,ymax=upper_LPC),
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

#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: plotBuilder_la.R                                              #
# Contains: plotBuilder_la function                                   #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal function to produces a longitudinal accuracy plot.
##'
##' @description This is an internally called function used to produces
##'   a longitudinal accuracy plot from fitted ans sampled values with
##'   or not non-parametric confidence intervals.
##'
##' @details returns a inital plot for the longitudinal accuracy
##'   correlation.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@usp.br}
##'
##' @keywords internal
plotBuilder_la<-function(CCC, Pearson, Cb, ENV.Cb, tk.plot,
                       tk.plot2, ldb, model, ci, arg){
if(ci==FALSE){
  if(ldb == 1) {
    data_plot<-data.frame("LA"=Cb,
                          "Time"=tk.plot)
    data_plot2<-data.frame("Cb"=CCC[[1]]$V1/Pearson[[1]]$V1,
                           "Time"=tk.plot2)
    LA<-data_plot$Cb
    Time<-data_plot$Time
    Plot<-ggplot(data_plot, aes(y=LA, x=Time))+
      geom_line(data=data_plot, colour=arg$colour, size=arg$size)+
      geom_point(data=data_plot2, aes(y=Cb, x=Time), shape=arg$shape)+
      scale_y_continuous(limits = arg$scale_y_continuous)+
<<<<<<< HEAD
      ggtitle(paste(levels(model$data$method)[2], "vs.",
                    levels(model$data$method)[1]))+
=======
      ggtitle(paste(levels(model$data$method)[2], "vs.", levels(model$data$method)[1]))+
>>>>>>> 0defa447a19649cc6f57a67efad6946bd7b187aa
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
      data_plot[[i]]<-data.frame("LA"=Cb[,i],
                                 "Time"=tk.plot)

      data_plot2[[i]]<-data.frame("Cb"=CCC[[i]]$V1/Pearson[[i]]$V1,
                                  "Time"=tk.plot2)
      LA<-data_plot[[i]]$Cb
      Time<-data_plot[[i]]$Time

      Plot[[i]]<-ggplot(data_plot[[i]], aes(y=LA, x=Time))+
<<<<<<< HEAD
        geom_line(data=data_plot[[i]], colour=arg$colour,
                  size=arg$size)+
        geom_point(data=data_plot2[[i]], aes(y=Cb, x=Time),
                   shape=arg$shape)+
        scale_y_continuous(limits = arg$scale_y_continuous)+
        ggtitle(paste(levels(model$data$method)[i+1], "vs.",
                      levels(model$data$method)[1]))+
=======
        geom_line(data=data_plot[[i]], colour=arg$colour, size=arg$size)+
        geom_point(data=data_plot2[[i]], aes(y=Cb, x=Time), shape=arg$shape)+
        scale_y_continuous(limits = arg$scale_y_continuous)+
        ggtitle(paste(levels(model$data$method)[i+1], "vs.", levels(model$data$method)[1]))+
>>>>>>> 0defa447a19649cc6f57a67efad6946bd7b187aa
        labs(x = paste0(arg$xlab))+
        labs(y = paste0(arg$ylab))+
        theme(plot.title = element_text(hjust = 0.5))
      if(arg$scale_y_continuous[2]==1){
<<<<<<< HEAD
        Plot[[i]]<-Plot[[i]]+geom_hline(yintercept = 1,
                                        linetype="dashed")
=======
        Plot[[i]]<-Plot[[i]]+geom_hline(yintercept = 1, linetype="dashed")
>>>>>>> 0defa447a19649cc6f57a67efad6946bd7b187aa
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
      invisible(print(ml))
    }
  }
}else{
if(ldb == 1) {
  data_plot<-data.frame("LA"=Cb,
                        "Time"=tk.plot,
                        "lower_LA"=t(ENV.Cb)[,1],
                        "upper_LA"=t(ENV.Cb)[,2]
  )
  data_plot2<-data.frame("Cb"=CCC[[1]]$V1/Pearson[[1]]$V1,
                         "Time"=tk.plot2)
  LA<-data_plot$Cb
  Time<-data_plot$Time
  lower_LA<-data_plot$lower_LA
  upper_LA<-data_plot$upper_LA

  Plot<-ggplot(data_plot, aes(y=LA, x=Time))+
    geom_line(data=data_plot, colour=arg$colour, size=arg$size)+
    geom_point(data=data_plot2, aes(y=Cb, x=Time), shape=arg$shape)+
    geom_ribbon(data=data_plot,aes(ymin=lower_LA,ymax=upper_LA),
                fill="grey70", alpha=0.3,show.legend = TRUE)+
    scale_y_continuous(limits = arg$scale_y_continuous)+
<<<<<<< HEAD
    ggtitle(paste(levels(model$data$method)[2], "vs.",
                  levels(model$data$method)[1]))+
=======
    ggtitle(paste(levels(model$data$method)[2], "vs.", levels(model$data$method)[1]))+
>>>>>>> 0defa447a19649cc6f57a67efad6946bd7b187aa
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
    data_plot[[i]]<-data.frame("LA"=Cb[,i],
                               "Time"=tk.plot,
                               "lower_LA"=t(ENV.Cb[[i]])[,1],
                               "upper_LA"=t(ENV.Cb[[i]])[,2])

    data_plot2[[i]]<-data.frame("Cb"=CCC[[i]]$V1/Pearson[[i]]$V1,
                                "Time"=tk.plot2)
    LA<-data_plot[[i]]$Cb
    Time<-data_plot[[i]]$Time
    lower_LA<-data_plot[[i]]$lower_LA
    upper_LA<-data_plot[[i]]$upper_LA

    Plot[[i]]<-ggplot(data_plot[[i]], aes(y=LA, x=Time))+
      geom_line(data=data_plot[[i]], colour=arg$colour, size=arg$size)+
      geom_point(data=data_plot2[[i]], aes(y=Cb, x=Time),
                 shape=arg$shape)+
      geom_ribbon(data=data_plot[[i]],aes(ymin=lower_LA,ymax=upper_LA),
                  fill="grey70", alpha=0.3,show.legend = TRUE)+
      scale_y_continuous(limits = arg$scale_y_continuous)+
<<<<<<< HEAD
      ggtitle(paste(levels(model$data$method)[i+1], "vs.",
                    levels(model$data$method)[1]))+
=======
      ggtitle(paste(levels(model$data$method)[i+1], "vs.", levels(model$data$method)[1]))+
>>>>>>> 0defa447a19649cc6f57a67efad6946bd7b187aa
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
      invisible(print(ml))
      }
    }
  }
}

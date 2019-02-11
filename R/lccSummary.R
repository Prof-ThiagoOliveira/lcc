#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: lccSummary.R                                                  #
# Contains: lccSummary function                                       #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 18/06/2018                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################
##' @title Internal function to summarize fitted and sampled values for \code{lcc} objects
##'
##' @description This is an internally called function used to summarize fitted and sampled values, and the concordance correlation coefficient between them for \code{lcc} objects.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@usp.br}
##'
##' @keywords internal
lccSummary<-function(model, q_f, diffbeta, tk,
                      tk.plot, tk.plot2, rho, ENV.LCC,
                      rho.pearson, ENV.LPC, Cb, ENV.Cb,
                      ldb, ci, components){
  if(components==FALSE){
    if(ci==FALSE){
    CCC<-CCC_lin(dataset=model$data, resp="y", subject="ind", method="FacA", time="time")
    if(ldb==1){
        LCC.data<-data.frame("Time"=tk.plot,"LCC"=rho)
        names(LCC.data)<-c("Time", paste0(expression("LCC: "),levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]))
        CCC.data<-data.frame(tk.plot2, CCC)
        names(CCC.data)<-c("Time", paste0(expression("CCC: "),levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]))
        if(length(CCC.data[,2])==length(LCC.data[,2])){
          GF<-CCC(CCC.data[,2],LCC.data[,2])
        }else{
          rho_gf <- lccWrapper(model = model, q_f = q_f, n.delta = 1,
                             tk = tk, diffbeta = as.numeric(diffbeta[[1]]))
          GF<-CCC(CCC.data[,2],rho_gf)
        }
        plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, "gof" = GF)
    }else{
      LCC.data<-list()
      GF<-list()
      rho_gf <- list()
      for(i in 1:ldb) {
        LCC.data[[i]]<-data.frame("Time"=tk.plot,"LCC"=rho[[i]])
        names(LCC.data[[i]])<-c("Time", paste0(expression("LCC: "), levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]))
        names(CCC[[i]])<-c(paste0(expression("CCC: "), levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]))
        CCC.data<-data.frame(tk.plot2, CCC)
        if(length(CCC.data[,2])==length(LCC.data[[i]][,2])){
        GF[[i]]<-CCC(as.numeric(unlist(CCC[[i]])),LCC.data[[i]][,2])
        names(GF[[i]])<-c(paste0(levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]))
        }else{
          rho_gf[[i]] <- lccWrapper(model = model, q_f = q_f, n.delta = 1,
                                                   tk = tk, diffbeta = as.numeric(diffbeta[[i]]))
          GF[[i]]<-CCC(as.numeric(unlist(CCC[[i]])),rho_gf[[i]])
          names(GF[[i]])<-c(paste0(levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]))
        }
      }
      plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, "gof" = GF)
    }
  }else{
    if(ldb==1){
    CCC<-CCC_lin(dataset=model$data, resp="y", subject="ind", method="FacA", time="time")
    LCC.data<-data.frame("Time"=tk.plot,"LCC"=rho, "Lower"=ENV.LCC[1,], "Upper"=ENV.LCC[2,])
      names(LCC.data)<-c("Time", paste0(expression("LCC: "), levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]), "Lower", "Upper")
    CCC.data<-data.frame(tk.plot2, CCC)
      names(CCC.data)<-c("Time", paste0(expression("CCC: "), levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]))
      if(length(CCC.data[,2])==length(LCC.data[,2])){
        GF<-CCC(CCC.data[,2],LCC.data[,2])
      }else{
        rho_gf <- lccWrapper(model = model, q_f = q_f, n.delta = 1,
                              tk = tk, diffbeta = as.numeric(diffbeta[[1]]))
        GF<-CCC(CCC.data[,2],rho_gf)
      }
       plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, "gof" = GF)
    }else{
      CCC<-CCC_lin(dataset=model$data, resp="y", subject="ind", method="FacA", time="time")
      LCC.data<-list()
      GF<-list()
      rho_gf<-list()
      for(i in 1:ldb) {
        LCC.data[[i]]<-data.frame("Time"=tk.plot,"LCC"=rho[[i]], "Lower"=ENV.LCC[[i]][1,], "Upper"=ENV.LCC[[i]][2,])
        names(LCC.data[[i]])<-c("Time", paste0(expression("LCC: "), levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]), "Lower", "Upper")
        names(CCC[[i]])<-c(paste0(expression("CCC: "), levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]))
        CCC.data<-data.frame(tk.plot2, CCC)
        if(length(CCC.data[,2])==length(LCC.data[[i]][,2])){
          GF[[i]]<-CCC(as.numeric(unlist(CCC[[i]])),LCC.data[[i]][,2])
          names(GF[[i]])<-c(paste0(levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]))
        }else{
          rho_gf[[i]] <- lccWrapper(model = model, q_f = q_f, n.delta = 1,
                                     tk = tk, diffbeta = as.numeric(diffbeta[[i]]))
          GF[[i]]<-CCC(as.numeric(unlist(CCC[[i]])),rho_gf[[i]])
          names(GF[[i]])<-c(paste0(levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]))
        }
      }
      plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, "gof" = GF)
      }
    }
  }else{
    if(ci==FALSE){
      CCC<-CCC_lin(dataset=model$data, resp="y", subject="ind", method="FacA", time="time")
      Pearson<-Pearson(dataset=model$data, resp="y", subject="ind", method="FacA", time="time")
       if(ldb==1){
        LA<-CCC[[1]]/Pearson[[1]]         
        LCC.data<-data.frame("Time"=tk.plot,"LCC"=rho, "LPC"=rho.pearson, "LA"=Cb)
        names(LCC.data)<-c("Time",
                           paste0(expression("LCC: "), levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]),
                           paste0(expression("LPC: "),levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]),
                           paste0(expression("LA: "), levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]))
        CCC.data<-data.frame(tk.plot2, CCC, Pearson, LA)
        names(CCC.data)<-c("Time", paste0(expression("CCC: "), levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]),
                           paste0(expression("Pearson: "), levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]),
                           paste0(expression("Cb: "), levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]))
        if(length(CCC.data[,2])==length(LCC.data[,2])){
          GF<-CCC(CCC.data[,2],LCC.data[,2])
        }else{
          rho_gf <- lccWrapper(model = model, q_f = q_f, n.delta = 1,
                                tk = tk, diffbeta = as.numeric(diffbeta[[1]]))
          GF<-CCC(CCC.data[,2],rho_gf)
        }
        plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, "gof" = GF)
      }else{
        LCC.data<-list()
        CCC.data<-list()
        LA<-list()
        GF<-list()
        rho_gf<-list()
        for(i in 1:ldb) {
          LA[[i]]<-CCC[[i]]/Pearson[[i]]
          LCC.data[[i]]<-data.frame("Time"=tk.plot,"LCC"=rho[[i]], "LPC"=rho.pearson[[i]], "LA"=Cb[[i]])
          names(LCC.data[[i]])<-c("Time",
                                  paste0(expression("LCC: "),levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]),
                                  paste0(expression("LPC: "),levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]),
                                  paste0(expression("LA: "),levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]))
          CCC.data[[i]]<-data.frame(tk.plot2, CCC[[i]], Pearson[[i]], LA[[i]])
          names(CCC.data[[i]])<-c("Time", paste0(expression("CCC: "),levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]),
                                          paste0(expression("Pearson: "),levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]),
                                          paste0(expression("Cb: "),levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]))
          if(length(CCC.data[[i]][,2])==length(LCC.data[[i]][,2])){
            GF[[i]]<-CCC(as.numeric(unlist(CCC[[i]])),LCC.data[[i]][,2])
            names(GF[[i]])<-c(paste0(levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]))
          }else{
            rho_gf[[i]] <- lccWrapper(model = model, q_f = q_f, n.delta = i,
                                       tk = tk, diffbeta = as.numeric(diffbeta[[i]]))
            GF[[i]]<-CCC(as.numeric(unlist(CCC[[i]])),rho_gf[[i]])
            names(GF[[i]])<-c(paste0(levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]))
          }
        }
        plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, "gof" = GF)
      }
    }else{
      if(ldb==1){
        CCC<-CCC_lin(dataset=model$data, resp="y", subject="ind", method="FacA", time="time")
        Pearson<-Pearson(dataset=model$data, resp="y", subject="ind", method="FacA", time="time")
        LA<-CCC[[1]]/Pearson[[1]]
        LCC.data<-data.frame("Time"=tk.plot,"LCC"=rho, "Lower"=ENV.LCC[1,], "Upper"=ENV.LCC[2,])
        LPC.data<-data.frame("Time"=tk.plot,"LCC"=rho.pearson, "Lower"=ENV.LPC[1,], "Upper"=ENV.LPC[2,])
        LA.data<-data.frame("Time"=tk.plot,"LA"=Cb, "Lower"=ENV.Cb[1,], "Upper"=ENV.Cb[2,])
        names(LCC.data)<-c("Time", paste0(expression("LCC: "),levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]), "Lower", "Upper")
        names(LPC.data)<-c("Time", paste0(expression("LPC: "),levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]), "Lower", "Upper")
        names(LA.data)<-c("Time", paste0(expression("LA: "), levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]), "Lower", "Upper")
        CCC.data<-data.frame(tk.plot2, CCC, Pearson, LA)
        names(CCC.data)<-c("Time", paste0(expression("CCC: "), levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]),
                                   paste0(expression("Pearson: "), levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]),
                                   paste0(expression("Cb: "), levels(model$data$FacA)[2], " vs. ", levels(model$data$FacA)[1]))
        fit<-list("LCC" = LCC.data, "LPC" = LPC.data, "LA" = LA.data)
        if(length(CCC.data[,2])==length(LCC.data[,2])){
          GF<-CCC(CCC.data[,2],LCC.data[,2])
        }else{
          rho_gf <- lccWrapper(model = model, q_f = q_f, n.delta = 1,
                                tk = tk, diffbeta = as.numeric(diffbeta[[1]]))
          GF<-CCC(CCC.data[,2],rho_gf)
        }
        plot.data<-list("fitted"=fit,"sampled" = CCC.data, "gof" = GF)
      }else{
        CCC<-CCC_lin(dataset=model$data, resp="y", subject="ind", method="FacA", time="time")
        Pearson<-Pearson(dataset=model$data, resp="y", subject="ind", method="FacA", time="time")
        LA<-list()
        CCC.data<-list()
        LCC.data<-list()
        LPC.data<-list()
        LA.data<-list()
        GF<-list()
        rho_gf<-list()
        for(i in 1:ldb) {
          LA[[i]]<-CCC[[i]]/Pearson[[i]]
          LCC.data[[i]]<-data.frame("Time"=tk.plot,"LCC"=rho[[i]], "Lower"=ENV.LCC[[i]][1,], "Upper"=ENV.LCC[[i]][2,])
          names(LCC.data[[i]])<-c("Time", paste0(expression("LCC: "), levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]), "Lower", "Upper")
          LPC.data[[i]]<-data.frame("Time"=tk.plot,"LPC"=rho.pearson[[i]], "Lower"=ENV.LPC[[i]][1,], "Upper"=ENV.LPC[[i]][2,])
          names(LPC.data[[i]])<-c("Time", paste0(expression("LPC: "), levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]), "Lower", "Upper")
          LA.data[[i]]<-data.frame("Time"=tk.plot,"LA"=Cb[[i]], "Lower"=ENV.Cb[[i]][1,], "Upper"=ENV.Cb[[i]][2,])
          names(LA.data[[i]])<-c("Time", paste0(expression("LA: "), levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]), "Lower", "Upper")
          CCC.data[[i]]<-data.frame(tk.plot2, CCC[[i]], Pearson[[i]], LA[[i]])
          names(CCC.data[[i]])<-c("Time", paste0(expression("CCC: "), levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]),
                                         paste0(expression("Pearson: "), levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]),
                                         paste0(expression("Cb: "), levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]))
          if(length(CCC.data[[i]][,2])==length(LCC.data[[i]][,2])){
            GF[[i]]<-CCC(as.numeric(unlist(CCC[[i]])),LCC.data[[i]][,2])
            names(GF[[i]])<-c(paste0(levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]))
          }else{
            rho_gf[[i]] <- lccWrapper(model = model, q_f = q_f, n.delta = i,
                                       tk = tk, diffbeta = as.numeric(diffbeta[[i]]))
            GF[[i]]<-CCC(as.numeric(unlist(CCC[[i]])),rho_gf[[i]])
            names(GF[[i]])<-c(paste0(levels(model$data$FacA)[i+1], " vs. ", levels(model$data$FacA)[1]))
          }
        }
        fit<-list("LCC" = LCC.data, "LPC" = LPC.data, "LA" = LA.data)
        plot.data<-list("fitted"=fit, "sampled" = CCC.data, "gof" = GF)
      }
    }
  }
  return(invisible(plot.data))
}

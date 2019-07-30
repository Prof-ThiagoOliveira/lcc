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
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal function to summarize fitted and sampled values for
##'   \code{lcc} objects
##'
##' @description This is an internally called function used to summarize
##'   fitted and sampled values, and the concordance correlation
##'   coefficient between them for \code{lcc} objects.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@usp.br}
##'
##' @importFrom stats predict
##'
##' @keywords internal
lccSummary<-function(model, q_f, diffbeta, tk,
                      tk.plot, tk.plot2, rho, ENV.LCC,
                      rho.pearson, ENV.LPC, Cb, ENV.Cb,
                      ldb, ci, components){
  if(components==FALSE){
    if(ci==FALSE){
<<<<<<< HEAD
      CCC<-CCC_lin(dataset=model$data, resp="resp", subject="subject",
                   method="method", time="time")
=======
    CCC<-CCC_lin(dataset=model$data, resp="resp", subject="subject", method="method", time="time")
>>>>>>> 0defa447a19649cc6f57a67efad6946bd7b187aa
    if(ldb==1){
        LCC.data<-data.frame("Time"=tk.plot,"LCC"=rho)
        names(LCC.data)<-
          c("Time", paste0(expression("LCC: "),
                           levels(model$data$method)[2], " vs. ",
                           levels(model$data$method)[1]))
        CCC.data<-data.frame(tk.plot2, CCC)
        names(CCC.data)<-
          c("Time", paste0(expression("CCC: "),
                           levels(model$data$method)[2],
                           " vs. ", levels(model$data$method)[1]))
        GF<-CCC(predict(model), model$data$resp)
<<<<<<< HEAD
        plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data,
                        "gof" = GF)
=======
        plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, "gof" = GF)
>>>>>>> 0defa447a19649cc6f57a67efad6946bd7b187aa
    }else{
      LCC.data<-list()
      for(i in 1:ldb) {
        LCC.data[[i]]<-data.frame("Time"=tk.plot,"LCC"=rho[[i]])
        names(LCC.data[[i]])<-
          c("Time", paste0(expression("LCC: "),
                           levels(model$data$method)[i+1],
                           " vs. ", levels(model$data$method)[1]))
        names(CCC[[i]])<-
          c(paste0(expression("CCC: "),
                   levels(model$data$method)[i+1], " vs. ",
                   levels(model$data$method)[1]))
        CCC.data<-data.frame(tk.plot2, CCC)
      }
      GF<-CCC(predict(model), model$data$resp)
      plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, "gof" = GF)
    }
  }else{
    if(ldb==1){
      CCC<-CCC_lin(dataset=model$data, resp="resp",
                   subject="subject", method="method", time="time")
      LCC.data<-data.frame("Time"=tk.plot,"LCC"=rho,
                           "Lower"=ENV.LCC[1,], "Upper"=ENV.LCC[2,])
      names(LCC.data)<-
        c("Time", paste0(expression("LCC: "),
                         levels(model$data$method)[2],
                         " vs. ", levels(model$data$method)[1]),
          "Lower", "Upper")
    CCC.data<-data.frame(tk.plot2, CCC)
      names(CCC.data)<-
        c("Time", paste0(expression("CCC: "),
                         levels(model$data$method)[2],
                         " vs. ", levels(model$data$method)[1]))
      GF<-CCC(predict(model), model$data$resp)
       plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, "gof" = GF)
    }else{
      CCC<-CCC_lin(dataset=model$data, resp="resp",
                   subject="subject", method="method", time="time")
      LCC.data<-list()
      for(i in 1:ldb) {
        LCC.data[[i]]<-data.frame("Time"=tk.plot,"LCC"=rho[[i]],
                                  "Lower"=ENV.LCC[[i]][1,],
                                  "Upper"=ENV.LCC[[i]][2,])
        names(LCC.data[[i]])<-
          c("Time", paste0(expression("LCC: "),
                           levels(model$data$method)[i+1],
                           " vs. ", levels(model$data$method)[1]),
            "Lower", "Upper")
        names(CCC[[i]])<-
          c(paste0(expression("CCC: "),
                   levels(model$data$method)[i+1],
                   " vs. ", levels(model$data$method)[1]))
        CCC.data<-data.frame(tk.plot2, CCC)
      }
      GF<-CCC(predict(model), model$data$resp)
      plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, "gof" = GF)
      }
    }
  }else{
    if(ci==FALSE){
      CCC<-CCC_lin(dataset=model$data, resp="resp",
                   subject="subject", method="method", time="time")
      Pearson<-Pearson(dataset=model$data, resp="resp",
                       subject="subject", method="method", time="time")
       if(ldb==1){
        LA<-CCC[[1]]/Pearson[[1]]
        LCC.data<-data.frame("Time"=tk.plot,"LCC"=rho,
                             "LPC"=rho.pearson, "LA"=Cb)
        names(LCC.data)<-
          c("Time", paste0(expression("LCC: "),
                           levels(model$data$method)[2], " vs. ",
                           levels(model$data$method)[1]),
            paste0(expression("LPC: "),levels(model$data$method)[2],
                   " vs. ", levels(model$data$method)[1]),
            paste0(expression("LA: "), levels(model$data$method)[2],
                   " vs. ", levels(model$data$method)[1]))
        CCC.data<-data.frame(tk.plot2, CCC, Pearson, LA)
        names(CCC.data)<-
          c("Time", paste0(expression("CCC: "),
                           levels(model$data$method)[2], " vs. ",
                           levels(model$data$method)[1]),
            paste0(expression("Pearson: "),
                   levels(model$data$method)[2], " vs. ",
                   levels(model$data$method)[1]),
            paste0(expression("Cb: "), levels(model$data$method)[2],
                   " vs. ", levels(model$data$method)[1]))
        GF<-CCC(predict(model), model$data$resp)
        plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, "gof" = GF)
      }else{
        LCC.data<-list()
        CCC.data<-list()
        LA<-list()
        for(i in 1:ldb) {
          LA[[i]]<-CCC[[i]]/Pearson[[i]]
          LCC.data[[i]]<-data.frame("Time"=tk.plot,"LCC"=rho[[i]],
                                    "LPC"=rho.pearson[[i]],
                                    "LA"=Cb[[i]])
          names(LCC.data[[i]])<-
            c("Time",
              paste0(expression("LCC: "),levels(model$data$method)[i+1],
                     " vs. ", levels(model$data$method)[1]),
              paste0(expression("LPC: "),levels(model$data$method)[i+1],
                     " vs. ", levels(model$data$method)[1]),
              paste0(expression("LA: "),levels(model$data$method)[i+1],
                     " vs. ", levels(model$data$method)[1]))
          CCC.data[[i]]<-data.frame(tk.plot2, CCC[[i]], Pearson[[i]],
                                    LA[[i]])
          names(CCC.data[[i]])<-
            c("Time", paste0(expression("CCC: "),
                             levels(model$data$method)[i+1],
                             " vs. ", levels(model$data$method)[1]),
              paste0(expression("Pearson: "),
                     levels(model$data$method)[i+1], " vs. ",
                     levels(model$data$method)[1]),
              paste0(expression("Cb: "),levels(model$data$method)[i+1],
                     " vs. ", levels(model$data$method)[1]))
        }
        GF<-CCC(predict(model), model$data$resp)
        plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, "gof" = GF)
      }
    }else{
      if(ldb==1){
        CCC<-CCC_lin(dataset=model$data, resp="resp", subject="subject",
                     method="method", time="time")
        Pearson<-Pearson(dataset=model$data, resp="resp",
                         subject="subject", method="method", time="time")
        LA<-CCC[[1]]/Pearson[[1]]
        LCC.data<-data.frame("Time"=tk.plot,"LCC"=rho,
                             "Lower"=ENV.LCC[1,], "Upper"=ENV.LCC[2,])
        LPC.data<-data.frame("Time"=tk.plot,"LCC"=rho.pearson,
                             "Lower"=ENV.LPC[1,], "Upper"=ENV.LPC[2,])
        LA.data<-data.frame("Time"=tk.plot,"LA"=Cb, "Lower"=ENV.Cb[1,],
                            "Upper"=ENV.Cb[2,])
        names(LCC.data)<-
          c("Time",
            paste0(expression("LCC: "),levels(model$data$method)[2],
                   " vs. ", levels(model$data$method)[1]),
            "Lower", "Upper")
        names(LPC.data)<-
          c("Time", paste0(expression("LPC: "),
                           levels(model$data$method)[2],
                           " vs. ", levels(model$data$method)[1]),
            "Lower", "Upper")
        names(LA.data)<-
          c("Time", paste0(expression("LA: "),
                           levels(model$data$method)[2],
                           " vs. ", levels(model$data$method)[1]),
            "Lower", "Upper")
        CCC.data<-data.frame(tk.plot2, CCC, Pearson, LA)
        names(CCC.data)<-
          c("Time",
            paste0(expression("CCC: "), levels(model$data$method)[2],
                   " vs. ", levels(model$data$method)[1]),
            paste0(expression("Pearson: "), levels(model$data$method)[2],
                   " vs. ", levels(model$data$method)[1]),
            paste0(expression("Cb: "), levels(model$data$method)[2],
                   " vs. ", levels(model$data$method)[1]))
        fit<-list("LCC" = LCC.data, "LPC" = LPC.data, "LA" = LA.data)
        GF<-CCC(predict(model), model$data$resp)
        plot.data<-list("fitted"=fit,"sampled" = CCC.data, "gof" = GF)
      }else{
<<<<<<< HEAD
        CCC<-CCC_lin(dataset=model$data, resp="resp", subject="subject",
                     method="method", time="time")
        Pearson<-Pearson(dataset=model$data, resp="resp",
                         subject="subject", method="method",
                         time="time")
=======
        CCC<-CCC_lin(dataset=model$data, resp="resp", subject="subject", method="method", time="time")
        Pearson<-Pearson(dataset=model$data, resp="resp", subject="subject", method="method", time="time")
>>>>>>> 0defa447a19649cc6f57a67efad6946bd7b187aa
        LA<-list()
        CCC.data<-list()
        LCC.data<-list()
        LPC.data<-list()
        LA.data<-list()
        for(i in 1:ldb) {
          LA[[i]]<-CCC[[i]]/Pearson[[i]]
          LCC.data[[i]]<-data.frame("Time"=tk.plot,"LCC"=rho[[i]],
                                    "Lower"=ENV.LCC[[i]][1,],
                                    "Upper"=ENV.LCC[[i]][2,])
          names(LCC.data[[i]])<-
            c("Time",
              paste0(expression("LCC: "),
                     levels(model$data$method)[i+1], " vs. ",
                     levels(model$data$method)[1]), "Lower", "Upper")
          LPC.data[[i]]<-data.frame("Time"=tk.plot,
                                    "LPC"=rho.pearson[[i]],
                                    "Lower"=ENV.LPC[[i]][1,],
                                    "Upper"=ENV.LPC[[i]][2,])
          names(LPC.data[[i]])<-
            c("Time", paste0(expression("LPC: "),
                             levels(model$data$method)[i+1], " vs. ",
                             levels(model$data$method)[1]),
              "Lower", "Upper")
          LA.data[[i]]<-data.frame("Time"=tk.plot,"LA"=Cb[[i]],
                                   "Lower"=ENV.Cb[[i]][1,],
                                   "Upper"=ENV.Cb[[i]][2,])
          names(LA.data[[i]])<-
            c("Time", paste0(expression("LA: "),
                             levels(model$data$method)[i+1], " vs. ",
                             levels(model$data$method)[1]),
              "Lower", "Upper")
          CCC.data[[i]]<-data.frame(tk.plot2, CCC[[i]],
                                    Pearson[[i]], LA[[i]])
          names(CCC.data[[i]])<-
            c("Time", paste0(expression("CCC: "),
                             levels(model$data$method)[i+1], " vs. ",
                             levels(model$data$method)[1]),
              paste0(expression("Pearson: "),
                     levels(model$data$method)[i+1], " vs. ",
                     levels(model$data$method)[1]),
              paste0(expression("Cb: "), levels(model$data$method)[i+1],
                     " vs. ", levels(model$data$method)[1]))
        }
        fit<-list("LCC" = LCC.data, "LPC" = LPC.data, "LA" = LA.data)
        GF<-CCC(predict(model), model$data$resp)
        plot.data<-list("fitted"=fit, "sampled" = CCC.data, "gof" = GF)
      }
    }
  }
  return(invisible(plot.data))
}

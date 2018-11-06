#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: lcc_internal.R                                                #
# Contains: lccInternal function                                      #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 18/06/2018                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################
##' @title Internal function to prepare \code{lcc} objects
##'
##' @description This is an internally called function used to prepare \code{lcc} objects for calculate the longitudinal concordance correlation, longitudinal Pearson correlation, longitudinal bias corrector, and plotting
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@usp.br} and Rafael de Andrade Moral, \email{rafael_moral@@yahoo.com.br}
##'
##' @keywords internal
lccInternal <- function(model, q_f, q_r, tk, covar = covar, pdmat, diffbeta,
                         time_lcc, ci, percentileMet, alpha, nboot, labels,
                         var.class, weights.form, show.warnings,
                         components, lme.control, method.init) {
  tk2<-tk
  tk.plot2 <- tk2
  tk.plot <- tk
  if(is.null(time_lcc)==FALSE){
    if(is.null(time_lcc$time)){
    tk.plot <- time_lcc(time=model$data$time, from=time_lcc$from, to=time_lcc$to, n=time_lcc$n)
    }else{
    tk.plot <- time_lcc(time=time_lcc$time, from=time_lcc$from, to=time_lcc$to, n=time_lcc$n)
    }
  }
  ldb <- length(diffbeta)
  nd<-length(summary(model)$modelStruct$varStruct)
  if(nd<=1){
    if(ci==FALSE) {
      if(ldb == 1) {
        rho <- lccWrapper(model = model, q_f = q_f, n.delta = 1,
                           tk = tk.plot, diffbeta = as.numeric(diffbeta[[1]]))
          if(components==TRUE){
          rho.pearson<-lpcWrapper(model = model, q_f = q_f, tk = tk.plot, n.delta = 1)
          Cb<-laWrapper(model = model, q_f = q_f, n.delta = 1,
                         tk = tk.plot, diffbeta = as.numeric(diffbeta[[1]]))
        }
        summary.lcc<-lccSummary(q_f = q_f, diffbeta = diffbeta, tk = tk, 
                                 tk.plot = tk.plot, tk.plot2 = tk.plot2,
                                 rho = rho, rho.pearson = rho.pearson, Cb=Cb,
                                 model = model, ldb=ldb, ci = FALSE,
                                 components = components)
      } else {
        rho <- list()
        for(i in 1:ldb)  rho[[i]] <- lccWrapper(model = model, q_f = q_f, n.delta = 1,
                                                 tk = tk.plot, diffbeta = as.numeric(diffbeta[[i]]))
        rho.ret <- data.frame(do.call(cbind.data.frame, rho))
        if(components==TRUE){
          rho.pearson <- list()
          Cb<-list()
          for(i in 1:ldb)  rho.pearson[[i]] <- lpcWrapper(model = model, q_f = q_f, n.delta = 1,
                                                           tk = tk.plot)
          for(i in 1:ldb) Cb[[i]]<- laWrapper(model = model, q_f = q_f, n.delta = 1,
                                                tk = tk.plot, diffbeta = as.numeric(diffbeta[[i]]))
          rho.pearson.ret <- data.frame(do.call(cbind.data.frame, rho.pearson))
          Cb.ret <- data.frame(do.call(cbind.data.frame, Cb))
        }
        summary.lcc<-lccSummary(q_f = q_f, diffbeta = diffbeta, tk = tk,
                                 tk.plot = tk.plot, tk.plot2 = tk.plot2,
                                 rho = rho.ret, rho.pearson = rho.pearson.ret,
                                 Cb = Cb.ret, model = model, ldb = ldb, ci = FALSE,
                                 components = components)
      }
    }else{
   CI<-ciBuilder(model = model, nboot = nboot, q_f = q_f, q_r = q_r, covar = covar,
                  pdmat = pdmat, var.class = var.class, weights.form = weights.form,
                  show.warnings = show.warnings, tk = tk.plot, diffbeta = diffbeta,
                  ldb = ldb, tk.plot=tk.plot, tk.plot2 = tk.plot2,
                  ci=TRUE, percentileMet = percentileMet, alpha = alpha, components = components,
                  lme.control = lme.control, method.init = method.init)
   ENV.LCC<-CI$ENV.LCC
   ENV.LPC<-CI$ENV.LPC
   ENV.Cb<-CI$ENV.Cb
   summary.lcc<-lccSummary(q_f = q_f, diffbeta = diffbeta, tk = tk, tk.plot = tk.plot, tk.plot2 = tk.plot2,
                            rho = CI$rho, rho.pearson = CI$LPC, Cb = CI$Cb, ldb = ldb,
                            model = model, ENV.LCC = CI$ENV.LCC, ENV.LPC = CI$ENV.LPC,
                            ENV.Cb = CI$ENV.Cb, ci = TRUE, components = components)
    }
  }else{
    if(ci==FALSE) {
      if(ldb == 1) {
        rho <- lccWrapper(model = model, q_f = q_f, n.delta = 1,
                           tk = tk.plot, diffbeta = as.numeric(diffbeta[[1]]))
        if(components==TRUE){
          rho.pearson <- lpcWrapper(model = model, q_f = q_f, n.delta = 1,
                             tk = tk.plot)
          Cb<-laWrapper(model = model, q_f = q_f, n.delta = 1,
                          tk = tk.plot, diffbeta = as.numeric(diffbeta[[1]]))
        }
        summary.lcc<-lccSummary(q_f = q_f, diffbeta = diffbeta, tk = tk, 
                                 tk.plot = tk.plot, tk.plot2 = tk.plot2,
                                 rho = rho, rho.pearson = rho.pearson, Cb = Cb,
                                 model = model, ldb = ldb,
                                 ci = FALSE, components = components)
      } else {
        rho <- list()
        for(i in 1:ldb)  rho[[i]] <- lccWrapper(model = model, q_f = q_f, n.delta = i,
                                                 tk = tk.plot, diffbeta = as.numeric(diffbeta[[i]]))
        rho.ret <- data.frame(do.call(cbind.data.frame, rho))
        if(components == TRUE){
          rho.pearson <- list()
          Cb<- list()
          for(i in 1:ldb)  rho.pearson[[i]] <- lpcWrapper(model = model, q_f = q_f, n.delta = i,
                                                   tk = tk.plot)
          for(i in 1:ldb) Cb[[i]]<-laWrapper(model = model, q_f = q_f, n.delta = i,
                                               tk = tk.plot, diffbeta = as.numeric(diffbeta[[i]]))
          rho.pearson.ret <- data.frame(do.call(cbind.data.frame, rho.pearson))
          Cb.ret <- data.frame(do.call(cbind.data.frame, Cb))
        }
        summary.lcc<-lccSummary(q_f = q_f, diffbeta = diffbeta, tk = tk,
                                 tk.plot = tk.plot, tk.plot2 = tk.plot2,
                                 rho = rho.ret, rho.pearson = rho.pearson.ret,
                                 Cb = Cb.ret, model = model, ldb = ldb,
                                 ci = FALSE, components = components)
      }
    }else{
      CI<-ciBuilder(model = model, nboot = nboot, q_f = q_f, q_r = q_r, covar = covar,
                     pdmat = pdmat, var.class = var.class, weights.form = weights.form,
                     show.warnings = show.warnings, tk = tk.plot, diffbeta = diffbeta,
                     ldb = ldb, tk.plot=tk.plot, tk.plot2 = tk.plot2,
                     ci=TRUE, percentileMet = percentileMet, alpha =alpha, components = components,
                     lme.control = lme.control, method.init = method.init)
      ENV.LCC<-CI$ENV.LCC
      ENV.LPC<-CI$ENV.LPC
      ENV.Cb<-CI$ENV.Cb
      summary.lcc<-lccSummary(q_f = q_f, diffbeta = diffbeta, tk = tk,
                               tk.plot = tk.plot, tk.plot2 = tk.plot2,
                               rho = CI$rho, rho.pearson = CI$LPC, Cb = CI$Cb,
                               ldb = ldb, model = model, ENV.LCC = CI$ENV.LCC,
                               ENV.LPC = CI$ENV.LPC,  ENV.Cb = CI$ENV.Cb,
                               ci = TRUE, components = components)
   }
  }
  
  if(ldb == 1){
    if(ci==FALSE){
      if(components == FALSE){
  internal_lcc<-list("Summary.lcc"=summary.lcc, "rho"=rho, "tk.plot" = tk.plot, 
    "tk.plot2" = tk.plot2, "ldb" = ldb, "ci" = ci, "components" = components, "nd"=nd)
      }else{
        internal_lcc<-list("Summary.lcc"=summary.lcc, "rho"=rho, "rho.pearson" = rho.pearson, 
          "Cb"=Cb, "tk.plot" = tk.plot, "tk.plot2" = tk.plot2, "ldb" = ldb, "ci" = ci, 
          "components" = components, "nd"=nd)
  }
  }else{
    if(components == FALSE){
    internal_lcc<-list("Summary.lcc"=summary.lcc, "rho"=CI$rho, "ENV.LCC"=CI$ENV.LCC, 
      "ENV.LPC" = CI$ENV.LPC, "ENV.LA" = CI$ENV.Cb, "tk.plot" = tk.plot, "tk.plot2" = tk.plot2, 
      "ldb" = ldb, "ci" = ci, "components" = components, "nd"=nd)
    }else{
      internal_lcc<-list("Summary.lcc"=summary.lcc, "rho"=CI$rho, "rho.pearson" = CI$LPC, 
        "Cb"=CI$Cb, "ENV.LCC"=CI$ENV.LCC, "ENV.LPC" = CI$ENV.LPC, "ENV.LA" = CI$ENV.Cb, 
        "tk.plot" = tk.plot, "tk.plot2" = tk.plot2, "ldb" = ldb, "ci" = ci,
        "components" = components, "nd"=nd)  
      }
    }
  }else{
    if(ci==FALSE){
      if(components == FALSE){
    internal_lcc<-list("Summary.lcc"=summary.lcc, "rho"=rho.ret, "tk.plot" = tk.plot, 
      "tk.plot2" = tk.plot2, "ldb" = ldb, "ci" = ci, "components" = components, "nd"=nd)
      }else{
        internal_lcc<-list("Summary.lcc"=summary.lcc, "rho"=rho.ret, "rho.pearson" = rho.pearson.ret, 
          "Cb"=Cb.ret, "tk.plot" = tk.plot,"tk.plot2" = tk.plot2, "ldb" = ldb, 
          "ci" = ci, "components" = components, "nd"=nd)
    }
  }else{
    if(components == FALSE){
    internal_lcc<-list("Summary.lcc"=summary.lcc, "rho"=CI$rho, "ENV.LCC"=CI$ENV.LCC, 
      "ENV.LPC" = CI$ENV.LPC, "ENV.LA" = CI$ENV.Cb, "tk.plot" = tk.plot, 
      "tk.plot2" = tk.plot2, "ldb" = ldb, "ci" = ci, "components" = components, "nd"=nd)
    }else{
      internal_lcc<-list("Summary.lcc"=summary.lcc, "rho"=CI$rho, "rho.pearson" = CI$LPC, 
        "Cb"=CI$Cb, "ENV.LCC"=CI$ENV.LCC, "ENV.LPC" = CI$ENV.LPC, "ENV.LA" = CI$ENV.Cb, 
        "tk.plot" = tk.plot, "tk.plot2" = tk.plot2, "ldb" = ldb, "ci" = ci, 
        "components" = components, "nd"=nd)
      }
    }
  }
  
  return(invisible(internal_lcc))
}

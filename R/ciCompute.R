#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: ciCompute.R                                                   #
# Contains: ciCompute function                                        #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 18/06/2018                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal function to compute the non-parametric bootstrap interval.
##'
##' @description This is an internally called function used to compute the non-parametric bootstrap interval.
##'
##' @usage NULL
##'
##' @details returns a matrix or list of matrix containing the non-parametric bootstrap interval.
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@usp.br}
##'
##' @importFrom stats quantile sd qnorm
##'
##' @keywords internal
ciCompute<-function(rho, rho.pearson, Cb, tk.plot, tk.plot2, ldb, model, ci,
                  percentileMet, LCC_Boot, LPC_Boot, Cb_Boot, alpha){
  ZFisher<-function(x){
    1/2*log((1+x)/(1-x))
  }
  Logit<-function(x){
    log(x/(1-x))
  }
    if(ldb == 1) {
    LCC_IC <- matrix(0, ncol=length(LCC_Boot), nrow=length(LCC_Boot[[1]]))
    if(percentileMet=="TRUE"){
      for(i in 1:length(LCC_Boot)) {
        if(is.null(LCC_Boot[[i]])==FALSE){
          LCC_IC[,i] <- LCC_Boot[[i]]
        }else(cat(i,"\n"))
      }
    ENV.LCC <- apply(LCC_IC, 1, quantile, probs=c(alpha/2,1-alpha/2))
    } else{
    for(i in 1:length(LCC_Boot)) {
      if(is.null(LCC_Boot[[i]])==FALSE){
        LCC_IC[,i] <- ZFisher(LCC_Boot[[i]])
      }else(cat(i,"\n"))
     }
    SE<-apply(LCC_IC, 1, sd)
    mean<-apply(LCC_IC, 1, mean)
    ENV.LCC<-matrix(NA, nrow = 2, ncol = length(SE))
    for(i in 1:length(SE)){
    ENV.LCC[,i]<-c(mean[i], mean[i])-c(qnorm(1-alpha/2)*SE[i],qnorm(alpha/2)*SE[i])
     }
    ENV.LCC<-(exp(2*ENV.LCC)-1)/(exp(2*ENV.LCC)+1)
    }
    LPC_IC <- matrix(0, ncol=length(LPC_Boot), nrow=length(LPC_Boot[[1]]))
    if(percentileMet=="TRUE"){
    for(i in 1:length(LPC_Boot)) {
      if(is.null(LPC_Boot[[i]])==FALSE){
        LPC_IC[,i] <- LPC_Boot[[i]]
      }else(cat(i,"\n"))
    }
    ENV.LPC <- apply(LPC_IC, 1, quantile, probs=c(alpha/2,1-alpha/2))
    } else{
      for(i in 1:length(LPC_Boot)) {
        if(is.null(LPC_Boot[[i]])==FALSE){
          LPC_IC[,i] <- ZFisher(LPC_Boot[[i]])
        }else(cat(i,"\n"))
      }
      SE<-apply(LPC_IC, 1, sd)
      mean<-apply(LPC_IC, 1, mean)
      ENV.LPC<-matrix(NA, nrow = 2, ncol = length(SE))
      for(i in 1:length(SE)){
        ENV.LPC[,i]<-c(mean[i], mean[i])-c(qnorm(1-alpha/2)*SE[i],qnorm(alpha/2)*SE[i])
      }
      ENV.LPC<-(exp(2*ENV.LPC)-1)/(exp(2*ENV.LPC)+1)
    }
    Cb_IC <- matrix(0, ncol=length(Cb_Boot), nrow=length(Cb_Boot[[1]]))
    if(percentileMet=="TRUE"){
    for(i in 1:length(Cb_Boot)) {
      if(is.null(Cb_Boot[[i]])==FALSE){
        Cb_IC[,i] <- Cb_Boot[[i]]
      }else(cat(i,"\n"))
    }
    ENV.Cb <- apply(Cb_IC, 1, quantile, probs=c(alpha/2,1-alpha/2))
    } else{
      for(i in 1:length(Cb_Boot)) {
        if(is.null(Cb_Boot[[i]])==FALSE){
          Cb_IC[,i] <- Logit(Cb_Boot[[i]])
        }else(cat(i,"\n"))
      }
      SE<-apply(Cb_IC, 1, sd)
      mean<-apply(Cb_IC, 1, mean)
      ENV.Cb<-matrix(NA, nrow = 2, ncol = length(SE))
      for(i in 1:length(SE)){
        ENV.Cb[,i]<-c(mean[i], mean[i])-c(qnorm(1-alpha/2)*SE[i],qnorm(alpha/2)*SE[i])
      }
      ENV.Cb<-exp(ENV.Cb)/(1+exp(ENV.Cb))
    }
    CI.LCC<-list("rho"=rho,"ENV.LCC"=ENV.LCC,"LPC"=rho.pearson,"ENV.LPC"=ENV.LPC,
                 "Cb"= Cb, "ENV.Cb" = ENV.Cb)
  }else{
    LCC_IC<-list(NA)
    ENV.LCC<-list(NA)
    SE_LCC<-list()
    mean_LCC<-list()
    for(i in 1:ldb){
      LCC_IC[[i]] <- matrix(0, ncol=length(LCC_Boot), nrow=length(LCC_Boot[[1]][[i]]))
      if(percentileMet=="TRUE"){
      for(j in 1:length(LCC_Boot)) {
        if(is.null(LCC_Boot[[j]])==FALSE){
          LCC_IC[[i]][,j] <- LCC_Boot[[j]][[i]]
          }else(cat(i,"\n"))
        }
      ENV.LCC[[i]] <- apply(LCC_IC[[i]], 1, quantile, probs=c(alpha/2,1-alpha/2))
      } else{
        for(j in 1:length(LCC_Boot)) {
          if(is.null(LCC_Boot[[j]])==FALSE){
            LCC_IC[[i]][,j] <- ZFisher(LCC_Boot[[j]][[i]])
          }else(cat(i,"\n"))
        }
        SE_LCC[[i]]<-apply(LCC_IC[[i]], 1, sd)
        mean_LCC[[i]]<-apply(LCC_IC[[i]], 1, mean)
        ENV.LCC[[i]]<-matrix(NA, nrow = 2, ncol = length(SE_LCC[[i]]))
        for(k in 1:length(SE_LCC[[i]])){
          ENV.LCC[[i]][,k]<-c(mean_LCC[[i]][k], mean_LCC[[i]][k])-c(qnorm(1-alpha/2)*SE_LCC[[i]][k],qnorm(alpha/2)*SE_LCC[[i]][k])
        }
        ENV.LCC[[i]]<-(exp(2*ENV.LCC[[i]])-1)/(exp(2*ENV.LCC[[i]])+1)
      }
    }
    LPC_IC<-list(NA)
    ENV.LPC<-list(NA)
    SE_LPC<-list()
    mean_LPC<-list()
    for(i in 1:ldb){
      LPC_IC[[i]] <- matrix(0, ncol=length(LPC_Boot), nrow=length(LPC_Boot[[1]][[i]]))
      if(percentileMet=="TRUE"){
      for(j in 1:length(LPC_Boot)) {
        if(is.null(LPC_Boot[[j]])==FALSE){
          LPC_IC[[i]][,j] <- LPC_Boot[[j]][[i]]
        }else(cat(i,"\n"))
      }
      ENV.LPC[[i]] <- apply(LPC_IC[[i]], 1, quantile, probs=c(alpha/2,1-alpha/2))
      } else{
        for(j in 1:length(LPC_Boot)) {
          if(is.null(LPC_Boot[[j]])==FALSE){
            LPC_IC[[i]][,j] <- ZFisher(LPC_Boot[[j]][[i]])
          }else(cat(i,"\n"))
        }
        SE_LPC[[i]]<-apply(LPC_IC[[i]], 1, sd)
        mean_LPC[[i]]<-apply(LPC_IC[[i]], 1, mean)
        ENV.LPC[[i]]<-matrix(NA, nrow = 2, ncol = length(SE_LPC[[i]]))
        for(k in 1:length(SE_LPC[[i]])){
          ENV.LPC[[i]][,k]<-c(mean_LPC[[i]][k], mean_LPC[[i]][k])-c(qnorm(1-alpha/2)*SE_LPC[[i]][k],qnorm(alpha/2)*SE_LPC[[i]][k])
        }
        ENV.LPC[[i]]<-(exp(2*ENV.LPC[[i]])-1)/(exp(2*ENV.LPC[[i]])+1)
      }
    }
    Cb_IC<-list(NA)
    ENV.Cb<-list(NA)
    SE_Cb<-list()
    mean_Cb<-list()
    for(i in 1:ldb){
      Cb_IC[[i]] <- matrix(0, ncol=length(Cb_Boot), nrow=length(Cb_Boot[[1]][[i]]))
      if(percentileMet=="TRUE"){
      for(j in 1:length(Cb_Boot)) {
        if(is.null(Cb_Boot[[j]])==FALSE){
          Cb_IC[[i]][,j] <- Cb_Boot[[j]][[i]]
        }else(cat(i,"\n"))
      }
      ENV.Cb[[i]] <- apply(Cb_IC[[i]], 1, quantile, probs=c(alpha/2,1-alpha/2))
      } else {
        for(j in 1:length(Cb_Boot)) {
          if(is.null(Cb_Boot[[j]])==FALSE){
            Cb_IC[[i]][,j] <- Logit(Cb_Boot[[j]][[i]])
          }else(cat(i,"\n"))
        }
        SE_Cb[[i]]<-apply(Cb_IC[[i]], 1, sd)
        mean_Cb[[i]]<-apply(Cb_IC[[i]], 1, mean)
        ENV.Cb[[i]]<-matrix(NA, nrow = 2, ncol = length(SE_Cb[[i]]))
        for(k in 1:length(SE_Cb[[i]])){
          ENV.Cb[[i]][,k]<-c(mean_Cb[[i]][k], mean_Cb[[i]][k])-c(qnorm(1-alpha/2)*SE_Cb[[i]][k],qnorm(alpha/2)*SE_Cb[[i]][k])
        }
        ENV.Cb[[i]]<-exp(ENV.Cb[[i]])/(1+exp(ENV.Cb[[i]]))
      }
    }
    CI.LCC<-list("rho"=rho,"ENV.LCC"=ENV.LCC,
                 "LPC"=rho.pearson,"ENV.LPC"=ENV.LPC,
                 "Cb"=Cb,"ENV.Cb"=ENV.Cb)
  }
  return(CI.LCC)
}

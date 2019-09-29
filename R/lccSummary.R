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

##' @rdname summary.lcc
##' @method summary lcc
##' @title  Summarize an lcc object
##'
##' @description Additional information about the fit of longitudinal
##'   concordance correlation, longitudinal Pearson correlation, and
##'   longitudinal accuracy represented by an object of class
##'   \code{\link[lcc]{lcc}}. The returned object has a
##'   \code{\link[base]{print}} method.
##'
##' @return an object inheriting from class \code{summary.lcc}
##'   including: \item{fitted}{the fitted values extracted from the
##'   \code{lcc} object.} \item{gof}{the goodness of fit measurement 
##'   is calculated using the concordance correlation coefficient between 
##'   fitted and observed values. Value of 1 denote perfect concordance.}
##'   \item{AIC}{the Akaike Information Criterion corresponding to object.}
##'   \item{BIC}{the Bayesian Information Criterion corresponding to object.}
##'   \item{logLik}{If \code{REML=FALSE}, returns the log-likelihood value 
##'   of the linear mixed-effects model; otherwise, the restricted 
##'   log-likelihood is returned}
##'
##' @param object an object inheriting from class
##'   \code{\link[lcc]{lcc}}, representing a fitted longitudinal
##'   concordance correlation function.
##'
##' @param type an optional character string specifying the type of
##'   output to be returned. If \code{type="model"}, prints the summary
##'   of the polynomial mixed-effects regression model. If
##'   \code{type="lcc"}, prints the summary of the fitted and sampled
##'   values for LCC, LPC, and LA as well as the concordance correlation
##'   coefficient between fitted LCC values and observed values as
##'   goodness of fit (gof) measurement. Defaults to \code{type="lcc"}.
##'
##' @param ... not used.
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@usp.br}, Rafael de Andrade Moral
##'
##' @references Lin, L. A Concordance Correlation Coefficient to
##'   Evaluate Reproducibility. \emph{Biometrics}, 45, n. 1, 255-268,
##'   1989.
##' @references Oliveira, T.P.; Hinde, J.; Zocchi S.S. Longitudinal
##'   Concordance Correlation Function Based on Variance Components: An
##'   Application in Fruit Color Analysis. \emph{Journal of
##'   Agricultural, Biological, and Environmental Statistics}, v. 23,
##'   n. 2, 233–254, 2018.
##'
##' @seealso \code{\link[lcc]{lcc}}.
##'
##' @examples
##'
##' data(hue)
##' ## Second degree polynomial model with random intercept, slope and
##' ## quadratic term
##' fm1<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean",
##'          method = "Method", time = "Time", qf = 2, qr = 2)
##' summary(fm1, type="model")
##' summary(fm1, type="lcc")
##' @export
summary.lcc <- function(object, type, ...){
  if(class(object)!="lcc") stop("Object must inherit from class \"lcc\"",
                             call.=FALSE)
  ret <- list()
  if(missing(type)) type="lcc"
  if(type=="model" | type=="lcc"){
    cat("Longitudinal concordance correlation model fit by ")
    cat( if(object[1]$model$method == "REML") "REML\n" else "maximum likelihood\n")
    ret$AIC <- AIC(object[1]$model)
    ret$BIC <- BIC(object[1]$model)
    ret$logLik <- c(object[1]$model$logLik)
    print(with(ret, data.frame(AIC, BIC, logLik, row.names = " ")), ...)
    cat("\n")
    if(type == "lcc"){
      if(class(object$Summary.lcc$comp) == "character"){
        if(is.null(object$plot_info$ENV.LCC)){
          cat(paste0(" gof: ", round(object$Summary.lcc$gof, 4)), "\n")
          cat("\n")
          cat(object$Summary.lcc$comp, "\n")
          ret$fitted <- object$Summary.lcc$fitted
          print(ret$fitted, ...)  
        }else{
          cat(paste0(" Lower and upper bound of ", (1-object$plot_info$alpha)*100,"%"), "bootstrap confidence interval", "\n")
          cat(" Number of bootstrap samples: ", object$plot_info$nboot, "\n")
          cat(paste0(" gof: ", round(object$Summary.lcc$gof, 4)), "\n")
          cat("\n")
          cat(object$Summary.lcc$comp, "\n")
          ret$fitted <- object$Summary.lcc$fitted
          print(ret$fitted, ...) 
        }
      }else{
        summ <- sum(sapply(object$Summary.lcc$comp, length))
        if(is.null(object$plot_info$ENV.LCC)){
          cat(paste0(" gof: ", round(object$Summary.lcc$gof, 4)), "\n")
          cat("\n")
          for(i in 1:summ){
            cat(object$Summary.lcc$comp[[i]], "\n")
            ret$fitted <- object$Summary.lcc$fitted
            print(ret$fitted[[i]]) 
            cat("\n")
          }
        }else{
          cat(paste0(" Lower and upper bound of ", (1-object$plot_info$alpha)*100,"%"), "bootstrap confidence interval", "\n")
          cat(" Number of bootstrap samples: ", object$plot_info$nboot, "\n")
          cat(paste0(" gof: ", round(object$Summary.lcc$gof, 4)), "\n")
          cat("\n")
          for(i in 1:summ){
            cat(object$Summary.lcc$comp[[i]], ": LCC", "\n")
            ret$fitted <- object$Summary.lcc$fitted
            print(ret$fitted$LCC[[i]])   
            cat("\n")
            cat(object$Summary.lcc$comp[[i]], ": LPC", "\n")
            print(ret$fitted$LPC[[i]])
            cat("\n")
            cat(object$Summary.lcc$comp[[i]], ": LA", "\n")
            print(ret$fitted$LA[[i]])
            cat("\n", "\n")
          }
        }
      }
    }
    if(type == "model"){
      cat("  Fixed Effects:", "\n")
      ret$model <- list(modelStruct = object[1]$model$modelStruct,
                        dims = object[1]$model$dims,
                        contrasts = object[1]$model$contrasts,
                        coefficients = object[1]$model$coefficients,
                        varFix = object[1]$model$varFix,
                        sigma = object[1]$model$sigma,
                        apVar = object[1]$model$apVar,
                        numIter = object[1]$model$numIter,
                        groups = object[1]$model$groups,
                        call = object[1]$model$call,
                        terms = object[1]$model$terms,
                        method = object[1]$model$method,
                        fitted = object[1]$model$fitted,
                        residuals = object[1]$model$residuals,
                        fixDF = object[1]$model$fixDF,
                        na.action = object[1]$model$na.action,
                        data = object[1]$model$data)
      print(fixef(object[1]$model), ...)
      cat("\n")
      dd <- object[1]$model$dims
      print(summary(object[1]$model$modelStruct), sigma = object[1]$model$sigma, ...)
      cat("Number of Observations:", dd[["N"]])
      cat("\nNumber of Groups: ")
      Ngrps <- dd$ngrps[1:dd$Q]
      if ((lNgrps <- length(Ngrps)) == 1) {	# single nesting
        cat(Ngrps,"\n")
      } else {				# multiple nesting
        sNgrps <- 1:lNgrps
        aux <- rep(names(Ngrps), sNgrps)
        aux <- split(aux, array(rep(sNgrps, lNgrps),
                                c(lNgrps, lNgrps))[!lower.tri(diag(lNgrps))])
        names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
        cat("\n")
        print(rev(Ngrps), ...)
      }
    }
    class(ret) <- "summary.lcc"
    return(invisible(ret))
  }else{stop("Available 'type' are lcc or model", call.=FALSE)}
}

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
      CCC<-CCC_lin(dataset=model$data, resp="resp", subject="subject",
                   method="method", time="time")
    if(ldb==1){
        comp <- paste0(levels(model$data$method)[2], " vs. ",
                       levels(model$data$method)[1])
        LCC.data<-data.frame("Time"=tk.plot,"LCC"=rho)
        CCC.data<-data.frame("Time" = tk.plot2, "CCC" = CCC)
        colnames(CCC.data) <- c("Time", "CCC")
        GF<-CCC(predict(model), model$data$resp)
        plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data,
                        "gof" = GF, "comp"=comp)
    }else{
      LCC.data<-list()
      comp <- list()
      for(i in 1:ldb) {
        comp[[i]] <- paste0(levels(model$data$method)[i+1],
                            " vs. ", levels(model$data$method)[1]) 
        LCC.data[[i]]<-data.frame("Time"=tk.plot,"LCC"=rho[[i]])
        CCC.data<-data.frame("Time" = tk.plot2, "CCC" = CCC)
        colnames(CCC.data) <- c("Time", "CCC")
      }
      GF<-CCC(predict(model), model$data$resp)
      plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, "gof" = GF,
                      "comp"=comp)
    }
  }else{
    if(ldb==1){
      comp = paste0(levels(model$data$method)[2],
                    " vs. ", levels(model$data$method)[1])
      CCC<-CCC_lin(dataset=model$data, resp="resp",
                   subject="subject", method="method", time="time")
      LCC.data<-data.frame("Time"=tk.plot,"LCC"=rho,
                           "Lower"=ENV.LCC[1,], "Upper"=ENV.LCC[2,])
      CCC.data<-data.frame("Time" = tk.plot2, "CCC" = CCC)
      colnames(CCC.data) <- c("Time", "CCC")
      GF<-CCC(predict(model), model$data$resp)
       plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, "gof" = GF,
                       "comp" = comp)
    }else{
      CCC<-CCC_lin(dataset=model$data, resp="resp",
                   subject="subject", method="method", time="time")
      LCC.data<-list()
      comp <- list()
      for(i in 1:ldb) {
        comp[[i]] <- paste0(levels(model$data$method)[i+1],
                            " vs. ", levels(model$data$method)[1])
        LCC.data[[i]]<-data.frame("Time" = tk.plot,"LCC"=rho[[i]],
                                  "Lower" = ENV.LCC[[i]][1,],
                                  "Upper" = ENV.LCC[[i]][2,])
        CCC.data<-data.frame("Time" = tk.plot2, "CCC" = CCC)
        colnames(CCC.data) <- c("Time", "CCC")
      }
      GF<-CCC(predict(model), model$data$resp)
      plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, 
                      "gof" = GF,"comp" = comp)
      }
    }
  }else{
    if(ci==FALSE){
      CCC<-CCC_lin(dataset=model$data, resp="resp",
                   subject="subject", method="method", time="time")
      Pearson<-Pearson(dataset=model$data, resp="resp",
                       subject="subject", method="method", time="time")
       if(ldb==1){
        comp <- paste0(levels(model$data$method)[2], " vs. ",
                       levels(model$data$method)[1]) 
        LA <- CCC[[1]]/Pearson[[1]]
        LCC.data<-data.frame("Time"=tk.plot,"LCC"=rho,
                             "LPC"=rho.pearson, "LA"=Cb)
        CCC.data<-data.frame("Time" = tk.plot2, "CCC" = CCC, 
                             "Pearson" = Pearson, "Cb" = LA)
        colnames(CCC.data) <- c("Time", "CCC", "Pearson", "Cb")
        GF<-CCC(predict(model), model$data$resp)
        plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, 
                        "gof" = GF, "comp"=comp)
      }else{
        LCC.data <- list()
        CCC.data <- list()
        LA <- list()
        comp <- list()
        for(i in 1:ldb) {
          comp[[i]] <- paste0(levels(model$data$method)[i+1],
                              " vs. ", levels(model$data$method)[1])
          LA[[i]] <- CCC[[i]]/Pearson[[i]]
          LCC.data[[i]] <- data.frame("Time"=tk.plot,"LCC"=rho[[i]],
                                      "LPC"=rho.pearson[[i]],
                                      "LA"=Cb[[i]])
          CCC.data[[i]]<-data.frame("Time" = tk.plot2, "CCC" = CCC[[i]], 
                                    "Pearson" = Pearson[[i]], "Cb" = LA[[i]])
          colnames(CCC.data[[i]]) <- c("Time", "CCC", "Pearson", "Cb")
        }
        GF<-CCC(predict(model), model$data$resp)
        plot.data<-list("fitted"=LCC.data,"sampled"=CCC.data, 
                        "gof" = GF, "comp" = comp)
      }
    }else{
      if(ldb==1){
        CCC<-CCC_lin(dataset=model$data, resp="resp", subject="subject",
                     method="method", time="time")
        Pearson<-Pearson(dataset=model$data, resp="resp",
                         subject="subject", method="method", time="time")
        LA<-CCC[[1]]/Pearson[[1]]
        comp <- paste0(levels(model$data$method)[2],
                       " vs. ", levels(model$data$method)[1])
        LCC.data<-data.frame("Time"=tk.plot,"LCC"=rho,
                             "Lower"=ENV.LCC[1,], "Upper"=ENV.LCC[2,])
        LPC.data<-data.frame("Time"=tk.plot,"LPC"=rho.pearson,
                             "Lower"=ENV.LPC[1,], "Upper"=ENV.LPC[2,])
        LA.data<-data.frame("Time"=tk.plot,"LA"=Cb, "Lower"=ENV.Cb[1,],
                            "Upper"=ENV.Cb[2,])
        CCC.data<-data.frame("Time" = tk.plot2, "CCC" = CCC, 
                             "Pearson" = Pearson, "Cb" = LA)
        colnames(CCC.data) <- c("Time", "CCC", "Pearson", "Cb")
        fit<-list("LCC" = LCC.data, "LPC" = LPC.data, "LA" = LA.data)
        GF<-CCC(predict(model), model$data$resp)
        plot.data<-list("fitted"=fit,"sampled" = CCC.data, 
                        "gof" = GF, "comp" = comp)
      }else{
        CCC<-CCC_lin(dataset=model$data, resp="resp", subject="subject",
                     method="method", time="time")
        Pearson<-Pearson(dataset=model$data, resp="resp",
                         subject="subject", method="method",
                         time="time")
        LA<-list()
        CCC.data <- list()
        LCC.data <- list()
        LPC.data <- list()
        LA.data <- list()
        comp <- list()
        for(i in 1:ldb) {
          comp[[i]] <- paste0(levels(model$data$method)[i+1], " vs. ",
                              levels(model$data$method)[1])
          LA[[i]]<-CCC[[i]]/Pearson[[i]]
          LCC.data[[i]]<-data.frame("Time"=tk.plot,"LCC"=rho[[i]],
                                    "Lower"=ENV.LCC[[i]][1,],
                                    "Upper"=ENV.LCC[[i]][2,])
          LPC.data[[i]]<-data.frame("Time"=tk.plot,
                                    "LPC"=rho.pearson[[i]],
                                    "Lower"=ENV.LPC[[i]][1,],
                                    "Upper"=ENV.LPC[[i]][2,])
          LA.data[[i]]<-data.frame("Time"=tk.plot,"LA"=Cb[[i]],
                                   "Lower"=ENV.Cb[[i]][1,],
                                   "Upper"=ENV.Cb[[i]][2,])
          CCC.data[[i]]<-data.frame("Time" = tk.plot2, "CCC" = CCC[[i]],
                                    "Pearson" = Pearson[[i]], "LA" = LA[[i]])
          colnames(CCC.data[[i]]) <- c("Time", "CCC", "Pearson", "Cb")
        }
        fit<-list("LCC" = LCC.data, "LPC" = LPC.data, "LA" = LA.data)
        GF<-CCC(predict(model), model$data$resp)
        plot.data<-list("fitted"=fit, "sampled" = CCC.data, 
                        "gof" = GF, "comp" = comp)
      }
    }
  }
  return(invisible(plot.data))
}

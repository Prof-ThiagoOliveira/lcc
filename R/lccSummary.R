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

##' @rdname print.summary.lcc
##' @method print summary.lcc
##' @title  Print the summary of an lcc object
##'
##' @description Information summarizing the fitted longitudinal
##'   concordance correlation is printed. This includes the AIC, BIC,
##'   and log-likelihood at convergence. If \code{type = "lcc"}, prints
##'   the fitted values while \code{type = "model"} prints the fixed
##'   effects estimates and their standard errors, standard deviations,
##'   correlations for the random effects, within-group correlation, and
##'   variance function parameters.
##'
##' @param x an object inheriting from class
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
##' @param verbose an optional logical value used to control the amount
##'   of printed output when \code{type = "model"}. Defaults to
##'   \code{FALSE}
##'
##' @param digits a non-null value for \code{digits} specifies the
##'   minimum number of significant digits to be printed in values. The
##'   default, \code{NULL}.
##'
##' @param ... not used.
##'
##' @importFrom stats AIC BIC
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@usp.br}
##'
##' @seealso \code{\link{summary.lcc}}, \code{\link{lccPlot}},
##'   \code{\link[lcc]{lcc}}
##'
##' @examples
##'
##' data(hue)
##' ## Second degree polynomial model with random intercept, slope and
##' ## quadratic term
##' fm1<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean",
##'          method = "Method", time = "Time", qf = 2, qr = 2)
##' print(summary(fm1, type="model"))
##' @export
print.summary.lcc <- function(x, type, verbose =  FALSE, digits = NULL, ...){
  object <- x
  if(missing(type)) type="lcc"
  if(type=="model" | type=="lcc"){
    if(type == "lcc"){
      cat("Longitudinal concordance correlation model fit by ")
      cat( if(object$model$method == "REML") "REML\n" else "maximum likelihood\n")
      AIC <- AIC(object$model)
      BIC <- BIC(object$model)
      logLik <- c(object$model$logLik)
      print(data.frame(AIC, BIC, logLik, row.names = " "), digits = digits, ...)
      cat("\n")
      gof <- object$gof
      cat(paste0(" gof: ", round(gof, 4)), "\n")
      cat("\n")
      if(class(object$comp) == "character"){
        if(is.null(object$info$ENV.LCC)){
          cat(object$comp, "\n")
          fitted <- object$fitted
          print(fitted, digits = digits,  ...)
        }else{
          cat(paste0(" Lower and upper bound of ", (1-object$plot_info$alpha)*100,"%"), "bootstrap confidence interval", "\n")
          cat(" Number of bootstrap samples: ", object$plot_info$nboot, "\n")
          cat("\n")
          cat(object$comp, "\n")
          fitted <- object$fitted
          print(fitted, digits = digits,  ...)
        }
      }else{
        summ <- sum(sapply(object$comp, length))
        if(is.null(object$info$ENV.LCC)){
          for(i in 1:summ){
            cat(object$comp[[i]], "\n")
            fitted <- object$fitted
            print(fitted[[i]],  digits = digits, ...)
            cat("\n")
          }
        }else{
          cat(paste0(" Lower and upper bound of ", (1-object$info$alpha)*100,"%"), "bootstrap confidence interval", "\n")
          cat(" Number of bootstrap samples: ", object$info$nboot, "\n")
          cat("\n")
          for(i in 1:summ){
            cat(object$comp[[i]], ": LCC", "\n")
            fitted <- object$fitted
            print(fitted$LCC[[i]],  digits = digits, ...)
            cat("\n")
            cat(object$comp[[i]], ": LPC", "\n")
            print(fitted$LPC[[i]],  digits = digits, ...)
            cat("\n")
            cat(object$comp[[i]], ": LA", "\n")
            print(fitted$LA[[i]],  digits = digits, ...)
            cat("\n", "\n")
          }
        }
      }
    }
    if(type == "model"){
      print(summary(object[1]$model),  verbose = verbose)
    }
  }else{stop("Available 'type' are lcc or model", call.=FALSE)}
}

##' @rdname summary.lcc
##' @title  Summarize an lcc object
##' @usage
##' \method{summary}{lcc}(object, type, adjustSigma, verbose, ...)
##'
##' @method summary lcc
##' @aliases summary.lcc
##' @description Additional information about the fit of longitudinal
##'   concordance correlation, longitudinal Pearson correlation, and
##'   longitudinal accuracy represented by an object of class
##'   \code{\link[lcc]{lcc}}. The returned object has a
##'   \code{\link[base]{print}} method.
##'
##' @return an object inheriting from class \code{summary.lcc}
##'   including: \item{fitted}{the fitted values extracted from the
##'   \code{lcc} object.} \item{gof}{the goodness of fit (gof) measurement
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
##'   coefficient between fitted values from the model and observed
##'   values as goodness of fit (gof) measurement. Defaults to
##'   \code{type="lcc"}.
##'
##' @param adjustSigma an optional logical value used when \code{type =
##'   model}. If TRUE and the estimation method used to obtain object
##'   was maximum likelihood, the residual standard error is multiplied
##'   by \code{sqrt(nobs/(nobs - npar))}.  See
##'   \code{\link[nlme]{summary.lme}} for more information.  Default is
##'   TRUE.
##'
##' @param verbose an optional logical value used to control the amount
##' of output in the \code{print.summary.lme} method when
##' \code{type = model} is used. Defaults to FALSE.
##'
##' @param ...  additional arguments affecting the summary produced.
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@usp.br}
##'
##' @importFrom stats AIC BIC
##'
##' @seealso \code{\link{AIC}}, \code{\link{BIC}},
##' \code{print.summary.lcc},  \code{\link[lcc]{lcc}}
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
summary.lcc <- function(object, type, adjustSigma = TRUE,
                        verbose = FALSE, ...)
{
  if(class(object)!="lcc") stop("Object must inherit from class \"lcc\"",
                                call.=FALSE)
  if(missing(type)) type="lcc"
  if(type=="model" | type=="lcc"){
    if(type == "lcc"){
      # Object lcc
      object$call <- object[1]$model$call
      object$fitted <- object$Summary.lcc$fitted
      object$sampled <- object$Summary.lcc$sampled
      object$gof <- object$Summary.lcc$gof
      object$data <- object$dataset
      object$AIC <- AIC(object[1]$model)
      object$BIC <- BIC(object[1]$model)
      object$logLik <- c(object[1]$model$logLik)
      object$info <- object$plot_info
      object$comp <- object$Summary.lcc$comp
      object <- object[names(object) != "Summary.lcc"]
      object <- object[names(object) != "plot_info"]
      object <- object[names(object) != "dataset"]
      #-----------------------------------------------------------------
      ## generating the final object
      #-----------------------------------------------------------------
      structure(object, type = type,  oClass = class(object),
                class = c("summary.lcc", class(object)))
    }else {
    #-------------------------------------------------------------------
    # Model
    #-------------------------------------------------------------------
      obj <- object[1]$model
      summary(obj,  adjustSigma = TRUE, verbose = FALSE)
    }
  }else {
  stop("Available 'type' are lcc or model", call.=FALSE)
  }
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

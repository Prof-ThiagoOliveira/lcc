#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: lcc.R                                                         #
# Contains: generic methods                                           #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 30/10/2019                                           #
# Last update: 01/11/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

#=======================================================================
# is.lcc function
#=======================================================================
##' Reports whether x is a lcc object
##' @rdname is.lcc
##' @param x An object to test
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@alumni.usp.br}
##' @keywords internal
##' @export
is.lcc <- function(x) inherits(x, "lcc")

#=======================================================================
# Print method
#=======================================================================
##' @rdname print.lcc
##' @method print lcc
##' @title Print Method for \code{lcc} Objects
##'
##' @description Prints detailed information about the fitted longitudinal
##'   concordance correlation model contained in an \code{lcc} object.
##'
##' @param x An object of class \code{lcc}, representing a fitted longitudinal
##'   concordance correlation model.
##' @param digits Minimum number of significant digits to be printed in values.
##'   Default is \code{NULL}, which uses the default precision.
##' @param ... Further arguments passed to \code{print}.
##'
##' @return The function is used for its side effect of printing and returns
##'   the input \code{lcc} object invisibly.
##'
##' @seealso \code{\link[lcc]{lcc}}, \code{\link[lcc]{summary.lcc}}
##' @examples
##' \dontrun{
##' fm1 <- lcc(data = hue, subject = "Fruit", resp = "H_mean",
##'           method = "Method", time = "Time", qf = 2, qr = 2)
##' print(fm1)
##' }
##'
##' @export
print.lcc <- function(x, digits = NULL, ...){
  cat("Longitudinal concordance correlation model fit by ")
  cat( if(x[1]$model$method == "REML") "REML\n" else "maximum likelihood\n")
  AIC <- AIC(x[1]$model)
  BIC <- BIC(x[1]$model)
  logLik <- c(x[1]$model$logLik)
  print(data.frame(AIC, BIC, logLik, row.names = " "), digits = digits, ...)
  cat("\n")
  print(x$Summary.lcc$fitted, digits =  digits,  ...)
  dd <- x$model$dims
  Ngrps <- dd$ngrps[1:dd$Q]
  cat("\n")
  cat("Number of Observations:", dd[["N"]])
  cat("\nNumber of Groups: ")
  if ((lNgrps <- length(Ngrps)) == 1) {	# single nesting
    cat(Ngrps,"\n")
  } else {				# multiple nesting
    sNgrps <- 1:lNgrps
    aux <- rep(names(Ngrps), sNgrps)
    aux <- split(aux, array(rep(sNgrps, lNgrps),
                            c(lNgrps, lNgrps))[!lower.tri(diag(lNgrps))])
    names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
    cat("\n")
    print(rev(Ngrps),  digits = digits, ...)
  }
  invisible(x)
}

#=======================================================================
# fitted method
#=======================================================================
##' @rdname fitted.lcc
##' @method fitted lcc
##' @title Extract Fitted Values from an \code{lcc} Object
##'
##' @description Extracts and prints the fitted values from an object of class
##'   \code{lcc}, as returned by modeling functions. The function allows selection
##'   of different types of fitted values based on longitudinal data analysis.
##'
##' @param object An object of class \code{lcc}, representing a fitted longitudinal
##'   concordance correlation model.
##' @param type The type of fitted values to extract: "lcc" for longitudinal
##'   concordance correlation, "lpc" for longitudinal Pearson correlation,
##'   or "la" for longitudinal accuracy. Defaults to "lcc".
##' @param digits Minimum number of significant digits to be printed.
##'   Default is \code{NULL}, which uses the default precision.
##' @param ... Additional arguments (currently not used).
##'
##' @return The function prints the fitted values and returns them as a data frame.
##'
##' @seealso \code{\link[lcc]{lcc}}, \code{\link[lcc]{summary.lcc}},
##'   \code{\link[lcc]{lccPlot}}
##' @examples
##' data(hue)
##' fm1 <- lcc(data = hue, subject = "Fruit", resp = "H_mean",
##'            method = "Method", time = "Time", qf = 2, qr = 2,
##'            components = TRUE)
##' fitted(fm1)
##' fitted(fm1, type = "lpc")
##' fitted(fm1, type = "la")
##'
##' @export
fitted.lcc <- function(object, type = "lcc", digits = NULL, ...){
  if (!inherits(object, "lcc")) {
    stop("Object must inherit from class 'lcc'", call. = FALSE)
  }
  if (!type %in% c("lcc", "lpc", "la")) {
    stop("Available 'type' are 'lcc', 'lpc', or 'la'", call. = FALSE)
  }
  if (object$plot_info$components == FALSE && type != "lcc") {
    stop(paste0("It is necessary to include components = TRUE in the lcc() function ",
                "to calculate the fitted values for type '", type, "'"), call. = FALSE)
  }
  
  cat(paste0("Fitted ", switch(type,
                               "lcc" = "longitudinal concordance correlation function",
                               "lpc" = "longitudinal Pearson correlation function",
                               "la"  = "longitudinal accuracy function"), "\n\n"))
  return(fittedBuilder(object = object, type = type))
}


#=======================================================================
# Summary method
#=======================================================================
#-----------------------------------------------------------------------
# Print method for summary.lcc
#-----------------------------------------------------------------------
##' @rdname print.summary.lcc
##' @title Print Summary of an \code{lcc} Object
##'
##' @description Provides a detailed summary of a fitted longitudinal
##'   concordance correlation model, including AIC, BIC, log-likelihood, and
##'   other relevant statistics. The function supports detailed output for
##'   different types of model fits.
##'
##' @param x An object of class \code{\link[lcc]{summary.lcc}},
##'   representing a summarized longitudinal concordance correlation function.
##' @param verbose Logical value to control the amount of printed output for
##'   model details. Defaults to \code{FALSE}.
##' @param digits Specifies the minimum number of significant digits to be
##'   printed in values. Default is \code{NULL}.
##' @param ... Further arguments passed to \code{print}.
##'
##' @seealso \code{\link{summary.lcc}}, \code{\link{lccPlot}},
##'   \code{\link[lcc]{lcc}}
##' @examples
##' \dontrun{
##' fm1 <- lcc(data = hue, subject = "Fruit", resp = "H_mean",
##'            method = "Method", time = "Time", qf = 2, qr = 2)
##' print(summary(fm1, type = "model"))
##' }
##'
##' @export
print.summary.lcc <- function(x, verbose =  FALSE, digits = NULL, ...){
  if (inherits(x, "lcc")) {
    cat("Longitudinal concordance correlation model fit by ")
    cat( if(x$model$method == "REML") "REML\n" else "maximum likelihood\n")
    AIC <- AIC(x$model)
    BIC <- BIC(x$model)
    logLik <- c(x$model$logLik)
    print(data.frame(AIC, BIC, logLik, row.names = " "), digits = digits, ...)
    cat("\n")
    gof <- x$gof
    cat(paste0(" gof: ", round(gof, 4)), "\n")
    cat("\n")
    if (inherits(x$comp, "character")) {
      if (is.null(x$info$ENV.LCC)) {
        cat(x$comp, "\n")
        fitted <- x$fitted
        print(fitted, digits = digits,  ...)
      }else{
        cat(paste0(" Lower and upper bound of ", (1-x$plot_info$alpha)*100,"%"), "bootstrap confidence interval", "\n")
        cat(" Number of bootstrap samples: ", x$plot_info$nboot, "\n")
        cat("\n")
        cat(x$comp, "\n")
        fitted <- x$fitted
        print(fitted, digits = digits,  ...)
      }
    }else{
      summ <- sum(sapply(x$comp, length))
      if(is.null(x$info$ENV.LCC)){
        for(i in 1:summ){
          cat(x$comp[[i]], "\n")
          fitted <- x$fitted
          print(fitted[[i]],  digits = digits, ...)
          cat("\n")
        }
      }else{
        cat(paste0(" Lower and upper bound of ", (1-x$info$alpha)*100,"%"), "bootstrap confidence interval", "\n")
        cat(" Number of bootstrap samples: ", x$info$nboot, "\n")
        cat("\n")
        for(i in 1:summ){
          cat(x$comp[[i]], ": LCC", "\n")
          fitted <- x$fitted
          print(fitted$LCC[[i]],  digits = digits, ...)
          cat("\n")
          cat(x$comp[[i]], ": LPC", "\n")
          print(fitted$LPC[[i]],  digits = digits, ...)
          cat("\n")
          cat(x$comp[[i]], ": LA", "\n")
          print(fitted$LA[[i]],  digits = digits, ...)
          cat("\n", "\n")
        }
      }
    }
  }else{
    if (inherits(x, "model")) {
      dd <- x$dims
      verbose <- verbose || attr(x, "verbose")
      cat( "Linear mixed-effects model fit by " )
      cat( if(x$method == "REML") "REML\n" else "maximum likelihood\n")
      ##  method <- x$method
      cat(" Data:", deparse( x$call$data ), "\n")
      if (!is.null(x$call$subset)) {
        cat("  Subset:", deparse(asOneSidedFormula(x$call$subset)[[2L]]),"\n")
      }
      print(data.frame(AIC = x$AIC, BIC = x$BIC, logLik = c(x$logLik),
                       row.names = " "), ...)
      if (verbose) { cat("Convergence at iteration:",x$numIter,"\n") }
      cat("\n")
      print(summary(x$modelStruct), sigma = x$sigma,
            reEstimates = x$coef$random, verbose = verbose, ...)
      cat("Fixed effects: ")
      fixF <- x$call$fixed
      if (inherits(fixF, "formula") || is.call(fixF)) {
        cat(deparse(x$call$fixed), "\n")
      } else {
        cat(deparse(lapply(fixF, function(el) as.name(deparse(el)))), "\n")
      }
      ## fixed effects t-table and correlations
      xtTab <- as.data.frame(x$tTable)
      wchPval <- match("p-value", names(xtTab))
      for(i in names(xtTab)[-wchPval]) {
        xtTab[, i] <- format(zapsmall(xtTab[, i]))
      }
      xtTab[,wchPval] <- format(round(xtTab[,wchPval], 4))
      if (any(wchLv <- (as.double(levels(xtTab[, wchPval])) == 0))) {
        levels(xtTab[, wchPval])[wchLv] <- "<.0001"
      }
      row.names(xtTab) <- dimnames(x$tTable)[[1L]]
      print(xtTab, ...)
      if (nrow(x$tTable) > 1) {
        corr <- x$corFixed
        class(corr) <- "correlation"
        print(corr, title = " Correlation:", ...)
      }
      cat("\nStandardized Within-Group Residuals:\n")
      print(x$residuals, ...)
      cat("\nNumber of Observations:",x$dims[["N"]])
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
      invisible(x)
    }else{
      stop("Available only for classes summary.lcc or summary.lme", call.=FALSE)
    }
  }
}

#-----------------------------------------------------------------------
# summary.lcc
#-----------------------------------------------------------------------

##' @rdname summary.lcc
##' @title  Summarize an \code{lcc} Object
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
##'   \code{type="model"}.
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
##' @param ...  not used.
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @importFrom stats AIC BIC asOneSidedFormula pt resid
##'
##' @seealso \code{\link{AIC}}, \code{\link{BIC}},
##' \code{print.summary.lcc},  \code{\link[lcc]{lcc}}
##'
##' @examples
##'
##' ## Second degree polynomial model with random intercept, slope and
##' ## quadratic term
##' fm1<-lcc(data = hue, subject = "Fruit", resp = "H_mean",
##'          method = "Method", time = "Time", qf = 2, qr = 2)
##' summary(fm1, type="model")
##' summary(fm1, type="lcc")
##' @export
summary.lcc <- function(object, type, adjustSigma = TRUE,
                        verbose = FALSE, ...)
{
  if (!inherits(object, "lcc")) stop("Object must inherit from class \"lcc\"",
                                call.=FALSE)
  if(missing(type)) type <- "model"
  if(type=="model" || type=="lcc"){
    if(type == "lcc") {
      # Object lcc
      object$fitted <- object$Summary.lcc$fitted
      object$sampled <- object$Summary.lcc$sampled
      object$gof <- object$Summary.lcc$gof
      object$data <- object$data
      object$AIC <- AIC(object[1]$model)
      object$BIC <- BIC(object[1]$model)
      object$logLik <- c(object[1]$model$logLik)
      object$info <- object$plot_info
      object$comp <- object$Summary.lcc$comp
      object <- object[names(object) != "Summary.lcc"]
      object <- object[names(object) != "plot_info"]
      object <- object[names(object) != "data"]
      #-----------------------------------------------------------------
      ## generating the final object
      #-----------------------------------------------------------------
      structure(object, type = "lcc",  oClass = class(object),
                class = c("summary.lcc", type))
    }else {
      #-------------------------------------------------------------------
      # Model
      #-------------------------------------------------------------------
      obj <- object[1]$model
      ##  variance-covariance estimates for fixed effects
      fixed <- fixef(obj)
      stdFixed <- sqrt(diag(as.matrix(obj$varFix)))
      obj$corFixed <- array(t(obj$varFix/stdFixed)/stdFixed,
                            dim(obj$varFix), list(names(fixed),names(fixed)))
      if (adjustSigma && obj$method == "ML")
        stdFixed <- stdFixed *
        sqrt(obj$dims$N/(obj$dims$N - length(stdFixed)))
      ## fixed effects coefficients, std. deviations and t-ratios
      fDF <- obj$fixDF[["X"]]
      tVal <- fixed/stdFixed
      obj$tTable <- cbind(Value = fixed, Std.Error = stdFixed, DF = fDF,
                          "t-value" = tVal, "p-value" = 2 * pt(-abs(tVal), fDF))
      ## residuals
      resd <- resid(obj, type = "pearson")
      if (length(resd) > 5) {
        resd <- quantile(resd, na.rm = TRUE) # might have NAs from na.exclude
        names(resd) <- c("Min","Q1","Med","Q3","Max")
      }
      obj$residuals <- resd
      ## generating the final obj
      aux <- logLik(obj)
      obj$BIC <- BIC(aux)
      obj$AIC <- AIC(aux)
      structure(obj, verbose = verbose, oClass = class(obj),
                class = c("summary.lcc", "model", class(obj)))
    }
  }else {
    stop("Available 'type' are 'lcc' or 'model'", call. = FALSE)
  }
}


#=======================================================================
# Plot method
#=======================================================================
##' @rdname plot.lcc
##' @method plot lcc
##' @title Diagnostic Plots for an \code{lcc} Object
##'
##' @description 
##' Generates a series of diagnostic plots for evaluating the fit of a linear
##' mixed-effects model represented by an \code{lcc} object. This function
##' provides six types of plots, including residual plots, fitted value
##' comparisons, and normal Q-Q plots. Users can select specific plots or display
##' all by default.
##'
##' @usage
##' \method{plot}{lcc}(x, which = c(1L:6L),
##'      caption = list("Residuals vs Fitted",
##'                     "Residuals vs Time",
##'                     "Residuals by Subject",
##'                     "Observed values vs Fitted values",
##'                     "Normal Q-Q Plot (Conditional residuals)",
##'                     "Normal Q-Q Plot (Random effects)"),
##'      sub.caption = NULL, main = NULL,
##'      panel = if(add.smooth) panel.smooth else points,
##'      add.smooth = TRUE, ask = TRUE,
##'      id.n = 3, labels.id = names(residuals(x)),
##'      label.pos = c(4, 2), cex.id = 0.75, cex.caption = 1,
##'      cex.oma.man = 1.25, ...)
##'
##' @param x An object of class \code{\link[lcc]{lcc}}, representing a
##'   fitted longitudinal concordance correlation function.
##' @param which A numeric vector specifying which plots to display.
##'   The valid range is c(1L:6L), corresponding to the plot types.
##' @param caption Captions for the plots, provided as a vector or list
##'   of valid graphics annotations. Default captions are provided for each plot.
##' @param sub.caption A common sub-title for all plots; defaults to
##'   \code{NULL}.
##' @param main The main title for the plots, displayed above the captions.
##' @param panel Panel function to be used for adding points to the plots.
##'   Defaults to \code{panel.smooth} if \code{add.smooth} is \code{TRUE}, 
##'   otherwise \code{points}.
##' @param add.smooth Logical; indicates whether a smoother should be added
##'   to most plots. Defaults to \code{TRUE}.
##' @param ask Logical; if \code{TRUE}, prompts the user before displaying
##'   each plot in a multi-plot layout. Defaults to \code{TRUE}.
##' @param id.n Number of extreme points to label in the first three plots.
##' @param labels.id Labels for the extreme points, defaulting to observation
##'   numbers if \code{NULL}.
##' @param label.pos Positioning of labels in the left and right halves
##'   of the graph, applicable for plots 1-3.
##' @param cex.id Magnification factor for point labels.
##' @param cex.caption Size of the plot captions.
##' @param cex.oma.man Size of the overall margin annotation (applies only
##'   if \code{sub.caption} is above the figures in multi-plot layouts).
##' @param ... Additional graphical parameters passed to \code{par}.
##'
##' @details
##' The Q-Q plots use normalized residuals. Standardized residuals are pre-multiplied
##' by the inverse square-root factor of the estimated error correlation matrix,
##' while random effects are adjusted using the estimated variances from matrix G.
##' Simulation envelopes in Q-Q plots are generated using the \code{hnp} package.
##'
##' The function is partly adapted from \code{\link[stats]{plot.lm}}.
##'
##' @seealso \code{\link{lccPlot}}, \code{\link[lcc]{lcc}},
##'   \code{mtext}, \code{text}, \code{plotmath}
##'
##' @examples
##' \dontrun{
##' fm1 <- lcc(data = hue, subject = "Fruit", resp = "H_mean",
##'            method = "Method", time = "Time", qf = 2, qr = 2)
##' plot(fm1)
##' }
##'
##' @importFrom hnp hnp
##' @importFrom nlme getVarCov ranef
##' @importFrom grDevices as.graphicsAnnot dev.flush dev.hold dev.interactive devAskNewPage extendrange n2mfrow
##' @importFrom graphics abline boxplot mtext panel.smooth par plot points strheight text title
##' @importFrom stats fitted residuals qqnorm qqline
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@alumni.usp.br}
##' @export
plot.lcc <- function(x, which = c(1L:6L),
                     caption = list("Residuals vs Fitted",
                                    "Residuals vs Time",
                                    "Residuals by Subject",
                                    "Observed values vs Fitted values",
                                    "Normal Q-Q Plot (Conditional residuals)",
                                    "Normal Q-Q Plot (Random effects)"),
                     sub.caption = NULL, main = NULL,
                     panel = if(add.smooth) panel.smooth else points,
                     add.smooth = TRUE, ask = TRUE,
                     id.n = 3, labels.id = names(residuals(x)),
                     label.pos = c(4, 2), cex.id = 0.75, cex.caption = 1,
                     cex.oma.man = 1.25, ...) {
  
  if (!is.lcc(x))
    stop("use only with 'lcc' objects", call. = FALSE)
  
  if (!is.numeric(which) || any(which < 1) || any(which > 6))
    stop("'which' must be in 1:6")
  
  model <- x$model
  r <- residuals(model)
  r_norm <- residuals(model, type = "normalized")
  yh <- fitted(model)
  time <- model$data$time
  Response <- model$data$resp
  Subject <- model$data$subject
  
  one.fig <- prod(par("mfcol")) == 1
  
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  
  plotFunc <- function(i) {
    switch(i,
           {
             plotGeneric("Fitted values", "Residuals", yh, r)
           },
           {
             plotGeneric("Time", "Residuals", time, r)
           },
           {
             boxplot(r ~ Subject, ylab = "Residuals", main = main, ...)
             addTitleAndCaption(i)
           },
           {
             plot(Response ~ yh, ylab = "Observed Values", xlab = "Fitted Values", main = main, ...)
             abline(0, 1)
             addTitleAndCaption(i)
           },
           {
             qqnorm(r_norm, main = main, ...)
             qqline(r_norm)
             addTitleAndCaption(i)
           },
           {
             plotRandomEffects()
             addTitleAndCaption(i)
           }
    )
  }
  
  plotRandomEffects <- function() {
    vars <- sqrt(diag(getVarCov(model)))
    ranefs <- as.matrix(ranef(model))
    re <- ranefs %*% diag(1 / vars)
    ncol.re <- ncol(re)
    if (ncol.re == 1) {
      hnp(re, scale = TRUE, halfnormal = FALSE, print.on = TRUE, main = main, ...)
    } else {
      par(mfrow = rev(n2mfrow(ncol.re)))
      for (j in 1:ncol.re) {
        hnp(re[, j], scale = TRUE, halfnormal = FALSE, print.on = TRUE, main = main, ...)
        mtext(paste0("b", j - 1, "i"), 1, -1.5, cex = cex.caption)
      }
    }
  }
  
  plotGeneric <- function(xlab, ylab, x, y) {
    ylim <- range(y, na.rm = TRUE)
    if (id.n > 0) 
      ylim <- extendrange(ylim, f = 0.08)
    
    dev.hold()
    plot(x, y, xlab = xlab, ylab = ylab, main = main, ylim = ylim, type = "n", ...)
    panel(x, y, ...)
    abline(h = 0, lty = 3, col = "gray")
    
    if (id.n > 0 && !is.null(labels.id)) {
      show.r <- sort.list(abs(y), decreasing = TRUE)[1:id.n]
      labelPoints(x, y, show.r)
    }
    
    dev.flush()
  }
  
  labelPoints <- function(x, y, indices) {
    y.id <- y[indices]
    y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
    text(x[indices], y.id, labels.id[indices], cex = cex.id, xpd = TRUE,
         pos = label.pos, offset = 0.25)
  }
  
  addTitleAndCaption <- function(i) {
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(getCaption(i), 3, 0.25, cex = cex.caption)
  }
  
  getCaption <- function(k) {
    if (length(caption) < k) NA_character_ else as.graphicsAnnot(caption[[k]])
  }
  
  for (i in which) {
    plotFunc(i)
  }
  
  par(mfrow = c(1, 1))
  invisible()
}

#=======================================================================
# coef function
# =======================================================================
##' @rdname coef.lcc
##' @title Extract Model Coefficients
##' @usage \method{coef}{lcc}(object, ...)
##' @method coef lcc
##' @aliases coef.lcc
##'
##' @description The fixed effects estimated and corresponding random
##'   effects estimates are obtained at subject levels less or equal to
##'   i. The resulting estimates are returned as a data frame, with rows
##'   corresponding to subject levels and columns to coefficients.
##'
##' @param object an object inheriting from class \code{lcc},
##'   representing a fitted longitudinal concordance correlation
##'   function.
##' @param ... optional arguments passed to the \code{coef.lme}
##'   function.
##'
##' @details See methods for \code{\link{nlme}} objects to get more
##'   details.
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @seealso \code{\link[lcc]{lcc}}, \code{\link{summary.lcc}},
##'   \code{\link{lccPlot}}, \code{\link{vcov.lcc}}
##'
##' @importFrom stats coef
##'
##' @examples
##'
##' \dontrun{
##' fm1<-lcc(data = hue, subject = "Fruit", resp = "H_mean",
##'          method = "Method", time = "Time", qf = 2, qr = 2)
##' coef(fm1)
##' }
##'
##' @export
coef.lcc <- function(object, ...) {
  if (!is.lcc(object))
    stop("The 'object' must be of class 'lcc'.", call. = FALSE)
  
  x <- coef(object$model)
  colNames <- colnames(x)
  
  # Consolidate the substitutions 
  replacements <- list(
    "fixed" = "Fixed",
    "Poly" = "Time",
    "fmla." = "",
    "\\(time, degree = qr, raw = TRUE\\)" = "",
    "rand" = "Random",
    "poly" = "Time",
    "method" = ""
  )
  
  for (pattern in names(replacements)) {
    colNames <- gsub(pattern, replacements[[pattern]], colNames)
  }
  
  colnames(x) <- colNames
  class(x) <- c("coef.lcc", "ranef.lcc", "data.frame")
  
  x
}

#=======================================================================
# vcov function
#=======================================================================
##' @rdname vcov.lcc
##' @title Extract Variance-Covariance Matrix of the Fixed Effects for an lcc Object
##' @usage \method{vcov}{lcc}(object, ...)
##' @method vcov lcc
##' @aliases vcov.lcc
##'
##' @description 
##' Extracts the variance-covariance matrix of the fixed effects from a fitted
##' \code{lcc} model object. This function provides insights into the variability
##' and covariance structure of the fixed effects in the model.
##'
##' @param object An object of class \code{lcc}, representing a fitted 
##'   longitudinal concordance correlation model.
##' @param ... Optional arguments passed to the \code{vcov.lme}
##'   function from the \code{nlme} package.
##'
##' @details 
##' The function specifically retrieves the variance-covariance matrix associated
##' with the fixed effects of the \code{lcc} object, which is useful for understanding
##' the relationship between these effects. For more details on variance-covariance
##' matrices, refer to the methods for \code{\link{nlme}} objects.
##'
##' @seealso \code{\link{summary.lcc}}, \code{\link{lccPlot}},
##'   \code{\link[lcc]{lcc}}, \code{\link{coef.lcc}}
##'
##' @examples
##' \dontrun{
##' fm1 <- lcc(data = hue, subject = "Fruit", resp = "H_mean",
##'            method = "Method", time = "Time", qf = 2, qr = 2)
##' vcov(fm1)
##' }
##'
##' @importFrom stats vcov
##' @export
vcov.lcc <- function(object, ...) {
  if (!is.lcc(object))
    stop("The provided object is not of class 'lcc'.", call. = FALSE)
  
  covMatrix <- vcov(object$model, ...)
  
  # Streamline the replacement process
  replaceTerms <- c("fixed" = "Fixed", "Poly" = "Time", "method" = "")
  namesToChange <- names(replaceTerms)
  
  for (term in namesToChange) {
    pattern <- term
    replacement <- replaceTerms[[term]]
    rownames(covMatrix) <- gsub(pattern, replacement, rownames(covMatrix))
    colnames(covMatrix) <- gsub(pattern, replacement, colnames(covMatrix))
  }
  
  covMatrix
}

#=======================================================================
# getVarCov function
#=======================================================================
##' @rdname getVarCov.lcc
##' @title Extract Variance Components from a Fitted lcc Model
##' @usage \method{getVarCov}{lcc}(obj, type = "random.effects", ...)
##' @method getVarCov lcc
##' @aliases getVarCov.lcc
##'
##' @description 
##' Retrieves the variance-covariance matrix of the specified component from a fitted
##' \code{lcc} model object. The function can extract different types of variance-covariance
##' matrices based on the specified component type.
##'
##' @param obj An object of class \code{lcc}, representing a fitted 
##'   longitudinal concordance correlation model.
##' @param type Specifies the type of variance-covariance matrix to extract.
##'   Options are \code{"random.effects"} for random-effects variance-covariance,
##'   \code{"conditional"} for conditional variance-covariance of the responses, and
##'   \code{"marginal"} for marginal variance-covariance of the responses.
##'   Default is \code{"random.effects"}.
##' @param ... Optional arguments passed to the underlying \code{getVarCov}
##'   function from the \code{nlme} package.
##'
##' @details 
##' This function is useful for detailed inspection of the variance components
##' in different aspects of the model. For more information on the types of variance-covariance
##' matrices and their interpretations, refer to the documentation of the \code{nlme} package.
##'
##' @seealso \code{\link[lcc]{lcc}}, \code{\link{summary.lcc}},
##'   \code{\link{coef.lcc}}, \code{\link{vcov.lcc}}
##'
##' @importFrom nlme getVarCov
##'
##' @examples
##' \dontrun{
##' fm1 <- lcc(data = hue, subject = "Fruit", resp = "H_mean",
##'            method = "Method", time = "Time", qf = 2, qr = 2)
##' getVarCov(fm1)
##' }
##'
##' @export
getVarCov.lcc <- function(obj, type = "random.effects", ...) {
  if (!is.lcc(obj))
    stop("The provided object is not of class 'lcc'.", call. = FALSE)
  
  varCovMatrix <- getVarCov(obj$model, type = type, ...)
  
  # Simplify the process of renaming
  replacePatterns <- list(
    "fmla." = "",
    "\\(time, degree = qr, raw = TRUE\\)" = "",
    "rand" = "Random",
    "poly" = "Time"
  )
  
  for (pattern in names(replacePatterns)) {
    replacement <- replacePatterns[[pattern]]
    rownames(varCovMatrix) <- gsub(pattern, replacement, rownames(varCovMatrix))
    colnames(varCovMatrix) <- gsub(pattern, replacement, colnames(varCovMatrix))
  }
  
  varCovMatrix
}

#=======================================================================
# Residuals
#=======================================================================
##' @rdname residuals.lcc
##' @title Extract Residuals from a Fitted lcc Model
##' @usage \method{residuals}{lcc}(object, type = "response", ...)
##' @method residuals lcc
##' @aliases residuals.lcc
##'
##' @description 
##' Extracts residuals from the fitted longitudinal concordance correlation
##' model represented by an \code{lcc} object. Different types of residuals can
##' be obtained based on the specified type.
##'
##' @param object An object of class \code{lcc}, representing a fitted 
##'   longitudinal concordance correlation function.
##' @param type A character string specifying the type of residuals to extract.
##'   Options are \code{"response"} for residuals obtained by subtracting
##'   the fitted values from the response (default), \code{"pearson"} for
##'   "response" residuals divided by the estimated within-group standard error,
##'   and \code{"normalized"} for normalized residuals. Partial matching is
##'   used, so only the first character of the type is necessary.
##' @param ... Optional arguments passed to the \code{residuals.lme}
##'   function from the \code{nlme} package.
##'
##' @details 
##' The function provides a convenient way to examine the differences between
##' observed and predicted values in the model. Understanding these residuals
##' can be crucial for model diagnostics and validation. For more information,
##' refer to the methods for \code{\link{nlme}} objects.
##'
##' @seealso \code{\link[lcc]{lcc}}, \code{\link{summary.lcc}},
##'   \code{\link{coef.lcc}}, \code{\link{vcov.lcc}}
##'
##' @examples
##' \dontrun{
##' fm1 <- lcc(data = hue, subject = "Fruit", resp = "H_mean",
##'            method = "Method", time = "Time", qf = 2, qr = 2)
##' residuals(fm1)
##' }
##'
##' @export
residuals.lcc <- function(object, type = "response", ...) {
  if (!is.lcc(object))
    stop("The provided object is not of class 'lcc'.", call. = FALSE)
  
  residuals(object$model, type = type, ...)
}

#=======================================================================
# AIC and BIC
#=======================================================================
##' @rdname AIC.lcc
##' @title Akaike and Bayesian Information Criteria for an \code{lcc} Object
##' @method AIC lcc
##' @aliases AIC.lcc
##' @description 
##' Calculates the Akaike Information Criterion (AIC) or the Bayesian Information Criterion (BIC)
##' for a fitted longitudinal concordance correlation model represented by an \code{lcc} object.
##'
##' @param object An object inheriting from class \code{lcc},
##'   representing a fitted longitudinal concordance correlation function.
##' @param k Numeric value used as a penalty coefficient for the number of
##'   parameters in the fitted model; the default \code{k = 2} corresponds
##'   to the classical AIC.
##' @param ... Optional arguments passed to the underlying \code{AIC}
##'   function from the \code{stats} package.
##'
##' @details 
##' The function computes AIC or BIC values as a measure of the relative quality of
##' statistical models for a given set of data. Lower AIC or BIC values indicate a
##' better model fit with fewer parameters. For more information, refer to the methods
##' for \code{\link{AIC}} objects.
##'
##' @seealso \code{\link[lcc]{lcc}}, \code{\link{summary.lcc}},
##'   \code{\link{coef.lcc}}, \code{\link{vcov.lcc}}
##'
##' @importFrom stats AIC logLik na.omit
##'
##' @examples
##' \dontrun{
##' fm1 <- lcc(data = hue, subject = "Fruit", resp = "H_mean",
##'            method = "Method", time = "Time", qf = 2, qr = 2)
##' AIC(fm1)
##' }
##'
##' @export

AIC.lcc <- function(object, ..., k = 2) {
  if (!is.lcc(object))
    stop("The provided object is not of class 'lcc'.", call. = FALSE)
  
  if (!missing(...)) {
    models <- list(object$model, ...)
    lls <- lapply(models, logLik)
    vals <- sapply(lls, function(el) {
      nObs <- attr(el, "nobs")
      c(AIC = -2 * as.numeric(el) + k * attr(el, "df"), 
        nObs = if (is.null(nObs)) NA_integer_ else nObs)
    })
    if (length(unique(na.omit(vals["nObs", ]))) > 1)
      warning("models are not all fitted to the same number of observations")
    
    val <- data.frame(AIC = vals["AIC", ], row.names = names(models))
    val
  } else {
    AIC(object$model, k = k)
  }
}

##' @rdname AIC.lcc
##' @title Bayesian Information Criterion for an \code{lcc} Object
##' @method BIC lcc
##' @aliases BIC.lcc
##'
##' @description 
##' Calculates the Bayesian Information Criterion (BIC) for a fitted
##' longitudinal concordance correlation model represented by an \code{lcc} object.
##' BIC is used for model selection, with lower values indicating a better model.
##'
##' @param object An object of class \code{lcc}, representing a fitted 
##'   longitudinal concordance correlation function.
##' @param ... Optional arguments passed to the underlying \code{BIC}
##'   function from the \code{stats} package.
##'
##' @details 
##' The function computes BIC as a measure of the trade-off between model fit and
##' complexity. It is particularly useful for comparing models with different numbers
##' of parameters. For more information, refer to the documentation for \code{\link{BIC}}.
##'
##' @seealso \code{\link[lcc]{lcc}}, \code{\link{summary.lcc}},
##'   \code{\link{coef.lcc}}, \code{\link{vcov.lcc}}, \code{\link{AIC.lcc}}
##'
##' @importFrom stats BIC logLik nobs na.omit
##'
##' @examples
##' \dontrun{
##' attach(simulated_hue)
##' fm6 <- lcc(data = simulated_hue, subject = "Fruit",
##'            resp = "Hue", method = "Method", time = "Time",
##'            qf = 2, qr = 1, components = TRUE,
##'            time_lcc = list(n=50, from=min(Time), to=max(Time)))
##' AIC(fm6)
##' BIC(fm6)
##' }
##'
##' @export

BIC.lcc <- function(object, ...) {
  if (!is.lcc(object))
    stop("The provided object is not of class 'lcc'.", call. = FALSE)
  
  if (!missing(...)) {
    models <- list(object$model, ...)
    lls <- lapply(models, logLik)
    vals <- sapply(lls, function(el) {
      nObs <- attr(el, "nobs")
      c(logLik = as.numeric(el), df = attr(el, "df"), 
        nObs = if (is.null(nObs)) NA_integer_ else nObs)
    })
    if (length(unique(na.omit(vals["nObs", ]))) > 1)
      warning("models are not all fitted to the same number of observations")
    
    val <- data.frame(df = vals["df", ], BIC = -2 * vals["logLik", ] + log(vals["nObs", ]) * vals["df", ],
                      row.names = names(models))
    val
  } else {
    BIC(object$model)
  }
}

#=======================================================================
# ranef
# =======================================================================
##' @rdname ranef.lcc
##' @title Extract Random Effects from an \code{lcc} Model
##' @usage \method{ranef}{lcc}(object, ...)
##' @method ranef lcc
##' @aliases ranef.lcc
##'
##' @description 
##' Extracts the estimated random effects from a fitted longitudinal concordance
##' correlation model represented by an \code{lcc} object. The function returns
##' a data frame with rows corresponding to different groups at a specified level
##' and columns representing the random effects.
##'
##' @param object An object inheriting from class \code{lcc},
##'   representing a fitted longitudinal concordance correlation function.
##' @param ... Optional arguments passed to the \code{ranef.lme}
##'   function from the \code{nlme} package.
##'
##' @details 
##' This function is useful for examining the random effects associated with
##' groups or subjects in the model. For a detailed explanation of these effects,
##' see the documentation for \code{\link{nlme}} objects.
##'
##' @seealso \code{\link[lcc]{lcc}}, \code{\link{coef.lcc}},
##'
##' @importFrom nlme ranef
##'
##' @examples
##' \dontrun{
##' fm1 <- lcc(data = hue, subject = "Fruit", resp = "H_mean",
##'            method = "Method", time = "Time", qf = 2, qr = 2)
##' ranef(fm1)
##' }
##' @export

ranef.lcc <- function(object, ...) {
  if (!is.lcc(object))
    stop("The provided object is not of class 'lcc'.", call. = FALSE)
  
  randomEffects <- ranef(object$model, ...)
  
  # Clean up the column names for better readability
  replacements <- list(
    "fmla." = "",
    "\\(time, degree = qr, raw = TRUE\\)" = "",
    "rand" = "Random",
    "poly" = "Time"
  )
  
  for (pattern in names(replacements)) {
    replacement <- replacements[[pattern]]
    colnames(randomEffects) <- gsub(pattern, replacement, colnames(randomEffects))
  }
  
  class(randomEffects) <- c("ranef.lcc", "data.frame")
  randomEffects
}

#=======================================================================
# logLik
#=======================================================================
##' @rdname logLik.lcc
##' @title Extract Log-Likelihood of an \code{lcc} Object
##' @usage \method{logLik}{lcc}(object, ..., REML)
##' @method logLik lcc
##' @aliases logLik.lcc
##'
##' @description If \code{REML=TRUE}, the default, returns the
##'   restricted log-likelihood value of the linear mixed-effects model;
##'   else the log-likelihood value
##'
##' @param object an object inheriting from class \code{lcc},
##'   representing a fitted longitudinal concordance correlation
##'   function.
##'
##' @param REML an optional logical value.  If \code{TRUE} the
##'   restricted log-likelihood is returned, else, if \code{FALSE}, the
##'   log-likelihood is returned.
##' @param ... further arguments passed to \code{\link{logLik}}.
##'
##' @details See methods for \code{\link{nlme}} objects to get more
##'   details.
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @importFrom stats logLik
##'
##' @seealso \code{\link[lcc]{lcc}}, \code{\link{summary.lcc}}
##'
##' @examples
##'
##' \dontrun{
##' fm1<-lcc(data = hue, subject = "Fruit", resp = "H_mean",
##'          method = "Method", time = "Time", qf = 2, qr = 2)
##' logLik(fm1)
##' }
##'
##' @export

logLik.lcc <- function(object, ..., REML) {
   if (!is.lcc(object))
     stop("use only with \"lcc\" objects" , call. = FALSE)
   logLik(object$model,  REML = REML, ...)
}

#=======================================================================
# ANOVA
# =======================================================================
##' @rdname anova.lcc
##' @title Compare Likelihoods of Fitted Models from an \code{lcc} Object
##' @usage \method{anova}{lcc}(object, ..., test = TRUE, type = c("sequential", "marginal"), 
##'   adjustSigma = TRUE, verbose = FALSE)
##' @method anova lcc
##' @aliases anova.lcc
##'
##' @description 
##' Compares the fit of different longitudinal concordance correlation models 
##' (lcc objects). When comparing multiple models, the function returns a 
##' data frame with degrees of freedom, log-likelihood, AIC, and BIC for each model. 
##' For a single model, it returns F-values and P-values for fixed terms in the model.
##'
##' @param object An object inheriting from class \code{lcc} or \code{lme}.
##' @param ... Other optional fitted model objects inheriting from classes "lcc" or "lme".
##' @param test Logical; if \code{TRUE}, performs likelihood ratio tests to compare models. 
##'   Defaults to \code{TRUE}.
##' @param type Character string specifying the type of sum of squares for F-tests. 
##'   Options are "sequential" or "marginal". Defaults to "sequential".
##' @param adjustSigma Logical; if \code{TRUE}, adjusts the residual standard error 
##'   for maximum likelihood estimation. Defaults to \code{TRUE}.
##' @param verbose Logical; if \code{TRUE}, prints additional model details. 
##'   Defaults to \code{FALSE}.
##'
##' @details 
##' This function is an adaptation from \code{\link{anova.lme}}. 
##' It assesses whether the addition of terms significantly improves model fit.
##'
##' @seealso \code{\link[lcc]{lcc}}, \code{\link{summary.lcc}}
##'
##' @importFrom stats formula pchisq pf terms
##' @importFrom utils head tail
##'
##' @examples
##' \dontrun{
##' fm1.aov <- lcc(data = hue, subject = "Fruit", resp = "H_mean", method = "Method", 
##'                time = "Time", qf = 2, qr = 1)
##' fm2.aov <- update(fm1.aov, qr = 2)
##' anova(fm1.aov, fm2.aov)
##' }
##'
##' @examples
##' \dontrun{
##' fm3.aov <- update(fm2.aov, REML = FALSE)
##' fm4.aov <- update(fm2.aov, REML = FALSE, qf = 3)
##' anova(fm3.aov, fm4.aov)
##' }
##'
##' @examples
##' \dontrun{
##' fm5.aov <- update(fm2.aov, var.class = varExp, weights.form = "time")
##' anova(fm1.aov, fm2.aov, fm5.aov)
##' }
##'
##' @export
anova.lcc <- function(object, ..., test = TRUE, type = c("sequential", "marginal"),
                       adjustSigma = TRUE, verbose = FALSE) {
  fixSig <- attr(object$model$modelStruct, "fixedSigma")
  fixSig <- !is.null(fixSig) && fixSig
  dots <- list(...)
  if ((rt <- length(dots) + 1L) == 1L) {
    if (!inherits(object$model, "lme")) {
      stop("object must inherit from class \"lme\" ")
    }
    vFix <- attr(object$model$fixDF, "varFixFact")
    if (adjustSigma && object$model$method == "ML")
      vFix <- sqrt(object$model$dims$N/(object$model$dims$N - ncol(vFix))) *
        vFix
    c0 <- solve(t(vFix), fixef(object$model))
    assign <- attr(object$model$fixDF, "assign")
    nTerms <- length(assign)
    type <- match.arg(type)
    Fval <- Pval <- double(nTerms)
    nDF <- integer(nTerms)
    dDF <- object$model$fixDF$terms
    for (i in 1:nTerms) {
      nDF[i] <- length(assign[[i]])
      if (type == "sequential") {
        c0i <- c0[assign[[i]]]
      }
      else {
        c0i <- c(qr.qty(qr(vFix[, assign[[i]], drop = FALSE]),
                        c0))[1:nDF[i]]
      }
      Fval[i] <- sum(c0i^2)/nDF[i]
      Pval[i] <- 1 - pf(Fval[i], nDF[i], dDF[i])
    }
    aod <- data.frame(numDF = nDF, denDF = dDF, `F-value` = Fval,
                      `p-value` = Pval, check.names = FALSE)
    rownames(aod) <- names(assign)
    attr(aod, "rt") <- rt
  }
  else {
    ancall <- sys.call()
    ancall$verbose <- ancall$test <- ancall$type <- NULL
    object <- list(object$model, ...)
    termsClass <- vapply(object, data.class, "")
    valid.cl <- c("lme", "lcc")
    if (!all(match(termsClass, valid.cl, 0))) {
      valid.cl <- paste0("\"", valid.cl, "\"")
      stop(gettextf("objects must inherit from classes %s, or %s",
        paste(head(valid.cl, -1), collapse = ", "),
        tail(valid.cl, 1)), domain = NA)
    }
    for (i in 1:length(object)) {
      if (is.lcc(object[[i]])) {
        object[[i]] <- object[[i]]$model
      }else {
        object[[i]] <- object[[i]]
      }
    }
    getResponseFormula <-
      function(object)
      {
        ## Return the response formula as a one sided formula
        form <- formula(object)
        if (!(inherits(form, "formula") && (length(form) == 3))) {
          stop("'form' must be a two-sided formula")
        }
        eval(parse(text = paste("~", deparse(form[[2]]))))
      }
    resp <- vapply(object, function(el) deparse(getResponseFormula(el)[[2L]]),
                   "")
    subs <- as.logical(match(resp, resp[1L], FALSE))
    if (!all(subs))
      warning("some fitted objects deleted because response differs from the first model")
    if (sum(subs) == 1)
      stop("first model has a different response from the rest")
    object <- object[subs]
    rt <- length(object)
    termsModel <- lapply(object, function(el) formula(el)[-2])
    estMeth <- vapply(object, function(el) if (is.null(val <- el[["method"]]))
      NA_character_
    else val, "")
    if (length(uEst <- unique(estMeth[!is.na(estMeth)])) >
      1) {
      stop("all fitted objects must have the same estimation method")
    }
    estMeth[is.na(estMeth)] <- uEst
    REML <- uEst == "REML"
    if (REML) {
      aux <- vapply(termsModel, function(el) {
        tt <- terms(el)
        val <- paste(sort(attr(tt, "term.labels")),
          collapse = "&")
        if (attr(tt, "intercept") == 1)
          paste(val, "(Intercept)", sep = "&")
        else val
      }, ".")
      if (length(unique(aux)) > 1) {
        warning("fitted objects with different fixed effects. REML comparisons are not meaningful.")
      }
    }
    termsCall <- lapply(object, function(el) {
      if (is.null(val <- el$call) && is.null(val <- attr(el,
        "call")))
        stop("objects must have a \"call\" component or attribute")
      val
    })
    termsCall <- vapply(termsCall, function(el) paste(deparse(el),
      collapse = ""), "")
    aux <- lapply(object, logLik, REML)
    if (length(unique(vapply(aux, attr, 1, "nall"))) > 1) {
      stop("all fitted objects must use the same number of observations")
    }
    dfModel <- vapply(aux, attr, 1, "df")
    logLik <- vapply(aux, c, 1.1)
    aod <- data.frame(call = termsCall, Model = 1:rt, df = dfModel,
      AIC = vapply(aux, AIC, 1), BIC = vapply(aux, BIC,
        1), logLik = logLik, check.names = FALSE)
    if (test) {
      ddf <- diff(dfModel)
      if (sum(abs(ddf)) > 0) {
        effects <- rep("", rt)
        for (i in 2:rt) {
          if (ddf[i - 1] != 0) {
            effects[i] <- paste(i - 1, i, sep = " vs ")
          }
        }
        pval <- rep(NA, rt - 1)
        ldf <- as.logical(ddf)
        lratio <- 2 * abs(diff(logLik))
        lratio[!ldf] <- NA
        pval[ldf] <- pchisq(lratio[ldf], abs(ddf[ldf]),
          lower.tail = FALSE)
        aod <- data.frame(aod, Test = effects, L.Ratio = c(NA,
          lratio), `p-value` = c(NA, pval), check.names = FALSE,
          stringsAsFactors = TRUE)
      }
    }
    ## local function for complete deparsing
    c_deparse <- function(...) paste(deparse(..., width.cutoff=500),
                                     collapse="")
    row.names(aod) <- vapply(as.list(ancall[-1L]), c_deparse,
      "")
    attr(aod, "rt") <- rt
    attr(aod, "verbose") <- verbose
  }
  class(aod) <- c("anova.lcc", "data.frame")
  aod
}


##' @rdname print.anova.lcc
##' @title  Print the Anova of an \code{lcc} Object
##' @usage
##' \method{print}{anova.lcc}(x, verbose, ...)
##' @method print anova.lcc
##' @aliases print.anova.lcc
##' @description Method print for the \code{anova.lcc}.
##'
##' @param x an object inheriting from class
##'   \code{\link[lcc]{anova.lcc}}, representing a fitted longitudinal
##'   concordance correlation function.
##'
##' @param verbose an optional logical value used to control the amount
##'   of printed output. If \code{TRUE}, the calling sequences for each fitted
##'   model object are printed with the rest of the output, being omitted
##'   if \code{verbose = FALSE}. Defaults to \code{FALSE}.
##'
##' @param ... further arguments passed to \code{\link{print}}.
##'
##' @details Modified from \code{\link{anova.lme}}. For more details see
##' methods for \code{\link{nlme}}.
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @seealso \code{\link{summary.lcc}}, \code{\link{lccPlot}},
##'   \code{\link[lcc]{lcc}}
##'
##' @examples
##'
##' \dontrun{
##' ## Second degree polynomial model with random intercept, slope and
##' ## quadratic term
##' fm1<-lcc(data = hue, subject = "Fruit", resp = "H_mean",
##'          method = "Method", time = "Time", qf = 2, qr = 2)
##' print(anova(fm1))
##' }
##'
##' @export
print.anova.lcc <- function(x, verbose = attr(x, "verbose"), ...) {
  rt <- attr(x, "rt")
  
  if (rt == 1) {
    # Adjust label if exists
    if (!is.null(lab <- attr(x, "label"))) {
      cat(lab)
    }
    
    # Format p-values and F-values
    x[, "p-value"] <- sapply(x[, "p-value"], function(p) ifelse(p < 0.0001, "<.0001", round(p, 4)))
    x[, "F-value"] <- format(zapsmall(x[, "F-value"]))
    
    # Print the data frame
    print(as.data.frame(x), ...)
  } else {
    if (verbose) {
      cat("Call:\n")
      objNames <- row.names(x)
      for (i in 1:rt) {
        cat(" ", objNames[i], ":\n", sep = "")
        cat("  ", as.character(x[i, "call"]), "\n")
      }
      cat("\n")
    }
    
    # Format other columns
    for (colName in names(x)[-1]) {
      x[, colName] <- formatColumn(x[, colName], colName)
    }
    
    # Print the data frame
    print(as.data.frame(x[,-1]), ...)
  }
  
  invisible(x)
}

#' Format Columns for Print
#'
#' This internal helper function is used to format the columns of a data frame
#' for printing, specifically for use within the `print.anova.lcc` function. It
#' applies special formatting rules based on the column name, such as rounding
#' and special handling of small p-values.
#'
#' @param column A vector representing a column from a data frame.
#' @param colName A string indicating the name of the column, which determines
#'   the formatting rules to be applied.
#'
#' @return A vector with the same length as `column`, where each element has been
#'   formatted according to the column-specific rules.
#'
#' @details The function specifically handles the following columns:
#'   - "p-value": Rounds the values to four decimal places, and represents
#'     values less than 0.0001 as "<.0001".
#'   - "AIC", "BIC", "logLik", "L.Ratio": Applies `zapsmall` for formatting.
#'   Other columns are returned without changes.
#'
#' @examples
#' data <- data.frame(
#'   pvalue = c(0.00005, 0.0234, 0.5),
#'   AIC = c(123.4567, 234.5678, 345.6789)
#' )
#' data$pvalue <- formatColumn(data$pvalue, "p-value")
#' data$AIC <- formatColumn(data$AIC, "AIC")
#'
#' @export
formatColumn <- function(column, colName) {
  formattedCol <- column
  if (colName == "p-value") {
    formattedCol <- sapply(column, function(p) ifelse(p < 0.0001, "<.0001", format(round(p, 4))))
  } else if (colName %in% c("AIC", "BIC", "logLik", "L.Ratio")) {
    formattedCol <- format(zapsmall(column))
  }
  formattedCol
}

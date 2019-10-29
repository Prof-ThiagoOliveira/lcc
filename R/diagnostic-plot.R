#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: diagnostic-plot.R                                             #
# Contains: plot.lcc function                                         #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 28/10/2019                                           #
# Last update: 28/10/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @rdname plot.lcc
##' @method plot lcc
##' @title Plot an \code{lcc} object.
##'
##' @usage 
##' \method{plot}{lcc}(x, which = c(1L:6L),
##'      caption = list("Relsiduals vs Fitted",
##'                     "Residuals vs Time",
##'                     "Residuals by Subject",
##'                     "Observed values vs Fitted values",
##'                     "Normal Q-Q Plot (Conditional residuals)",
##'                     "Normal Q-Q Plot (Random effects)"),
##'      sub.caption =  NULL,  main = NULL,
##'      panel = if(add.smooth) panel.smooth else points,
##'      add.smooth = getOption("add.smooth"),
##'      ask = prod(par("mfcol")) < length(which) && dev.interactive(),
##'      id.n = 3, labels.id = names(residuals(x)),
##'      label.pos = c(4, 2), cex.id = 0.75, cex.caption = 1,
##'      cex.oma.man = 1.25, ...)
##'
##' @method plot lcc
##'
##' @aliases plot.lcc
##'
##' @description Diagnostic plots for conditional error and random
##' effects from the linear mixed-effects fit are obtained. Six plots
##' plots (selectable by 'which') are currently available: a plot of
##' residuals against fitted values, a plot of residuals against
##' time variable, a boxplot of residuals by subject,
##' a plot of observerd values against fitted values, a normal Q-Q plot
##' with simulation envelopes based on conditional error,  and a normal
##' Q-Q plot with simulation envelopes based on the random effects. By
##' default, all plots are provided.
##'
##' @param x an object inheriting from class \code{\link[lcc]{lcc}},
##'   representing a fitted longitudinal concordance correlation
##'   function.
##' @param which if a subset of the plots is required, specify a subset
##'   of the numbers from 1 to 6.
##' @param caption captions to appear above the plots. Vector or list of
##'   valid graphics annotations is required. All captions can be
##'   supressed using '""' or \code{NA}.
##' @param sub.caption common sub-title (at bottom). Default to
##'   \code{NULL}.
##' @param main The main title (on top) above the caption.
##' @param panel panel function. If \code{add.smooth = TRUE},
##'   \code{panel.smooth} is used rather than \code{points}.
##' @param add.smooth logical indicating if smoother should be added to
##'   most plots; see also \code{panel} above. Default to \code{TRUE}.
##' @param ask logical; if \code{TRUE}, the default, the user is _ask_ed
##'   before each plot, see \code{\link[graphics]{par}}.
##' @param id.n number of points to be labelled is the first three
##'   plots, starting with the most extreme.
##' @param labels.id vector of labels, from which the labels for extreme
##'   points will be chosen. Default to \code{NULL} (uses observation
##'   numbers).
##' @param label.pos positioning of labels, for the left half and right
##'   half of the graph respectively, for plots 1-3.
##' @param cex.id magnification of point label.
##' @param cex.caption controls the size of \code{caption}.
##' @param cex.oma.man controls the size of the \code{sub.caption} only
##'   if that is _above_ the figures when there is more than one.
##' @param ... further graphical parameters from 'par'.
##'
##' @details The Q-Q plot uses the normalized residuals. The
##'   standardized residuals is pre-multiplied by the inverse
##'   square-root factor of the estimated error correlation matrix while
##'   the random effects is pre-multiplied by the inverse square root of
##'   the estimated variances obtained from matrix G. The simulate
##'   envelopes are obtained from package hnp (Moral et al.,  2018).
##'   
##'   Code partially adapted from \code{\link[stats]{plot.lm}}.
##'
##' @importFrom hnp hnp
##'
##' @importFrom nlme getVarCov ranef
##' 
##' @importFrom grDevices as.graphicsAnnot dev.flush dev.hold dev.interactive devAskNewPage extendrange n2mfrow
##' 
##' @importFrom graphics abline boxplot mtext panel.smooth par plot points strheight text title
##' 
##' @importFrom stats fitted residuals
##' 
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@usp.br}
##'
##' @seealso \code{\link{lccPlot}}, \code{\link[lcc]{lcc}},
##'   \code{mtext}, \code{text}, \code{plotmath}
##'
##' @examples
##' data(hue)
##' ## Second degree polynomial model with random intercept, slope and
##' ## quadratic term
##' fm1 <- lcc(dataset = hue, subject = "Fruit", resp = "H_mean",
##'            method = "Method", time = "Time", qf = 2, qr = 2)
##' plot(fm1)
##' @export

plot.lcc <- function(x, which = c(1L:6L),
           caption = list("Relsiduals vs Fitted",
                          "Residuals vs Time",
                          "Residuals by Subject",
                          "Observed values vs Fitted values",
                          "Normal Q-Q Plot (Conditional residuals)",
                          "Normal Q-Q Plot (Random effects)"),
           sub.caption =  NULL,  main = NULL,
           panel = if(add.smooth) panel.smooth else points,
           add.smooth = getOption("add.smooth"),
           ask = prod(par("mfcol")) < length(which) && dev.interactive(),
           id.n = 3, labels.id = names(residuals(x)),
           label.pos = c(4, 2), cex.id = 0.75, cex.caption = 1,
           cex.oma.man = 1.25, ...)
  {
    if (!inherits(x, "lcc"))
      stop("use only with \"lcc\" objects")
    if(!is.numeric(which) || any(which < 1) || any(which > 6))
      stop("'which' must be in 1:6")
    #-------------------------------------------------------------------
    show <- rep(FALSE, 6)
    show[which] <- TRUE
    #-------------------------------------------------------------------
    # information from lme model
    #-------------------------------------------------------------------
    model <- x$model
    r <- residuals(model)
    r_norm <- residuals(model, type = "normalized")
    yh <- fitted(model)
    n <- length(r)
    time <- model$data$time
    #-------------------------------------------------------------------
     if(id.n > 0L) { ## label the largest residuals
       if(is.null(labels.id))
         labels.id <- paste(1L:n)
       iid <- 1L:id.n
       show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
       text.id <- function(x, y, ind, adj.x = TRUE) {
         if (is.factor(x)) {
           labpos <-
             if(adj.x) label.pos[1+as.numeric(y > mean(range(y)))] else 3
         }else {
           labpos <-
             if(adj.x) label.pos[1+as.numeric(x > mean(range(x)))] else 3
         }
         text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE,
              pos = labpos, offset = 0.25)
       }
    }
    #-------------------------------------------------------------------
    getCaption <- function(k) # allow caption = "" , plotmath etc
      if(length(caption) < k) NA_character_ else as.graphicsAnnot(caption[[k]])
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    #-------------------------------------------------------------------
    # Individual plots
    #-------------------------------------------------------------------
    if (show[1L]) {
      l.fit <- "Fitted values"
      ylim <- range(r, na.rm=TRUE)
      if(id.n > 0)
        ylim <- extendrange(r = ylim, f = 0.08)
      dev.hold()
      plot(yh, r, xlab = l.fit, ylab = "Residuals", main = main,
           ylim = ylim, type = "n", ...)
      panel(yh, r, ...)
      if (one.fig)
        title(sub = sub.caption, ...)
      mtext(getCaption(1), 3, 0.25, cex = cex.caption)
      if(id.n > 0) {
        y.id <- r[show.r]
        y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
        text.id(yh[show.r], y.id, show.r)
      }
      abline(h = 0, lty = 3, col = "gray")
      dev.flush()
    }
    #===================================================================
    if (show[2L]) {
      l.fit <- "Time"
      ylim <- range(r, na.rm=TRUE)
      if(id.n > 0)
        ylim <- extendrange(r = ylim, f = 0.08)
      dev.hold()
      plot(time, r, xlab = l.fit, ylab = "Residuals", main = main,
           ylim = ylim, type = "n", ...)
      panel(time, r, ...)
      if (one.fig)
      title(sub = sub.caption, ...)
      mtext(getCaption(2), 3, 0.25, cex = cex.caption)
      if(id.n > 0) {
        y.id <- r[show.r]
        y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
        text.id(time[show.r], y.id, show.r)
      }
      abline(h = 0, lty = 3, col = "gray")
      dev.flush()
    }
    #===================================================================
    if (show[3L]) {
      Subject <- model$data$subject
      boxplot(r ~ Subject,  ylab = "Residuals", main = main, ...)
      if (one.fig)
        title(sub = sub.caption, ...)
      mtext(getCaption(3), 3, 0.25, cex = cex.caption)
      if(id.n > 0) {
        y.id <- r[show.r]
        y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
        text.id(Subject[show.r], y.id, show.r)
      }
      abline(h = 0, lty = 3, col = "blue")
    }
    #===================================================================
    if (show[4L]) {
      Response <- model$data$resp
      plot(Response ~  yh, ylab = "Observed Values",
           xlab = "Fitted Values", main = main, ...)
      if (one.fig)
        title(sub = sub.caption, ...)
      mtext(getCaption(4), 3, 0.25, cex = cex.caption)
      abline(0, 1)
    }
    #===================================================================
    if (show[5L]) {
        hnp(r_norm, scale = TRUE, halfnormal = FALSE, print.on = TRUE,
                 main = main, ...)
        if (one.fig)
          title(sub = sub.caption, ...)
        mtext(getCaption(5), 3, 0.25, cex = cex.caption)
    }
    #===================================================================
    if (show[6L]) {
      vars <- sqrt(diag(getVarCov(model)))
      ranefs <- re <- as.matrix(ranef(model))
      re <- ranefs %*% diag(1 / vars)
      ncol.re <- ncol(re)
      if (ncol.re == 1) {
        if (is.null(main)) {
          main <- "Random effect (Intercept)"
        }
        hnp(re, scale = TRUE, halfnormal = FALSE, print.on = TRUE,
                   main = main, ...)
      }else {
        par(mfrow = rev(n2mfrow(ncol.re)))
        for (i in 1:ncol.re) {
          hnp(re[, i], scale = TRUE, halfnormal = FALSE, print.on = TRUE,
                   main = main, ...)
          if (one.fig)
            title(sub = sub.caption, ...)
          mtext(getCaption(6), 3, 0.25, cex = cex.caption)
          mtext(paste0("b", i - 1, "i"), 1, -1.5, cex = cex.caption)
        }
      }
    }
    par(mfrow = c(1,1))
    invisible()
  }

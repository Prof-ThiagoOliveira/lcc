#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: lccModel.R                                                    #
# Contains: lccModel function                                         #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal function to fit a linear mixed-effects model
##'
##' @description Internal helper to fit a linear mixed-effects model,
##'   following the formulation described in Laird and Ware (1982).
##'   See \code{\link[nlme]{lme}} for details.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @importFrom nlme lmeControl pdSymm
##' @importFrom stats model.matrix
##' @importFrom utils capture.output
##'
##' @references Laird, N. M.; Ware, J. H. (1982). Random-effects models
##'   for longitudinal data. \emph{Biometrics}, 38, 963–974.
##'
##' @references Pinheiro, J. C.; Bates, D. M. (1996). Unconstrained
##'   parametrizations for variance–covariance matrices.
##'   \emph{Statistics and Computing}, 6, 289–296.
##'
##' @references Pinheiro, J. C.; Bates, D. M. (2000). Mixed-effects
##'   models in S and S-PLUS. Springer.
##'
##' @keywords internal
lccModel <- function(dataset, resp, subject, method, time, qf, qr,
                     interaction, covar, gs = NULL, var.class = NULL,
                     weights.form, lme.control = NULL, method.init,
                     pdmat) {
  if (is.null(lme.control)) {
    lme.control <- lmeControl()
  }
  
  ## ------------------------------------------------------------------
  ## 1. Prepare data and polynomial basis
  ## ------------------------------------------------------------------
  Data <- dataBuilder(
    dataset = dataset,
    resp    = resp,
    subject = subject,
    method  = method,
    time    = time,
    gs      = gs
  )
  
  Poly <- with(Data, poly(time, degree = qf, raw = TRUE))
  
  ## Base fixed-effects design
  build_fixed_base <- function() {
    if (isTRUE(interaction)) {
      model.matrix(~ method * Poly, Data)
    } else {
      model.matrix(~ method + Poly, Data)
    }
  }
  
  fixed <- build_fixed_base()
  
  ## ------------------------------------------------------------------
  ## 2. Optional covariates in the fixed part
  ## ------------------------------------------------------------------
  if (!is.null(covar) && length(covar) > 0L) {
    pos <- pmatch(covar, names(Data))
    
    if (any(nap <- is.na(pos))) {
      stop(
        sprintf(
          ngettext(
            length(nap),
            "unrecognized 'covar' variable named %s ignored",
            "unrecognized 'covar' variable named %s ignored"
          ),
          paste(sQuote(covar[nap]), collapse = ", ")
        ),
        call. = FALSE
      )
      pos   <- pos[!nap]
      covar <- covar[!nap]
    }
    
    if (length(pos) > 0L) {
      COVAR <- vector("list", length(pos))
      for (i in seq_along(pos)) {
        mm <- model.matrix(~ Data[, pos[i]])[, -1, drop = FALSE]
        ## Preserve original naming convention (even though it is a bit odd)
        colnames(mm) <- paste(covar[[1]][i], levels(Data[, pos[i]])[-1])
        COVAR[[i]] <- mm
      }
      
      Data_covar <- do.call(cbind, COVAR)
      fixed      <- cbind(fixed, as.matrix(Data_covar))
    }
  }
  
  Data$fixed <- fixed
  
  ## ------------------------------------------------------------------
  ## 3. Random-effects design (qr = 0: random intercept; qr > 0: poly)
  ## ------------------------------------------------------------------
  if (qr > 0L) {
    Data$fmla.rand <- model.matrix(
      ~ poly(time, degree = qr, raw = TRUE),
      Data
    )
  }
  
  if (qr == 0L) {
    random_struct <- list(subject = pdSymm(form = ~ 1))
  } else {
    if (!is.function(pdmat)) {
      stop(
        "Available only for pdSymm, pdLogChol, pdDiag, pdIdent, pdCompSymm, and pdNatural.",
        call. = FALSE
      )
    }
    random_struct <- list(subject = pdmat(form = ~ fmla.rand - 1))
  }
  
  ## ------------------------------------------------------------------
  ## 4. Variance function (var.class) and weights.form
  ## ------------------------------------------------------------------
  weights_arg <- NULL
  if (!is.null(var.class)) {
    .form <- switch(
      weights.form,
      "time"       = ~ time,
      "method"     = ~ 1 | method,
      "time.ident" = ~ 1 | time,
      "both"       = ~ time | method
    )
    weights_arg <- var.class(form = .form)
  }
  
  ## ------------------------------------------------------------------
  ## 5. Call nlme::lme()
  ## ------------------------------------------------------------------
  lme_args <- list(
    fixed   = resp ~ fixed - 1,
    data    = Data,
    random  = random_struct,
    control = lme.control,
    method  = method.init
  )
  if (!is.null(weights_arg)) {
    lme_args$weights <- weights_arg
  }
  
  model.lme <- try(do.call(lme, lme_args), silent = TRUE)
  
  ## ------------------------------------------------------------------
  ## 6. Warnings and return object
  ## ------------------------------------------------------------------
  warning.count <- 0L
  if (inherits(model.lme, "try-error")) {
    warning.count <- 1L
    mes <- paste(capture.output(cat(model.lme[1])), collapse = " ")
  } else if (is.character(model.lme$apVar)) {
    warning.count <- 1L
    mes <- model.lme$apVar
  } else {
    mes <- NULL
  }
  
  lcc.fit <- list(
    model       = model.lme,
    summary     = summary(model.lme),
    q_f         = qf,
    data        = Data,
    wcount      = warning.count,
    lme.control = lme.control,
    message     = mes,
    method.init = method.init
  )
  class(lcc.fit) <- "lcc.fit"
  
  lcc.fit
}

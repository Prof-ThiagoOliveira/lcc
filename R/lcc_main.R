#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: lcc_main.R                                                    #
# Contains: lcc function                                              #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 22/11/2025                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Longitudinal concordance correlation (LCC) from polynomial mixed
##'   effects regression models using fixed effects and variance components
##'
##' @description The \code{lcc} function computes fitted values and
##'   non-parametric bootstrap confidence intervals for the longitudinal
##'   concordance correlation (LCC), longitudinal Pearson correlation (LPC),
##'   and longitudinal accuracy (LA).
##'
##'   These statistics are estimated from a polynomial mixed-effects model
##'   with flexible variance-covariance structures for random effects and
##'   variance functions that can model heteroscedastic within-subject
##'   errors, with or without time as a covariate.
##' @param data an object of class \code{data.frame}.
##'
##' @param resp character string. Name of the response variable in the
##'   data set.
##'
##' @param subject character string. Name of the subject variable in the
##'   data set.
##'
##' @param method character string. Name of the method variable in the
##'   data set. The first level of \code{method} is used as the
##'   gold-standard method.
##'
##' @param time character string. Name of the time variable in the data
##'   set.
##'
##' @param interaction logical. Indicates whether to estimate the
##'   interaction between \code{method} and \code{time}. If \code{TRUE}
##'   (the default), both main effects and their interaction are
##'   estimated. If \code{FALSE}, only the main effects of time and
##'   method are estimated.
##'
##' @param qf integer. Degree of the polynomial time trend, usually
##'   1, 2, or 3 (degree 0 is not allowed). Default is \code{qf = 1}.
##'
##' @param qr integer. Degree of the random-effects polynomial in time
##'   used to model subject-to-subject variation. Note that
##'   \code{qr = 0} specifies a random intercept (form
##'   \code{~ 1 | subject}); \code{qr = 1} specifies random intercept
##'   and slope (form \code{~ time | subject}). If \code{qr = qf = q},
##'   with \eqn{q \ge 1}, random effects at the subject level are added
##'   to all terms of the time polynomial regression (form
##'   \code{~ poly(time, q, raw = TRUE) | subject}). Default is
##'   \code{qr = 0}.
##'
##' @param covar character vector. Names of covariates to be included
##'   in the model as fixed effects. Defaults to \code{NULL}, meaning
##'   that no additional covariates are included.
##'
##' @param gs character string. Name of the level of \code{method} that
##'   represents the gold-standard. Defaults to the first level of
##'   \code{method}.
##'
##' @param pdmat standard classes of positive-definite matrix
##'   structures defined in \code{\link[nlme]{pdClasses}}. The
##'   available positive-definite matrix structures in \code{lcc} are
##'   \code{pdSymm} (the default), \code{pdLogChol}, \code{pdDiag},
##'   \code{pdIdent}, \code{pdCompSymm}, and \code{pdNatural}.
##'
##' @param var.class standard classes of variance functions used to
##'   model the variance structure of within-subject errors using
##'   covariates; see \code{\link[nlme]{varClasses}}. Defaults to
##'   \code{NULL}, which corresponds to homoscedastic within-subject
##'   errors. Available standard classes include:
##'   \describe{
##'     \item{\code{varIdent}:}{allows different variances according to
##'       the levels of a stratification variable.}
##'     \item{\code{varExp}:}{exponential function of the variance
##'       covariate; see \code{\link[nlme]{varExp}}.}
##'   }
##'
##' @param weights.form character string. A one-sided formula
##'   specifying a variance covariate and, optionally, a grouping
##'   factor for the variance parameters in \code{var.class}. If
##'   \code{var.class = varIdent}, the options \code{"method"} (form
##'   \code{~ 1 | method}) or \code{"time.ident"} (form
##'   \code{~ 1 | time}) must be used in \code{weights.form}. If
##'   \code{var.class = varExp}, the options \code{"time"} (form
##'   \code{~ time}) or \code{"both"} (form \code{~ time | method})
##'   must be used in \code{weights.form}.
##'
##' @param time_lcc list or \code{NULL}. Regular sequence for the time
##'   variable, merged with specific or experimental time values used
##'   for LCC, LPC, and LA predictions. Defaults to \code{NULL}. The
##'   list may contain the following components:
##'   \describe{
##'     \item{\code{time}:}{a vector of specific or experimental time
##'       values. The experimental time values are used by default.}
##'     \item{\code{from}:}{the starting (minimum) value of the time
##'       variable.}
##'     \item{\code{to}:}{the end (maximum) value of the time
##'       variable.}
##'     \item{\code{n}:}{integer specifying the desired length of the
##'       sequence. Values of \code{n} between 30 and 50 are usually
##'       adequate.}
##'   }
##'
##' @param ci logical. If \code{TRUE}, non-parametric bootstrap
##'   confidence intervals are calculated for the LCC, LPC, and LA
##'   statistics and printed in the output. Default is \code{FALSE}.
##'
##' @param ci.method character. Confidence-interval construction method:
##'   \code{"normal"} (default), \code{"percentile"}, or \code{"bca"}.
##'
##' @param boot.scheme character. Bootstrap resampling scheme. Defaults
##'   to \code{"np_case"} (subject-level case bootstrap). Other options:
##'   \describe{
##'     \item{\code{"np_case"}:}{resample subjects with replacement; keep
##'       original responses.}
##'     \item{\code{"np_case_resid_gr"}:}{case bootstrap; replace the
##'       response with fitted values plus residuals resampled from the
##'       global residual pool.}
##'     \item{\code{"np_case_resid_ir"}:}{case bootstrap; replace the
##'       response with fitted values plus residuals resampled within
##'       each subject.}
##'     \item{\code{"np_re_resid_gr"}:}{resample subject-specific fitted
##'       trajectories (includes random effects), then add residuals
##'       sampled from the pooled residuals.}
##'     \item{\code{"np_re_resid_ir"}:}{resample subject-specific fitted
##'       trajectories, then add residuals resampled within each
##'       subject.}
##'     \item{\code{"sp_case_pr"}:}{semiparametric case bootstrap;
##'       resample subjects, use their fitted trajectories, then add
##'       Gaussian noise with variance equal to the estimated residual
##'       variance.}
##'     \item{\code{"p_re_pr"}:}{fully parametric; simulate random
##'       effects from the estimated covariance and residuals from a
##'       Gaussian with the estimated residual variance, then generate
##'       responses via \eqn{X_i \hat{\beta} + Z_i u_i^* + \varepsilon_i^*}.}
##'   }
##'
##' @param alpha significance level. Default is \code{0.05}.
##'
##' @param nboot integer. Number of bootstrap samples. Default is
##'   \code{5000}.
##'
##' @param show.warnings logical. Controls the display of convergence
##'   warnings in the bootstrap samples. If \code{TRUE}, the indices of
##'   bootstrap samples with convergence errors are shown. If
##'   \code{FALSE} (the default), only the total number of convergence
##'   errors is reported.
##'
##' @param components logical. If \code{TRUE}, estimates and confidence
##'   intervals for LPC and LA are printed in the output. If
##'   \code{FALSE} (the default), only estimates and confidence
##'   intervals for the LCC statistic are provided.
##'
##' @param REML logical. If \code{TRUE} (the default), the model is fit
##'   by maximising the restricted log-likelihood. If \code{FALSE}, the
##'   full log-likelihood is maximised.
##'
##' @param lme.control list. Control values for the estimation
##'   algorithm, replacing the defaults of
##'   \code{\link[nlme]{lmeControl}} in the \pkg{nlme} package. Defaults
##'   to \code{NULL}. The returned list is passed as the \code{control}
##'   argument to \code{\link[nlme]{lme}}.
##'
##' @param numCore integer or \code{NULL}. Number of cores used in
##'   parallel during bootstrap computation. If \code{NULL} (default),
##'   the function attempts to use one fewer than the available cores;
##'   otherwise it uses the supplied value.
##'
##' @return An object of class \code{lcc}. The output is a list with
##'   the following components:
##'   \item{model}{summary of the polynomial mixed-effects regression
##'     model.}
##'   \item{Summary.lcc}{fitted values for LCC, or for LCC, LPC, and LA
##'     if \code{components = TRUE}; the concordance correlation
##'     coefficient (CCC) between methods at each sampled value of
##'     \code{time}, and the CCC between mixed-effects model
##'     predictions and observed data as a goodness-of-fit measure
##'     (gof).}
##'   \item{data}{the input data set.}
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@alumni.usp.br},
##'   Rafael de Andrade Moral,
##'   John Hinde
##'
##' @seealso \code{\link{summary.lcc}}, \code{\link{fitted.lcc}},
##'   \code{\link{print.lcc}}, \code{\link{lccPlot}},
##'   \code{\link{plot.lcc}}, \code{\link{coef.lcc}},
##'   \code{\link{ranef.lcc}}, \code{\link{vcov.lcc}},
##'   \code{\link{getVarCov.lcc}}, \code{\link{residuals.lcc}},
##'   \code{\link{AIC.lcc}}
##'
##' @references Lin, L. A concordance correlation coefficient to
##'   evaluate reproducibility. \emph{Biometrics}, 45(1), 255–268,
##'   1989.
##'
##' @references Oliveira, T. P.; Hinde, J.; Zocchi, S. S.
##'   Longitudinal concordance correlation function based on variance
##'   components: an application in fruit colour analysis.
##'   \emph{Journal of Agricultural, Biological, and Environmental
##'   Statistics}, 23(2), 233–254, 2018.
##'
##' @references Oliveira, T. P.; Moral, R. A.; Zocchi, S. S.;
##'   Demetrio, C. G. B.; Hinde, J. lcc: an R package to estimate the
##'   concordance correlation, Pearson correlation, and accuracy over
##'   time. \emph{PeerJ}, 8:e9850, 2020. DOI:10.7717/peerj.9850
##'
##' @keywords nlme ggplot2
##'
##' @examples
##' data(hue)
##' ## Second degree polynomial model with random intercept, slope and
##' ## quadratic term
##' fm1 <- lcc(data = hue, subject = "Fruit", resp = "H_mean",
##'            method = "Method", time = "Time", qf = 2, qr = 2)
##' print(fm1)
##' summary(fm1)
##' summary(fm1, type = "model")
##' lccPlot(fm1) +
##'   ylim(0, 1) +
##'   geom_hline(yintercept = 1, linetype = "dashed") +
##'   scale_x_continuous(breaks = seq(1, max(hue$Time), 2))
##'
##' ## Estimating longitudinal Pearson correlation and longitudinal
##' ## accuracy
##' fm2 <- update(fm1, components = TRUE)
##' summary(fm2)
##' lccPlot(fm2) +
##'   ylim(0, 1) +
##'   geom_hline(yintercept = 1, linetype = "dashed") +
##'   scale_x_continuous(breaks = seq(1, max(hue$Time), 2)) +
##'   theme_bw()
##'
##' \dontrun{
##' ## A grid of points as the Time variable for prediction
##' fm3 <- update(
##'   fm2,
##'   time_lcc = list(
##'     from = min(hue$Time),
##'     to   = max(hue$Time),
##'     n    = 40
##'   )
##' )
##' summary(fm3)
##' lccPlot(fm3) +
##'   ylim(0, 1) +
##'   geom_hline(yintercept = 1, linetype = "dashed") +
##'   scale_x_continuous(breaks = seq(1, max(hue$Time), 2)) +
##'   theme_bw()
##' }
##'
##' ## Including an exponential variance function using time as a
##' ## covariate
##' fm4 <- update(
##'   fm2,
##'   time_lcc    = list(from = min(hue$Time),
##'                      to   = max(hue$Time),
##'                      n    = 30),
##'   var.class   = varExp,
##'   weights.form = "time"
##' )
##' summary(fm4, type = "model")
##' fitted(fm4)
##' fitted(fm4, type = "lpc")
##' fitted(fm4, type = "la")
##' lccPlot(fm4) +
##'   geom_hline(yintercept = 1, linetype = "dashed")
##' lccPlot(fm4, type = "lpc") +
##'   geom_hline(yintercept = 1, linetype = "dashed")
##' lccPlot(fm4, type = "la") +
##'   geom_hline(yintercept = 1, linetype = "dashed")
##'
##' \dontrun{
##' ## Non-parametric confidence interval with 500 bootstrap samples
##' fm5 <- update(fm1, ci = TRUE, nboot = 500)
##' summary(fm5)
##' lccPlot(fm5) +
##'   geom_hline(yintercept = 1, linetype = "dashed")
##' }
##'
##' ## Comparing bootstrap schemes and CI methods (small nboot for example)
##' \dontrun{
##' set.seed(123)
##' fm_np_norm <- update(fm1, ci = TRUE, nboot = 100,
##'                      boot.scheme = "np_case", ci.method = "normal")
##' fm_np_pct  <- update(fm1, ci = TRUE, nboot = 100,
##'                      boot.scheme = "np_case", ci.method = "percentile")
##' fm_np_bca  <- update(fm1, ci = TRUE, nboot = 100,
##'                      boot.scheme = "np_case", ci.method = "bca")
##' fm_re_res  <- update(fm1, ci = TRUE, nboot = 100,
##'                      boot.scheme = "np_re_resid_gr", ci.method = "normal")
##'
##' lccPlot(fm_np_norm) + ggtitle("np_case / normal")
##' lccPlot(fm_np_pct)  + ggtitle("np_case / percentile")
##' lccPlot(fm_np_bca)  + ggtitle("np_case / bca")
##' lccPlot(fm_re_res)  + ggtitle("np_re_resid_gr / normal")
##' }
##'
##' ## Considering three methods of colour evaluation
##' \dontrun{
##' data(simulated_hue)
##' attach(simulated_hue)
##' fm6 <- lcc(
##'   data     = simulated_hue,
##'   subject  = "Fruit",
##'   resp     = "Hue",
##'   method   = "Method",
##'   time     = "Time",
##'   qf       = 2,
##'   qr       = 1,
##'   components = TRUE,
##'   time_lcc = list(
##'     n    = 50,
##'     from = min(Time),
##'     to   = max(Time)
##'   )
##' )
##' summary(fm6)
##' lccPlot(fm6, scales = "free")
##' lccPlot(fm6, type = "lpc", scales = "free")
##' lccPlot(fm6, type = "la", scales = "free")
##' detach(simulated_hue)
##' }
##'
##' ## Including an additional covariate in the linear predictor
##' ## (randomised block design)
##' \dontrun{
##' data(simulated_hue_block)
##' attach(simulated_hue_block)
##' fm7 <- lcc(
##'   data      = simulated_hue_block,
##'   subject   = "Fruit",
##'   resp      = "Hue",
##'   method    = "Method",
##'   time      = "Time",
##'   qf        = 2,
##'   qr        = 1,
##'   components = TRUE,
##'   covar     = c("Block"),
##'   time_lcc  = list(
##'     n    = 50,
##'     from = min(Time),
##'     to   = max(Time)
##'   )
##' )
##' summary(fm7)
##' lccPlot(fm7, scales = "free")
##' detach(simulated_hue_block)
##' }
##'
##' ## Testing the interaction effect between time and method
##' fm8 <- update(fm1, interaction = FALSE)
##' anova(fm1, fm8)
##'
##' \dontrun{
##' ## Using parallel computing with 3 cores, and set.seed(123) to
##' ## verify model reproducibility
##' set.seed(123)
##' fm9 <- lcc(
##'   data     = hue,
##'   subject  = "Fruit",
##'   resp     = "H_mean",
##'   method   = "Method",
##'   time     = "Time",
##'   qf       = 2,
##'   qr       = 2,
##'   ci       = TRUE,
##'   nboot    = 30,
##'   numCore  = 3
##' )
##'
##' ## Repeating the same model with the same seed
##' set.seed(123)
##' fm10 <- lcc(
##'   data     = hue,
##'   subject  = "Fruit",
##'   resp     = "H_mean",
##'   method   = "Method",
##'   time     = "Time",
##'   qf       = 2,
##'   qr       = 2,
##'   ci       = TRUE,
##'   nboot    = 30,
##'   numCore  = 3
##' )
##'
##' ## Verifying that fitted values and confidence intervals are
##' ## identical
##' identical(fm9$Summary.lcc$fitted, fm10$Summary.lcc$fitted)
##' }
##'
##' @export
lcc <- function(data, resp, subject, method, time,
                interaction   = TRUE,
                qf            = 1,
                qr            = 0,
                covar         = NULL,
                gs            = NULL,
                pdmat         = pdSymm,
                var.class     = NULL,
                weights.form  = NULL,
                time_lcc      = NULL,
                ci            = FALSE,
                boot.scheme   = "np_case",
                ci.method     = "normal",
                alpha         = 0.05,
                nboot         = 5000,
                show.warnings = FALSE,
                components    = FALSE,
                REML          = TRUE,
                lme.control   = NULL,
                numCore       = NULL) {
  
  # keep original call
  lcc_call <- match.call()
  
  if (is.null(numCore)) {
    cores <- parallel::detectCores(logical = FALSE)
    available <- if (is.finite(cores)) max(1L, cores - 1L) else 1L
    numCore <- min(available, 2L)
  } else if (numCore > 8L) {
    warn_general("Limiting 'numCore' to 8 to comply with parallel backend limits.")
    numCore <- 8L
  }
  
  # Resolve new arguments with backward compatibility
  ci.method <- check_choice(
    ci.method,
    choices = c("normal", "percentile", "bca"),
    arg = "ci.method"
  )
  boot.scheme <- check_choice(
    boot.scheme,
    choices = c(
      "np_case", "np_case_resid_gr", "np_case_resid_ir",
      "np_re_resid_gr", "np_re_resid_ir", "sp_case_pr", "p_re_pr"
    ),
    arg = "boot.scheme"
  )
  
  #-------------------------------------------------------------------
  # 1. Init: checks + resolve pdmat, var.class, REML
  #-------------------------------------------------------------------
  init_env <- init(
    var.class    = var.class,
    weights.form = weights.form,
    REML         = REML,
    qf           = qf,
    qr           = qr,
    pdmat        = pdmat,
    dataset      = data,
    resp         = resp,
    subject      = subject,
    method       = method,
    time         = time,
    gs           = gs,
    numCore      = numCore
  )
  
  pdmat      <- init_env$pdmat
  MethodREML <- init_env$MethodREML
  var.class  <- init_env$var.class
  
  #-------------------------------------------------------------------
  # 2. Fit model (lccModel)
  #-------------------------------------------------------------------
  model.info <- try(
    lccModel(
      dataset     = data,
      resp        = resp,
      subject     = subject,
      pdmat       = pdmat,
      method      = method,
      time        = time,
      qf          = qf,
      qr          = qr,
      interaction = interaction,
      covar       = covar,
      gs          = gs,
      var.class   = var.class,
      weights.form = weights.form,
      lme.control = lme.control,
      method.init = MethodREML
    ),
    silent = TRUE
  )
  
  # robust error handling for lccModel failure
  if (inherits(model.info, "try-error")) {
    cond <- attr(model.info, "condition")
    msg  <- if (!is.null(cond)) conditionMessage(cond) else as.character(model.info[1L])
    abort_lcc("Error in 'lccModel': {msg}")
  }
  
  #-------------------------------------------------------------------
  # 3. Check convergence via wcount
  #-------------------------------------------------------------------
  if (identical(model.info$wcount, "1")) {
    abort_lcc(model.info$message)
  }
  
  #-------------------------------------------------------------------
  # 4. Extract model + basic info
  #-------------------------------------------------------------------
  model        <- model.info$model
  q_f          <- qf
  q_r          <- qr
  lme.control  <- model.info$lme.control
  MethodREML   <- model.info$method.init
  tk           <- sort(unique(model.info$data$time))
  
  #-------------------------------------------------------------------
  # 5. Build fixed-effect contrasts (diffbeta) per method
  #-------------------------------------------------------------------
  lev.lab    <- levels(model.info$data$method)
  lev.method <- length(lev.lab)
  
  # create pattern "method<level>" repeated per polynomial degree q_f
  x <- y <- NULL  # for R CMD check (NSE in transform)
  lev.lab_df <- unique(merge(rep("method", q_f), lev.lab))
  lev.lab_df <- transform(lev.lab_df, newcol = paste(x, y, sep = ""))
  
  fx  <- fixef(model)
  pat <- lapply(seq_len(lev.method - 1L), function(i) {
    grep(lev.lab_df$newcol[i + 1L], names(fx))
  })
  
  # list of -beta_k for each non-reference method
  betas <- lapply(pat, function(idx) -fx[idx])
  
  #-------------------------------------------------------------------
  # 6. Internal calculations (LCC, LPC, LA, CI, etc.)
  #-------------------------------------------------------------------
  lcc_int <- lccInternal(
    model        = model,
    q_f          = q_f,
    q_r          = q_r,
    interaction  = interaction,
    tk           = tk,
    covar        = covar,
    pdmat        = pdmat,
    diffbeta     = betas,
    time_lcc     = time_lcc,
    ci           = ci,
    boot.scheme  = boot.scheme,
    ci.method    = ci.method,
    alpha        = alpha,
    nboot        = nboot,
    labels       = lev.lab_df,
    var.class    = var.class,
    weights.form = weights.form,
    show.warnings = show.warnings,
    components   = components,
    lme.control  = lme.control,
    method.init  = MethodREML,
    numCore      = numCore
  )
  
  #-------------------------------------------------------------------
  # 7. Build final lcc object
  #-------------------------------------------------------------------
  out <- list(
    "model"      = model,
    "Summary.lcc" = lcc_int[[1L]],
    "data"       = data,
    "plot_info"  = lcc_int[-1L],
    "call"       = lcc_call
  )
  class(out) <- "lcc"
  
  invisible(out)
}

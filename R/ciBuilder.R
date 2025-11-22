#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: ciBuilder.R                                                   #
# Contains: ciBuilder function                                        #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 22/11/2025                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal Function to Prepare the \code{ciCompute} Function.
##'
##' @description This is an internally called function used to prepare
##'   the \code{\link[lcc]{ciCompute}} function.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @keywords internal
ciBuilder <- function(model, nboot, q_f, q_r, interaction, covar, pdmat,
                      var.class, weights.form, show.warnings, tk,
                      diffbeta, ldb, tk.plot, tk.plot2, ci, percentileMet,
                      alpha, components, lme.control, method.init,
                      numCore) {
  #-------------------------------------------------------------------
  # 1. Bootstrap samples
  #-------------------------------------------------------------------
  Boot <- bootstrapSamples(
    nboot        = nboot,
    model        = model,
    q_f          = q_f,
    q_r          = q_r,
    interaction  = interaction,
    covar        = covar,
    pdmat        = pdmat,
    var.class    = var.class,
    weights.form = weights.form,
    show.warnings = show.warnings,
    lme.control  = lme.control,
    method.init  = method.init,
    numCore      = numCore,
    tk           = tk,
    ldb          = ldb,
    components   = components
  )
  
  LCC_Boot <- Boot$LCC_Boot
  LPC_Boot <- Boot$LPC_Boot
  Cb_Boot  <- Boot$Cb_Boot
  
  #-------------------------------------------------------------------
  # 3. Point estimates: LCC (rho)
  #-------------------------------------------------------------------
  if (ldb == 1L) {
    rho <- lccWrapper(
      model   = model,
      q_f     = q_f,
      n.delta = 1L,
      tk      = tk,
      diffbeta = as.numeric(diffbeta[[1L]])
    )
  } else {
    rho_list <- vector("list", ldb)
    for (i in seq_len(ldb)) {
      rho_list[[i]] <- lccWrapper(
        model   = model,
        q_f     = q_f,
        n.delta = 1L,
        tk      = tk,
        diffbeta = as.numeric(diffbeta[[i]])
      )
    }
    rho <- as.data.frame(do.call(cbind, rho_list))
  }
  
  #-------------------------------------------------------------------
  # 4. If we only want LCC, we're done after lcc_intervals()
  #-------------------------------------------------------------------
  if (!isTRUE(components)) {
    CI <- lcc_intervals(
      rho          = rho,
      tk.plot      = tk.plot,
      tk.plot2     = tk.plot2,
      ldb          = ldb,
      model        = model,
      ci           = ci,
      percentileMet = percentileMet,
      LCC_Boot     = LCC_Boot,
      alpha        = alpha
    )
    return(invisible(CI))
  }
  
  #-------------------------------------------------------------------
  # 5. Components = TRUE: LPC and LA
  #-------------------------------------------------------------------
  if (ldb == 1L) {
    rho.pearson <- lpcWrapper(
      model   = model,
      q_f     = q_f,
      tk      = tk,
      n.delta = 1L
    )
    
    Cb <- laWrapper(
      model   = model,
      q_f     = q_f,
      n.delta = 1L,
      tk      = tk,
      diffbeta = as.numeric(diffbeta[[1L]])
    )
  } else {
    rho_pearson_list <- vector("list", ldb)
    Cb_list          <- vector("list", ldb)
    
    for (i in seq_len(ldb)) {
      rho_pearson_list[[i]] <- lpcWrapper(
        model   = model,
        q_f     = q_f,
        tk      = tk,
        n.delta = 1L
      )
      
      Cb_list[[i]] <- laWrapper(
        model   = model,
        q_f     = q_f,
        n.delta = 1L,
        tk      = tk,
        diffbeta = as.numeric(diffbeta[[i]])
      )
    }
    
    rho.pearson <- as.data.frame(do.call(cbind, rho_pearson_list))
    Cb          <- as.data.frame(do.call(cbind, Cb_list))
  }
  
  #-------------------------------------------------------------------
  # 6. Delegate construction of CIs to ciCompute()
  #-------------------------------------------------------------------
  CI <- ciCompute(
    rho          = rho,
    rho.pearson  = rho.pearson,
    Cb           = Cb,
    tk.plot      = tk.plot,
    tk.plot2     = tk.plot2,
    ldb          = ldb,
    model        = model,
    ci           = ci,
    percentileMet = percentileMet,
    LCC_Boot     = LCC_Boot,
    LPC_Boot     = LPC_Boot,
    Cb_Boot      = Cb_Boot,
    alpha        = alpha
  )
  
  invisible(CI)
}


##' @title Internal Function to Compute the Non-Parametric Bootstrap
##'   Interval.
##'
##' @description This is an internally called function used to compute
##'   the non-parametric bootstrap interval.
##'
##' @usage NULL
##'
##' @details returns a matrix or list of matrix containing the
##'   non-parametric bootstrap interval.
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @importFrom stats quantile sd qnorm
##'
##' @keywords internal
ciCompute <- function(rho, rho.pearson, Cb, tk.plot, tk.plot2, ldb, model,
                      ci, percentileMet, LCC_Boot, LPC_Boot, Cb_Boot, alpha) {
  ## Fisher transform and inverse
  ZFisher     <- function(x) 0.5 * log((1 + x) / (1 - x))
  ZFisher_inv <- function(x) (exp(2 * x) - 1) / (exp(2 * x) + 1)
  
  ## Arcsin-sqrt transform for accuracy
  Arcsin     <- function(x) asin(sqrt((1 + x) / 2))
  Arcsin_inv <- function(x) 2 * (sin(x))^2 - 1
  
  percentile <- isTRUE(percentileMet) || identical(percentileMet, "TRUE")
  
  build_ci_metric <- function(boot_list, alpha, transform, inv_transform) {
    .build_ci_from_boot(
      boot_list     = boot_list,
      alpha         = alpha,
      transform     = if (!percentile) transform else NULL,
      inv_transform = if (!percentile) inv_transform else NULL,
      percentile    = percentile
    )
  }
  
  if (ldb == 1L) {
    ## Each *_Boot is list over bootstrap samples, numeric over time
    ENV.LCC <- build_ci_metric(LCC_Boot, alpha, ZFisher, ZFisher_inv)
    ENV.LPC <- build_ci_metric(LPC_Boot, alpha, ZFisher, ZFisher_inv)
    ENV.Cb  <- build_ci_metric(Cb_Boot,  alpha, Arcsin,  Arcsin_inv)
  } else {
    ## Each *_Boot[[b]] is a list over method combinations (ldb)
    ENV.LCC <- vector("list", ldb)
    ENV.LPC <- vector("list", ldb)
    ENV.Cb  <- vector("list", ldb)
    
    for (i in seq_len(ldb)) {
      LCC_i <- lapply(LCC_Boot, function(x) if (!is.null(x)) x[[i]] else NULL)
      LPC_i <- lapply(LPC_Boot, function(x) if (!is.null(x)) x[[i]] else NULL)
      Cb_i  <- lapply(Cb_Boot,  function(x) if (!is.null(x)) x[[i]] else NULL)
      
      ENV.LCC[[i]] <- build_ci_metric(LCC_i, alpha, ZFisher, ZFisher_inv)
      ENV.LPC[[i]] <- build_ci_metric(LPC_i, alpha, ZFisher, ZFisher_inv)
      ENV.Cb[[i]]  <- build_ci_metric(Cb_i,  alpha, Arcsin,  Arcsin_inv)
    }
  }
  
  CI <- list(
    rho      = rho,
    LPC      = rho.pearson,
    Cb       = Cb,
    ENV.LCC  = ENV.LCC,
    ENV.LPC  = ENV.LPC,
    ENV.Cb   = ENV.Cb,
    tk.plot  = tk.plot,
    tk.plot2 = tk.plot2,
    ldb      = ldb
  )
  
  CI
}

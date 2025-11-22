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
# Last update: 29/07/2019                                             #
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
  # 1. Bootstrap samples (models + fixed-effects differences)
  #-------------------------------------------------------------------
  Models <- bootstrapSamples(
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
    numCore      = numCore
  )
  
  boot_models  <- Models$Boot_Model
  boot_diffbet <- Models$Diffbetas
  
  #-------------------------------------------------------------------
  # 2. Bootstrap LCC (always needed)
  #-------------------------------------------------------------------
  LCC_Boot <- lccBootstrap(
    model_boot = boot_models,
    diff_boot  = boot_diffbet,
    ldb        = ldb,
    nboot      = nboot,
    tk         = tk,
    q_f        = q_f
  )
  
  #-------------------------------------------------------------------
  # 3. Point estimates: LCC (rho), always required
  #    Note: original code always uses n.delta = 1 here,
  #    so we keep that behaviour unchanged.
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
    # Same structure as original: data.frame with one column per method
    rho <- as.data.frame(do.call(cbind, rho_list))
  }
  
  #-------------------------------------------------------------------
  # 4. If we only want LCC, we’re done after lcc_intervals()
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
    # CI has components: list(rho = ..., ENV.LCC = ...)
    return(invisible(CI))
  }
  
  #-------------------------------------------------------------------
  # 5. Components = TRUE: need LPC and LA as well
  #-------------------------------------------------------------------
  LPC_Boot <- lpcBootstrap(
    model_boot = boot_models,
    ldb        = ldb,
    nboot      = nboot,
    tk         = tk,
    q_f        = q_f
  )
  Cb_Boot  <- laBootstrap(
    model_boot = boot_models,
    diff_boot  = boot_diffbet,
    ldb        = ldb,
    nboot      = nboot,
    tk         = tk,
    q_f        = q_f
  )
  
  #-------------------------------------------------------------------
  # 6. Point estimates: LPC (rho.pearson) and LA (Cb)
  #    Keep the original behaviour: n.delta = 1 in ciBuilder.
  #-------------------------------------------------------------------
  if (ldb == 1L) {
    rho_pearson <- lpcWrapper(
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
    
    rho_pearson <- as.data.frame(do.call(cbind, rho_pearson_list))
    Cb          <- as.data.frame(do.call(cbind, Cb_list))
  }
  
  #-------------------------------------------------------------------
  # 7. Compute CIs for LCC, LPC, and LA in one place
  #-------------------------------------------------------------------
  CI <- ciCompute(
    rho          = rho,
    rho.pearson  = rho_pearson,
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
                      ci, percentileMet, LCC_Boot, LPC_Boot, Cb_Boot,
                      alpha) {
  #-------------------------------------------------------------------
  # Transformations
  #-------------------------------------------------------------------
  ZFisher <- function(x) 0.5 * log((1 + x) / (1 - x))
  ZFisher_inv <- function(x) (exp(2 * x) - 1) / (exp(2 * x) + 1)
  
  Arcsin <- function(x) asin(sqrt(x))
  Arcsin_inv <- function(x) sign(x) * sin(x)^2
  
  # percentileMet may be logical or "TRUE"/"FALSE"
  percentile <- isTRUE(percentileMet) || identical(percentileMet, "TRUE")
  
  #-------------------------------------------------------------------
  # Helper to build CI for one metric and one method (or single-method)
  #-------------------------------------------------------------------
  build_ci_metric <- function(boot_list, alpha, transform, inv_transform) {
    .build_ci_from_boot(
      boot_list     = boot_list,
      alpha         = alpha,
      transform     = if (!percentile) transform else NULL,
      inv_transform = if (!percentile) inv_transform else NULL,
      percentile    = percentile
    )
  }
  
  #-------------------------------------------------------------------
  # ldb == 1: each *_Boot is a list of numeric vectors
  #-------------------------------------------------------------------
  if (ldb == 1L) {
    ENV.LCC <- build_ci_metric(
      boot_list     = LCC_Boot,
      alpha         = alpha,
      transform     = ZFisher,
      inv_transform = ZFisher_inv
    )
    
    ENV.LPC <- build_ci_metric(
      boot_list     = LPC_Boot,
      alpha         = alpha,
      transform     = ZFisher,
      inv_transform = ZFisher_inv
    )
    
    ENV.Cb <- build_ci_metric(
      boot_list     = Cb_Boot,
      alpha         = alpha,
      transform     = Arcsin,
      inv_transform = Arcsin_inv
    )
    
  } else {
    #---------------------------------------------------------------
    # ldb > 1: each *_Boot[[b]] is a list over methods (length ldb)
    #---------------------------------------------------------------
    ENV.LCC <- vector("list", ldb)
    ENV.LPC <- vector("list", ldb)
    ENV.Cb  <- vector("list", ldb)
    
    for (i in seq_len(ldb)) {
      # Extract i-th method across bootstrap replicates
      boot_LCC_i <- lapply(LCC_Boot, function(x) if (!is.null(x)) x[[i]] else NULL)
      boot_LPC_i <- lapply(LPC_Boot, function(x) if (!is.null(x)) x[[i]] else NULL)
      boot_Cb_i  <- lapply(Cb_Boot,  function(x) if (!is.null(x)) x[[i]] else NULL)
      
      ENV.LCC[[i]] <- build_ci_metric(
        boot_list     = boot_LCC_i,
        alpha         = alpha,
        transform     = ZFisher,
        inv_transform = ZFisher_inv
      )
      
      ENV.LPC[[i]] <- build_ci_metric(
        boot_list     = boot_LPC_i,
        alpha         = alpha,
        transform     = ZFisher,
        inv_transform = ZFisher_inv
      )
      
      ENV.Cb[[i]] <- build_ci_metric(
        boot_list     = boot_Cb_i,
        alpha         = alpha,
        transform     = Arcsin,
        inv_transform = Arcsin_inv
      )
    }
  }
  
  #-------------------------------------------------------------------
  # Return structure unchanged
  #-------------------------------------------------------------------
  CI.LCC <- list(
    "rho"     = rho,
    "ENV.LCC" = ENV.LCC,
    "LPC"     = rho.pearson,
    "ENV.LPC" = ENV.LPC,
    "Cb"      = Cb,
    "ENV.Cb"  = ENV.Cb
  )
  
  CI.LCC
}

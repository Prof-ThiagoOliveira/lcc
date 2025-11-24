
#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: lcc_workflow.R                                                #
# Contains: lccInternal function                                      #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 22/11/2025                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################
##' @title Internal Function to Prepare \code{lcc} Objects
##'
##' @description This is an internally called function used to prepare
##'   \code{lcc} objects for calculate the longitudinal concordance
##'   correlation, longitudinal Pearson correlation, longitudinal bias
##'   corrector, and plotting
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@alumni.usp.br} and Rafael de Andrade Moral,
##'   \email{rafael_moral@@yahoo.com.br}
##'
##' @param boot.scheme character bootstrap scheme passed from \code{lcc};
##'   defaults to \code{"np_case"}.
##' @param ci.method character CI method passed from \code{lcc};
##'   defaults to \code{"normal"}.
##' @keywords internal
lccInternal <- function(model, q_f, q_r, tk, interaction, covar,
                        pdmat, diffbeta, time_lcc, ci,
                        alpha, nboot, labels, var.class, weights.form,
                        show.warnings, components, lme.control,
                        method.init, numCore,
                        boot.scheme = "np_case",
                        ci.method   = "normal") {
  #-------------------------------------------------------------------
  # Time grids for prediction and plotting
  # tk.plot: prediction grid (possibly from time_lcc)
  # tk.plot2: original time grid (always from model)
  #-------------------------------------------------------------------
  tk.plot  <- tk
  tk.plot2 <- tk
  
  if (!is.null(time_lcc)) {
    if (is.null(time_lcc$time)) {
      tk.plot <- time_lcc(
        time = model$data$time,
        from = time_lcc$from,
        to   = time_lcc$to,
        n    = time_lcc$n
      )
    } else {
      tk.plot <- time_lcc(
        time = time_lcc$time,
        from = time_lcc$from,
        to   = time_lcc$to,
        n    = time_lcc$n
      )
    }
  }
  
  # Number of method comparisons (pairs vs reference)
  ldb <- length(diffbeta)
  
  # Variance-structure information (used to choose n.delta)
  varcomp   <- summary(model)
  varStruct <- varcomp$modelStruct$varStruct
  nd        <- length(varStruct)
  
  # Preallocate objects
  rho             <- NULL
  rho.ret         <- NULL
  rho.pearson     <- NULL
  rho.pearson.ret <- NULL
  Cb              <- NULL
  Cb.ret          <- NULL
  CI              <- NULL
  
  #===================================================================
  # Case 1: no bootstrap confidence intervals
  #===================================================================
  if (!ci) {
    # Precompute longitudinal kernel once for this model and tk.plot
    G       <- nlme::getVarCov(model)
    q_r_eff <- nrow(G) - 1L
    pre     <- .precompute_longitudinal(
      model = model, tk = tk.plot, q_f = q_f, q_r = q_r_eff
    )
    
    #---------------------------------------------------------------
    # ldb == 1 (two methods)
    #---------------------------------------------------------------
    if (ldb == 1L) {
      beta1 <- as.numeric(diffbeta[[1L]])
      
      # LCC
      rho_list <- .compute_LCC(pre, diffbeta = beta1)
      rho      <- rho_list[[1L]]  # for ldb == 1, n.delta is always 1
      
      if (components) {
        # LPC
        rho_pearson_list <- .compute_LPC(pre)
        rho.pearson      <- rho_pearson_list[[1L]]
        
        # Accuracy component (Cb)
        LA_list <- .compute_LA(pre, diffbeta = beta1)
        Cb      <- LA_list[[1L]]
      }
      
      summary.lcc <- lccSummary(
        model       = model,
        q_f         = q_f,
        diffbeta    = diffbeta,
        tk          = tk,
        tk.plot     = tk.plot,
        tk.plot2    = tk.plot2,
        rho         = rho,
        ENV.LCC     = NULL,
        rho.pearson = rho.pearson,
        ENV.LPC     = NULL,
        Cb          = Cb,
        ENV.Cb      = NULL,
        ldb         = ldb,
        ci          = FALSE,
        components  = components
      )
      
    } else {
      #-------------------------------------------------------------
      # ldb > 1 (more than two methods)
      #   nd <= 1: n.delta = 1 for all pairs
      #   nd  > 1: n.delta = i for i-th pair
      #-------------------------------------------------------------
      rho_list <- vector("list", ldb)
      
      for (i in seq_len(ldb)) {
        beta_i  <- as.numeric(diffbeta[[i]])
        rho_all <- .compute_LCC(pre, diffbeta = beta_i)
        
        if (length(rho_all) == 1L || sum(is.na(rho_all[[2L]])) != 0) {
          # fallback as in lccWrapper()
          rho_list[[i]] <- rho_all[[1L]]
        } else {
          n_delta <- if (nd <= 1L) 1L else i
          rho_list[[i]] <- rho_all[[n_delta]]
        }
      }
      rho.ret <- as.data.frame(do.call(cbind, rho_list))
      
      if (components) {
        # LPC: computed once, then indexed by n.delta
        rho_pearson_all <- .compute_LPC(pre)
        rho.pearson.ret <- vector("list", ldb)
        Cb_list         <- vector("list", ldb)
        
        for (i in seq_len(ldb)) {
          n_delta <- if (nd <= 1L) 1L else i
          
          # LPC
          rho.pearson.ret[[i]] <- rho_pearson_all[[n_delta]]
          
          # Cb
          beta_i  <- as.numeric(diffbeta[[i]])
          LA_list <- .compute_LA(pre, diffbeta = beta_i)
          Cb_list[[i]] <- LA_list[[n_delta]]
        }
        
        rho.pearson.ret <- as.data.frame(do.call(cbind, rho.pearson.ret))
        Cb.ret          <- as.data.frame(do.call(cbind, Cb_list))
      }
      
      summary.lcc <- lccSummary(
        model       = model,
        q_f         = q_f,
        diffbeta    = diffbeta,
        tk          = tk,
        tk.plot     = tk.plot,
        tk.plot2    = tk.plot2,
        rho         = rho.ret,
        ENV.LCC     = NULL,
        rho.pearson = rho.pearson.ret,
        ENV.LPC     = NULL,
        Cb          = Cb.ret,
        ENV.Cb      = NULL,
        ldb         = ldb,
        ci          = FALSE,
        components  = components
      )
    }
    
  } else {
    #=================================================================
    # Case 2: bootstrap confidence intervals
    #=================================================================
    CI <- ciBuilder(
      model        = model,
      nboot        = nboot,
      q_f          = q_f,
      q_r          = q_r,
      interaction  = interaction,
      covar        = covar,
      pdmat        = pdmat,
      var.class    = var.class,
      weights.form = weights.form,
      show.warnings = show.warnings,
      tk           = tk.plot,
      diffbeta     = diffbeta,
      ldb          = ldb,
      tk.plot      = tk.plot,
      tk.plot2     = tk.plot2,
      ci           = TRUE,
      alpha        = alpha,
      components   = components,
      lme.control  = lme.control,
      method.init  = method.init,
      numCore      = numCore,
      boot.scheme  = boot.scheme,
      ci.method    = ci.method
    )
    
    summary.lcc <- lccSummary(
      model       = model,
      q_f         = q_f,
      diffbeta    = diffbeta,
      tk          = tk,
      tk.plot     = tk.plot,
      tk.plot2    = tk.plot2,
      rho         = CI$rho,
      ENV.LCC     = CI$ENV.LCC,
      rho.pearson = CI$LPC,
      ENV.LPC     = CI$ENV.LPC,
      Cb          = CI$Cb,
      ENV.Cb      = CI$ENV.Cb,
      ldb         = ldb,
      ci          = TRUE,
      components  = components
    )
  }
  
  #===================================================================
  # Build internal object used by lcc() and lccPlot()
  #===================================================================
  internal_lcc <- list(
    "Summary.lcc" = summary.lcc,
    "tk.plot"     = tk.plot,
    "tk.plot2"    = tk.plot2,
    "ldb"         = ldb,
    "ci"          = ci,
    "components"  = components,
    "nd"          = nd,
    "qf"          = q_f,
    "qr"          = q_r,
    "interaction" = interaction,
    "covar"       = covar
  )
  
  if (ldb == 1L) {
    if (!ci) {
      internal_lcc$rho <- rho
      if (components) {
        internal_lcc$rho.pearson <- rho.pearson
        internal_lcc$Cb          <- Cb
      }
    } else {
      internal_lcc$rho      <- CI$rho
      internal_lcc$ENV.LCC  <- CI$ENV.LCC
      internal_lcc$ENV.LPC  <- CI$ENV.LPC
      internal_lcc$ENV.LA   <- CI$ENV.Cb
      internal_lcc$alpha    <- alpha
      internal_lcc$nboot    <- nboot
      if (components) {
        internal_lcc$rho.pearson <- CI$LPC
        internal_lcc$Cb          <- CI$Cb
      }
    }
  } else {
    if (!ci) {
      internal_lcc$rho <- rho.ret
      if (components) {
        internal_lcc$rho.pearson <- rho.pearson.ret
        internal_lcc$Cb          <- Cb.ret
      }
    } else {
      internal_lcc$rho      <- CI$rho
      internal_lcc$ENV.LCC  <- CI$ENV.LCC
      internal_lcc$ENV.LPC  <- CI$ENV.LPC
      internal_lcc$ENV.LA   <- CI$ENV.Cb
      internal_lcc$alpha    <- alpha
      internal_lcc$nboot    <- nboot
      if (components) {
        internal_lcc$rho.pearson <- CI$LPC
        internal_lcc$Cb          <- CI$Cb
      }
    }
  }
  
  invisible(internal_lcc)
}

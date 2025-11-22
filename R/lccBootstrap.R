########################################################################
# Package: lcc                                                         #
#                                                                      #
# File: lccBootstrap.R                                                 #
#                                                                      #
# Contains: dataBootstrap, bootstrapSamples, lccBootstrap,             #
# lpcBootstrap, laBootstrap                                            #
#                                                                      #
# Written by Thiago de Paula Oliveira                                  #
# copyright (c) 2017-18, Thiago P. Oliveira                            #
#                                                                      #
# First version: 11/10/2017                                            #
# Last update: 29/07/2019                                              #
# License: GNU General Public License version 2 (June, 1991) or later  #
#                                                                      #
########################################################################

##' @title Internal functions to estimate fixed effects and variance
##'   components.
##'
##' @description This is an internally called functions used to estimate
##'   fixed effects and variance components for each bootstrap sample.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @importFrom nlme fixef
##'
##' @keywords internal
bootstrapSamples <- function(nboot, model, q_f, q_r, interaction, covar,
                             var.class, pdmat, weights.form, show.warnings,
                             tk, diffbeta, ldb, components,
                             lme.control, method.init, numCore) {
  ## Pre-allocate
  LCC_Boot <- vector("list", nboot)
  LPC_Boot <- if (components) vector("list", nboot) else NULL
  Cb_Boot  <- if (components) vector("list", nboot) else NULL
  
  warnings <- 0L
  
  Data  <- model$data
  Data2 <- split(Data, Data$subject)
  
  ## Subject-level bootstrap
  split_by_subject <- function() {
    id <- sample(names(Data2), length(Data2), replace = TRUE)
    do.call("rbind", Data2[id])
  }
  
  ## Fixed-effects pattern (same as in lccInternal)
  base_fx  <- names(nlme::fixef(model))
  base_lev <- levels(Data$method)
  pat <- vector("list", ldb)
  for (i in seq_len(ldb)) {
    nams     <- base_fx[grepl(paste0(base_lev[i + 1L], "|poly"), base_fx)]
    pat[[i]] <- match(nams, base_fx)
  }
  compute_betas <- function(fx) {
    lapply(pat, function(idx) fx[idx])
  }
  
  ## Variance-structure / n_delta
  varStruct <- model$modelStruct$varStruct
  nd_vs <- if (is.null(varStruct)) 1L else
    length(nlme::coef(varStruct, unconstrained = FALSE))
  use_delta_by_level <- nd_vs > 1L
  
  ## One bootstrap iteration
  one_bootstrap <- function(i) {
    Data_boot <- split_by_subject()
    rownames(Data_boot) <- NULL
    
    fit <- lccModel(
      dataset      = Data_boot,
      resp         = "resp",
      subject      = "subject",
      pdmat        = pdmat,
      method       = "method",
      time         = "time",
      qf           = q_f,
      qr           = q_r,
      interaction  = interaction,
      covar        = covar,
      gs           = NULL,
      var.class    = var.class,
      weights.form = weights.form,
      lme.control  = lme.control,
      method.init  = method.init
    )
    
    ## If fit failed, fall back to original model
    if (fit$wcount == 1L) {
      boot_mod <- model
    } else {
      boot_mod <- fit$model
    }
    
    ## Precompute kernel for this bootstrap model
    G_i   <- nlme::getVarCov(boot_mod)
    q_r_i <- nrow(G_i) - 1L
    
    pre <- .precompute_longitudinal(
      model = boot_mod,
      tk    = tk,
      q_f   = q_f,
      q_r   = q_r_i
    )
    
    fx_boot       <- nlme::fixef(boot_mod)
    diffbeta_boot <- compute_betas(fx_boot)
    
    ## ----- LCC -----
    if (ldb == 1L) {
      rho_all <- .compute_LCC(
        pre      = pre,
        diffbeta = as.numeric(diffbeta_boot[[1L]])
      )
      rho_boot <- rho_all[[1L]]
    } else {
      rho_list <- vector("list", ldb)
      for (j in seq_len(ldb)) {
        n_delta   <- if (use_delta_by_level) j else 1L
        rho_all_j <- .compute_LCC(
          pre      = pre,
          diffbeta = as.numeric(diffbeta_boot[[j]])
        )
        if (length(rho_all_j) == 1L || sum(is.na(rho_all_j[[2L]])) != 0L) {
          rho_list[[j]] <- rho_all_j[[1L]]
        } else {
          rho_list[[j]] <- rho_all_j[[n_delta]]
        }
      }
      rho_boot <- rho_list
    }
    
    ## ----- Components -----
    rho.pearson_boot <- NULL
    Cb_boot          <- NULL
    
    if (components) {
      rho_pearson_all <- .compute_LPC(pre)
      
      if (ldb == 1L) {
        rho.pearson_boot <- rho_pearson_all[[1L]]
      } else {
        rhoP_list <- vector("list", ldb)
        for (j in seq_len(ldb)) {
          n_delta <- if (use_delta_by_level) j else 1L
          rhoP_list[[j]] <- rho_pearson_all[[n_delta]]
        }
        rho.pearson_boot <- rhoP_list
      }
      
      if (ldb == 1L) {
        LA_all <- .compute_LA(
          pre      = pre,
          diffbeta = as.numeric(diffbeta_boot[[1L]])
        )
        Cb_boot <- LA_all[[1L]]
      } else {
        Cb_list <- vector("list", ldb)
        for (j in seq_len(ldb)) {
          n_delta  <- if (use_delta_by_level) j else 1L
          LA_all_j <- .compute_LA(
            pre      = pre,
            diffbeta = as.numeric(diffbeta_boot[[j]])
          )
          Cb_list[[j]] <- LA_all_j[[n_delta]]
        }
        Cb_boot <- Cb_list
      }
    }
    
    list(
      LCC    = rho_boot,
      LPC    = rho.pearson_boot,
      Cb     = Cb_boot,
      wcount = fit$wcount
    )
  }
  
  ## -------------------------------------------------------------------
  ## Serial vs parallel execution + progress bar
  ## -------------------------------------------------------------------
  if (numCore <= 1L) {
    ## Serial
    pb <- utils::txtProgressBar(min = 0, max = nboot, style = 3)
    on.exit(close(pb), add = TRUE)
    
    for (i in seq_len(nboot)) {
      res <- one_bootstrap(i)
      
      LCC_Boot[[i]] <- res$LCC
      if (isTRUE(components)) {
        LPC_Boot[[i]] <- res$LPC
        Cb_Boot[[i]]  <- res$Cb
      }
      warnings <- warnings + as.integer(res$wcount)
      
      utils::setTxtProgressBar(pb, i)
    }
    
  } else {
    ## Parallel with doSNOW + foreach; keep existing progress behaviour
    cl <- parallel::makeCluster(numCore)
    doSNOW::registerDoSNOW(cl)
    
    pb <- utils::txtProgressBar(min = 0, max = nboot, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    results <- foreach::foreach(
      i = seq_len(nboot),
      .options.snow = opts,
      .packages = c("nlme", "lcc")
    ) %dorng% {
      one_bootstrap(i)
    }
    
    close(pb)
    parallel::stopCluster(cl)
    
    for (i in seq_along(results)) {
      res <- results[[i]]
      LCC_Boot[[i]] <- res$LCC
      if (isTRUE(components)) {
        LPC_Boot[[i]] <- res$LPC
        Cb_Boot[[i]]  <- res$Cb
      }
      warnings <- warnings + as.integer(res$wcount)
    }
  }
  
  if (show.warnings) {
    cat("\n  Convergence error in", warnings, "out of",
        nboot, "bootstrap samples.\n")
  }
  
  out <- list(
    LCC_Boot = LCC_Boot,
    LPC_Boot = LPC_Boot,
    Cb_Boot  = Cb_Boot
  )
  class(out) <- "lcc.bootstrap"
  out
}

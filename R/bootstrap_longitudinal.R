########################################################################
# Package: lcc                                                         #
#                                                                      #
# File: bootstrap_longitudinal.R                                       #
#                                                                      #
# Contains: bootstrapSamples, lccBootstrap,                            #
#            lpcBootstrap, laBootstrap                                 #
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
##' @importFrom foreach foreach %dopar%
##' @importFrom doRNG %dorng%
##' @importFrom doSNOW registerDoSNOW
##' @importFrom parallel makeCluster stopCluster
##' @keywords internal
bootstrapSamples <- function(nboot, model, q_f, q_r, interaction, covar,
                             var.class, pdmat, weights.form, show.warnings,
                             tk, diffbeta, ldb, components,
                             lme.control, method.init, numCore) {
  ## Pre-allocate
  n_tk     <- length(tk)
  LCC_Boot <- if (ldb == 1L) matrix(NA_real_, nrow = n_tk, ncol = nboot) else vector("list", nboot)
  LPC_Boot <- if (components) {
    if (ldb == 1L) matrix(NA_real_, nrow = n_tk, ncol = nboot) else vector("list", nboot)
  } else NULL
  Cb_Boot  <- if (components) {
    if (ldb == 1L) matrix(NA_real_, nrow = n_tk, ncol = nboot) else vector("list", nboot)
  } else NULL
  
  warnings <- 0L
  fail_indices <- integer(0L)
  
  Data           <- model$data
  subj_idx_list  <- split(seq_len(nrow(Data)), Data$subject)
  subj_names     <- names(subj_idx_list)
  n_subj         <- length(subj_idx_list)
  
  ## Subject-level bootstrap via indices to reduce copying
  split_by_subject <- function() {
    id  <- sample.int(n_subj, n_subj, replace = TRUE)
    idx <- unlist(subj_idx_list[id], use.names = FALSE)
    Data[idx, , drop = FALSE]
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
    length(coef(varStruct, unconstrained = FALSE))
  use_delta_by_level <- nd_vs > 1L
  
  ## Precompute polynomial bases for tk/q_f/q_r to reuse in .precompute_longitudinal
  basis_template <- list(
    tk   = tk,
    q_f  = q_f,
    q_r  = q_r,
    Tk_f = outer(tk, 0:q_f, `^`),
    Tk_r = outer(tk, 0:q_r, `^`)
  )
  
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
      q_r   = q_r_i,
      basis = if (q_r_i == q_r) basis_template else NULL
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
      
      if (ldb == 1L) {
        LCC_Boot[, i] <- res$LCC
        if (isTRUE(components)) {
          LPC_Boot[, i] <- res$LPC
          Cb_Boot[, i]  <- res$Cb
        }
      } else {
        LCC_Boot[[i]] <- res$LCC
        if (isTRUE(components)) {
          LPC_Boot[[i]] <- res$LPC
          Cb_Boot[[i]]  <- res$Cb
        }
      }
      if (res$wcount > 0L) {
        warnings <- warnings + as.integer(res$wcount)
        fail_indices <- c(fail_indices, i)
      }
      
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
      if (ldb == 1L) {
        LCC_Boot[, i] <- res$LCC
        if (isTRUE(components)) {
          LPC_Boot[, i] <- res$LPC
          Cb_Boot[, i]  <- res$Cb
        }
      } else {
        LCC_Boot[[i]] <- res$LCC
        if (isTRUE(components)) {
          LPC_Boot[[i]] <- res$LPC
          Cb_Boot[[i]]  <- res$Cb
        }
      }
      if (res$wcount > 0L) {
        warnings <- warnings + as.integer(res$wcount)
        fail_indices <- c(fail_indices, i)
      }
    }
  }
  
  if (show.warnings) {
    cat("\n  Convergence error in", warnings, "out of",
        nboot, "bootstrap samples.\n")
    if (length(fail_indices)) {
      cat("  Failed sample indices:", paste(fail_indices, collapse = ", "), "\n")
    }
  }
  
  out <- list(
    LCC_Boot = LCC_Boot,
    LPC_Boot = LPC_Boot,
    Cb_Boot  = Cb_Boot
  )
  class(out) <- "lcc.bootstrap"
  out
}
# lccBootstrap / lpcBootstrap / laBootstrap                           #
#######################################################################

##' @keywords internal
lccBootstrap <- function(model_boot, diff_boot, ldb, nboot, tk, q_f) {
  CCC_Boot <- vector("list", nboot)
  
  nd <- length(summary(model_boot[[1L]])$modelStruct$varStruct)
  use_delta_by_level <- nd > 1L
  
  for (i in seq_len(nboot)) {
    model_i <- model_boot[[i]]
    G_i     <- getVarCov(model_i)
    q_r_i   <- nrow(G_i) - 1L
    
    pre_i <- .precompute_longitudinal(model_i, tk, q_f = q_f, q_r = q_r_i)
    
    if (ldb == 1L) {
      CCC_Boot[[i]] <- {
        rho_list <- .compute_LCC(pre_i, diffbeta = as.numeric(diff_boot[[i]][[1L]]))
        # same selection rule as lccWrapper
        if (length(rho_list) == 1L || sum(is.na(rho_list[[2L]])) != 0) {
          rho_list[[1L]]
        } else {
          rho_list[[1L]]  # when ldb==1, n.delta is always 1
        }
      }
    } else {
      CCC_i <- vector("list", ldb)
      for (j in seq_len(ldb)) {
        n_delta <- if (use_delta_by_level) j else 1L
        rho_list <- .compute_LCC(pre_i, diffbeta = as.numeric(diff_boot[[i]][[j]]))
        if (length(rho_list) == 1L || sum(is.na(rho_list[[2L]])) != 0) {
          CCC_i[[j]] <- rho_list[[1L]]
        } else {
          CCC_i[[j]] <- rho_list[[n_delta]]
        }
      }
      CCC_Boot[[i]] <- CCC_i
    }
  }
  
  CCC_Boot
}

##' @keywords internal
lpcBootstrap <- function(model_boot, ldb, nboot, tk, q_f) {
  LPC_Boot <- vector("list", nboot)
  
  if (ldb == 1L) {
    for (i in seq_len(nboot)) {
      model_i <- model_boot[[i]]
      G_i     <- getVarCov(model_i)
      q_r_i   <- nrow(G_i) - 1L
      
      pre_i <- .precompute_longitudinal(model_i, tk, q_f = q_f, q_r = q_r_i)
      rho_pearson_list <- .compute_LPC(pre_i)
      
      LPC_Boot[[i]] <- rho_pearson_list[[1L]]
    }
    return(LPC_Boot)
  }
  
  # ldb > 1
  nd <- length(summary(model_boot[[1L]])$modelStruct$varStruct)
  use_delta_by_level <- nd > 1L
  
  for (i in seq_len(nboot)) {
    model_i <- model_boot[[i]]
    G_i     <- getVarCov(model_i)
    q_r_i   <- nrow(G_i) - 1L
    
    pre_i <- .precompute_longitudinal(model_i, tk, q_f = q_f, q_r = q_r_i)
    rho_pearson_list <- .compute_LPC(pre_i)
    
    LPC_i <- vector("list", ldb)
    for (j in seq_len(ldb)) {
      n_delta <- if (use_delta_by_level) j else 1L
      LPC_i[[j]] <- rho_pearson_list[[n_delta]]
    }
    LPC_Boot[[i]] <- LPC_i
  }
  
  LPC_Boot
}

##' @keywords internal
laBootstrap <- function(model_boot, diff_boot, ldb, nboot, tk, q_f) {
  Cb_Boot <- vector("list", nboot)
  
  if (ldb == 1L) {
    for (i in seq_len(nboot)) {
      model_i <- model_boot[[i]]
      G_i     <- getVarCov(model_i)
      q_r_i   <- nrow(G_i) - 1L
      
      pre_i <- .precompute_longitudinal(model_i, tk, q_f = q_f, q_r = q_r_i)
      LA_list <- .compute_LA(pre_i, diffbeta = as.numeric(diff_boot[[i]][[1L]]))
      
      Cb_Boot[[i]] <- LA_list[[1L]]
    }
    return(Cb_Boot)
  }
  
  # ldb > 1
  nd <- length(summary(model_boot[[1L]])$modelStruct$varStruct)
  use_delta_by_level <- nd > 1L
  
  for (i in seq_len(nboot)) {
    model_i <- model_boot[[i]]
    G_i     <- getVarCov(model_i)
    q_r_i   <- nrow(G_i) - 1L
    
    pre_i <- .precompute_longitudinal(model_i, tk, q_f = q_f, q_r = q_r_i)
    
    Cb_i <- vector("list", ldb)
    for (j in seq_len(ldb)) {
      n_delta <- if (use_delta_by_level) j else 1L
      LA_list <- .compute_LA(pre_i, diffbeta = as.numeric(diff_boot[[i]][[j]]))
      Cb_i[[j]] <- LA_list[[n_delta]]
    }
    Cb_Boot[[i]] <- Cb_i
  }
  
  Cb_Boot
}

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

##' @title Internal Functions to Generate Bootstrap Samples Based on
##'   Dataset.
##'
##' @description This is an internally called functions used to generate
##'   bootstrap samples.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira,
##'   \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @importFrom utils txtProgressBar setTxtProgressBar capture.output
##'
##' @importFrom foreach foreach %dopar%
##'
##' @importFrom doRNG %dorng%
##'
##' @importFrom doSNOW registerDoSNOW
##'
##' @importFrom parallel makeCluster stopCluster
##'
##' @keywords internal
dataBootstrap <- function(model) {
  data   <- model$data
  subj   <- data$subject
  n_subj <- length(unique(subj))
  
  # split once per call
  split_by_subj <- split(data, subj)
  
  # sample subject indices with replacement
  sample_idx <- sample.int(n_subj, n_subj, replace = TRUE)
  
  frames <- lapply(seq_along(sample_idx), function(i) {
    idx <- sample_idx[i]
    df  <- split_by_subj[[idx]]
    df$subject <- factor(i, levels = seq_len(n_subj))
    df
  })
  
  do.call(rbind, frames)
}


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
                             var.class, pdmat, weights.form,
                             show.warnings, lme.control, method.init,
                             numCore,
                             tk, ldb, components = TRUE) {
  
  #-------------------------------------------------------------------
  # Preallocate outputs (no lme objects stored)
  #-------------------------------------------------------------------
  LCC_Boot <- vector("list", nboot)
  LPC_Boot <- if (isTRUE(components)) vector("list", nboot) else NULL
  Cb_Boot  <- if (isTRUE(components)) vector("list", nboot) else NULL
  warnings <- 0L
  
  #-------------------------------------------------------------------
  # Subject-level split for bootstrap resampling
  #-------------------------------------------------------------------
  orig_data    <- model$data
  subj         <- orig_data$subject
  n_subj       <- length(unique(subj))
  split_by_subj <- split(orig_data, subj)
  subj_seq     <- seq_len(n_subj)
  
  bootstrap_dataset <- function() {
    sample_idx <- sample.int(n_subj, n_subj, replace = TRUE)
    frames <- lapply(seq_along(sample_idx), function(i) {
      idx <- sample_idx[i]
      df  <- split_by_subj[[idx]]
      df$subject <- factor(i, levels = subj_seq)
      df
    })
    do.call(rbind, frames)
  }
  
  #-------------------------------------------------------------------
  # Precompute pattern of fixed-effect indices per method (unchanged)
  #-------------------------------------------------------------------
  base_fx  <- nlme::fixef(model)
  base_lev <- levels(orig_data$method)
  
  x <- y <- NULL
  lev_lab_df <- unique(merge(rep("method", q_f), base_lev))
  lev_lab_df <- transform(lev_lab_df, newcol = paste(x, y, sep = ""))
  
  pat <- lapply(
    seq(2L, length(base_lev)),
    function(j) grep(lev_lab_df$newcol[j], names(base_fx))
  )
  
  compute_betas <- function(fx) {
    lapply(pat, function(idx) -fx[idx])
  }
  
  ldb_local <- length(pat)
  if (!missing(ldb) && !is.null(ldb) && ldb_local != ldb) {
    warning(
      "bootstrapSamples(): length(pat) (", ldb_local,
      ") does not match ldb (", ldb, "). Using length(pat)."
    )
  }
  
  #-------------------------------------------------------------------
  # Number of variance-structure parameters for delta vs deltal choice
  #-------------------------------------------------------------------
  varStruct <- summary(model)$modelStruct$varStruct
  nd        <- if (is.null(varStruct)) 0L else length(varStruct)
  use_delta_by_level <- nd > 1L
  
  #-------------------------------------------------------------------
  # Single bootstrap iteration: fit + precompute + LCC/LPC/LA
  #-------------------------------------------------------------------
  one_bootstrap <- function(i) {
    fit <- lccModel(
      dataset      = bootstrap_dataset(),
      resp         = "resp",
      subject      = "subject",
      covar        = covar,
      method       = "method",
      time         = "time",
      qf           = q_f,
      qr           = q_r,
      interaction  = interaction,
      pdmat        = pdmat,
      var.class    = var.class,
      weights.form = weights.form,
      lme.control  = lme.control,
      method.init  = method.init
    )
    
    if (fit$wcount == 1L) {
      boot_mod <- model
      if (show.warnings) {
        cat("\n  Estimation problem on bootstrap sample", i, "\n")
      }
    } else {
      boot_mod <- fit$model
    }
    
    # Precompute once for this bootstrap sample
    G_i   <- getVarCov(boot_mod)
    q_r_i <- nrow(G_i) - 1L
    
    pre_i <- .precompute_longitudinal(
      model = boot_mod,
      tk    = tk,
      q_f   = q_f,
      q_r   = q_r_i
    )
    
    fx    <- nlme::fixef(boot_mod)
    betas <- compute_betas(fx)
    
    #----------------------------
    # LCC (CCC) per bootstrap
    #----------------------------
    if (ldb_local == 1L) {
      rho_list <- .compute_LCC(pre_i, diffbeta = as.numeric(betas[[1L]]))
      CCC_i <- if (length(rho_list) == 1L ||
                   sum(is.na(rho_list[[2L]])) != 0) {
        rho_list[[1L]]
      } else {
        rho_list[[1L]]  # when ldb==1, n.delta is always 1
      }
    } else {
      CCC_i <- vector("list", ldb_local)
      for (j in seq_len(ldb_local)) {
        n_delta <- if (use_delta_by_level) j else 1L
        rho_list <- .compute_LCC(pre_i, diffbeta = as.numeric(betas[[j]]))
        if (length(rho_list) == 1L ||
            sum(is.na(rho_list[[2L]])) != 0) {
          CCC_i[[j]] <- rho_list[[1L]]
        } else {
          CCC_i[[j]] <- rho_list[[n_delta]]
        }
      }
    }
    
    #----------------------------
    # LPC + Cb (components)
    #----------------------------
    LPC_i <- NULL
    Cb_i  <- NULL
    
    if (isTRUE(components)) {
      # LPC (Pearson)
      rho_pearson_list <- .compute_LPC(pre_i)
      if (ldb_local == 1L) {
        LPC_i <- rho_pearson_list[[1L]]
      } else {
        LPC_i <- vector("list", ldb_local)
        for (j in seq_len(ldb_local)) {
          n_delta <- if (use_delta_by_level) j else 1L
          LPC_i[[j]] <- rho_pearson_list[[n_delta]]
        }
      }
      
      # Cb (accuracy)
      if (ldb_local == 1L) {
        LA_list <- .compute_LA(pre_i, diffbeta = as.numeric(betas[[1L]]))
        Cb_i    <- LA_list[[1L]]
      } else {
        Cb_i <- vector("list", ldb_local)
        for (j in seq_len(ldb_local)) {
          n_delta <- if (use_delta_by_level) j else 1L
          LA_list <- .compute_LA(pre_i, diffbeta = as.numeric(betas[[j]]))
          Cb_i[[j]] <- LA_list[[n_delta]]
        }
      }
    }
    
    list(
      LCC    = CCC_i,
      LPC    = LPC_i,
      Cb     = Cb_i,
      wcount = fit$wcount
    )
  }
  
  #-------------------------------------------------------------------
  # Serial vs parallel execution
  #-------------------------------------------------------------------
  if (numCore <= 1L) {
    
    for (i in seq_len(nboot)) {
      res <- one_bootstrap(i)
      LCC_Boot[[i]] <- res$LCC
      if (isTRUE(components)) {
        LPC_Boot[[i]] <- res$LPC
        Cb_Boot[[i]]  <- res$Cb
      }
      warnings <- warnings + as.integer(res$wcount)
    }
    
  } else {
    cl <- parallel::makeCluster(numCore, type = "SOCK")
    doSNOW::registerDoSNOW(cl)
    
    pb <- txtProgressBar(max = nboot, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    results <- foreach::foreach(
      i = seq_len(nboot),
      .options.snow = opts,
      .packages = c("nlme")
    ) %dorng% {
      one_bootstrap(i)
    }
    
    close(pb)
    parallel::stopCluster(cl)
    
    for (i in seq_len(nboot)) {
      LCC_Boot[[i]] <- results[[i]]$LCC
      if (isTRUE(components)) {
        LPC_Boot[[i]] <- results[[i]]$LPC
        Cb_Boot[[i]]  <- results[[i]]$Cb
      }
      warnings <- warnings + as.integer(results[[i]]$wcount)
    }
  }
  
  cat("\n  Convergence error in", warnings, "out of",
      nboot, "bootstrap samples.\n")
  
  out <- list(
    LCC_Boot = LCC_Boot,
    LPC_Boot = LPC_Boot,
    Cb_Boot  = Cb_Boot
  )
  class(out) <- "lcc.bootstrap"
  out
}

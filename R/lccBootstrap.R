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
                             numCore) {
  
  # Preallocate outputs
  Boot_model <- vector("list", nboot)
  Diff       <- vector("list", nboot)
  warnings   <- 0L
  
  # -------------------------------------------------------------------
  # Precompute subject-level split for bootstrap resampling
  # -------------------------------------------------------------------
  orig_data   <- model$data
  subj        <- orig_data$subject
  n_subj      <- length(unique(subj))
  split_by_subj <- split(orig_data, subj)
  subj_seq    <- seq_len(n_subj)
  
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
  
  # -------------------------------------------------------------------
  # Precompute pattern of fixed-effect indices per method
  # (these do not change across bootstrap samples)
  # -------------------------------------------------------------------
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
  
  # -------------------------------------------------------------------
  # Single bootstrap iteration (used in both serial and parallel code)
  # -------------------------------------------------------------------
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
    
    fx    <- nlme::fixef(boot_mod)
    betas <- compute_betas(fx)
    
    list(model = boot_mod, betas = betas, wcount = fit$wcount)
  }
  
  # -------------------------------------------------------------------
  # Without parallelization
  # -------------------------------------------------------------------
  if (numCore == 1L) {
    pb <- txtProgressBar(
      title = "Processing the bootstrap confidence intervals",
      style = 3, min = 0, max = nboot
    )
    
    for (i in seq_len(nboot)) {
      res <- one_bootstrap(i)
      Boot_model[[i]] <- res$model
      Diff[[i]]       <- res$betas
      warnings        <- warnings + as.integer(res$wcount)
      
      setTxtProgressBar(
        pb, i,
        label = paste(round(i / nboot * 100, 0), "% done")
      )
    }
    close(pb)
    
  } else {
    # -----------------------------------------------------------------
    # With parallelization
    # -----------------------------------------------------------------
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
    
    # unpack results
    for (i in seq_len(nboot)) {
      Boot_model[[i]] <- results[[i]]$model
      Diff[[i]]       <- results[[i]]$betas
      warnings        <- warnings + as.integer(results[[i]]$wcount)
    }
  }
  
  cat("\n  Convergence error in", warnings, "out of",
      nboot, "bootstrap samples.\n")
  
  out <- list(Boot_Model = Boot_model, Diffbetas = Diff)
  class(out) <- "lcc.bootstrap"
  out
}
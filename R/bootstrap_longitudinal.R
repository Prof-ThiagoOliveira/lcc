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
##' @param boot.scheme character indicating the bootstrap resampling
##'   scheme; defaults to "np_case" (subject-level case bootstrap). Other
##'   options:
##'   \itemize{
##'     \item "np_case": resample subjects with replacement; keep original responses.
##'     \item "np_case_resid_gr": case bootstrap; replace response by fitted values plus residuals sampled from the pooled residuals.
##'     \item "np_case_resid_ir": case bootstrap; replace response by fitted values plus residuals sampled within each subject.
##'     \item "np_re_resid_gr": resample subject-specific fitted trajectories (includes random effects); add residuals sampled from the pooled residuals.
##'     \item "np_re_resid_ir": resample subject-specific fitted trajectories; add residuals resampled within each subject.
##'     \item "sp_case_pr": semiparametric case bootstrap; resample subjects, use their fitted trajectories, then add Gaussian noise with variance equal to the estimated residual variance.
##'     \item "p_re_pr": fully parametric; simulate random effects from the estimated covariance and residuals from Gaussian noise, then generate responses via X beta + Z u + eps.
##'   }
#' @importFrom MASS mvrnorm
bootstrapSamples <- function(nboot, model, q_f, q_r, interaction, covar,
                             var.class, pdmat, weights.form, show.warnings,
                             tk, diffbeta, ldb, components,
                             lme.control, method.init, numCore,
                             boot.scheme = "np_case") {
  if (numCore > 1L && !requireNamespace("lcc", quietly = TRUE)) {
    warning("Package 'lcc' must be installed to run bootstrap in parallel; falling back to serial.")
    numCore <- 1L
  }
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
  split_by_subject <- function(return_idx = FALSE) {
    id  <- sample.int(n_subj, n_subj, replace = TRUE)
    idx <- unlist(subj_idx_list[id], use.names = FALSE)
    db  <- Data[idx, , drop = FALSE]
    rownames(db) <- NULL
    if (return_idx) list(data = db, idx = idx) else db
  }
  
  ## Fixed-effects pattern (align with original diffbeta)
  base_fx  <- names(nlme::fixef(model))
  base_lev <- levels(Data$method)
  pat <- vector("list", ldb)
  for (i in seq_len(ldb)) {
    if (!is.null(names(diffbeta[[i]]))) {
      pat[[i]] <- match(names(diffbeta[[i]]), base_fx)
    } else {
      nams     <- base_fx[grepl(paste0(base_lev[i + 1L], "|poly"), base_fx)]
      pat[[i]] <- match(nams, base_fx)
    }
  }
  compute_betas <- function(fx) {
    lapply(seq_along(pat), function(i) {
      vals <- fx[pat[[i]]]
      if (!is.null(names(diffbeta[[i]]))) {
        vals <- vals[names(diffbeta[[i]])]
      }
      if (length(vals) != length(diffbeta[[i]])) {
        vals <- vals[seq_len(min(length(vals), length(diffbeta[[i]])))]
        if (length(vals) < length(diffbeta[[i]])) {
          vals <- c(vals, rep(NA_real_, length(diffbeta[[i]]) - length(vals)))
        }
      }
      vals
    })
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

  fitted_orig <- stats::fitted(model)
  eps_hat     <- residuals(model)
  eps_by_subj <- split(eps_hat, Data$subject)
  ranef_mat   <- as.matrix(nlme::ranef(model))
  pred_level1 <- predict(model, level = 1)
  pred_level0 <- predict(model, level = 0)
  
  # Case bootstrap: resample whole subjects with replacement; keep original
  # response/trajectory (nonparametric case/cluster bootstrap)
  generate_np_case <- function() {
    split_by_subject(return_idx = FALSE)
  }
  
  # Case bootstrap + resampled residuals added to level-1 fitted values
  # global = TRUE  -> residuals pooled across subjects
  # global = FALSE -> residuals resampled within each subject
  generate_np_case_resid <- function(global = TRUE) {
    res <- split_by_subject(return_idx = TRUE)
    db  <- res$data
    idx <- res$idx
    
    fitted_db <- fitted_orig[idx]
    
    if (global) {
      res_star <- sample(eps_hat, length(idx), replace = TRUE)
    } else {
      subj_boot <- db$subject
      res_star  <- numeric(length(idx))
      for (s in unique(subj_boot)) {
        pos  <- which(subj_boot == s)
        pool <- eps_by_subj[[as.character(s)]]
        res_star[pos] <- sample(pool, length(pos), replace = TRUE)
      }
    }
    
    db$resp <- fitted_db + res_star
    rownames(db) <- NULL
    db
  }
  
  # Nonparametric random-effects + residual bootstrap:
  # resample subjects (and their BLUPs) then add resampled residuals
  generate_np_re_resid <- function(global = TRUE) {
    out_list <- vector("list", n_subj)
    for (k in seq_len(n_subj)) {
      sid  <- sample(levels(Data$subject), 1L, replace = TRUE)
      rows <- subj_idx_list[[sid]]

      fitted_sid <- pred_level1[rows]

      if (global) {
        res_sid <- sample(eps_hat, length(rows), replace = TRUE)
      } else {
        pool    <- eps_by_subj[[as.character(sid)]]
        res_sid <- sample(pool, length(rows), replace = TRUE)
      }

      tmp        <- Data[rows, , drop = FALSE]
      tmp$resp   <- fitted_sid + res_sid
      out_list[[k]] <- tmp
    }
    db <- do.call(rbind, out_list)
    rownames(db) <- NULL
    db
  }
  
  # Semi-parametric: case resampling + Gaussian residual noise around level-1 fits
  generate_sp_case_pr <- function() {
    res <- split_by_subject(return_idx = TRUE)
    db  <- res$data
    idx <- res$idx

    fitted_db <- pred_level1[idx]

    res_star <- stats::rnorm(nrow(db), mean = 0, sd = model$sigma)
    db$resp  <- fitted_db + res_star
    db
  }
  
  # Fully parametric: simulate random effects and residuals from fitted covariances
  generate_p_re_pr <- function() {
    out_list <- vector("list", n_subj)
    G_hat <- as.matrix(nlme::getVarCov(model))
    beta_hat <- nlme::fixef(model)
    for (k in seq_len(n_subj)) {
      sid <- sample(levels(Data$subject), 1L, replace = TRUE)
      rows <- subj_idx_list[[sid]]
      subj_dat <- Data[rows, , drop = FALSE]

      X_i <- as.matrix(subj_dat$fixed)
      Z_i <- if (!is.null(subj_dat$fmla.rand)) as.matrix(subj_dat$fmla.rand) else matrix(1, nrow = nrow(subj_dat), ncol = 1)

      u_i <- as.numeric(MASS::mvrnorm(1, mu = rep(0, ncol(Z_i)), Sigma = G_hat))
      eps_i <- stats::rnorm(nrow(subj_dat), mean = 0, sd = model$sigma)
      subj_dat$resp <- as.numeric(X_i %*% beta_hat + Z_i %*% u_i + eps_i)
      out_list[[k]] <- subj_dat
    }
    db <- do.call(rbind, out_list)
    rownames(db) <- NULL
    db
  }
  
  ## One bootstrap iteration
  one_bootstrap <- function(i) {
    Data_boot <- switch(
      boot.scheme,
      "np_case"          = generate_np_case(),
      "np_case_resid_gr" = generate_np_case_resid(global = TRUE),
      "np_case_resid_ir" = generate_np_case_resid(global = FALSE),
      "np_re_resid_gr"   = generate_np_re_resid(global = TRUE),
      "np_re_resid_ir"   = generate_np_re_resid(global = FALSE),
      "sp_case_pr"       = generate_sp_case_pr(),
      "p_re_pr"          = generate_p_re_pr(),
      stop("Unknown 'boot.scheme' in bootstrapSamples()", call. = FALSE)
    )
    
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
  if (ldb == 1L) {
    if (!is.null(LCC_Boot)) stopifnot(all(LCC_Boot >= -1, LCC_Boot <= 1, na.rm = TRUE))
    if (!is.null(LPC_Boot)) stopifnot(all(LPC_Boot >= -1, LPC_Boot <= 1, na.rm = TRUE))
    if (!is.null(Cb_Boot))  stopifnot(all(Cb_Boot  >=  0, Cb_Boot  <= 1, na.rm = TRUE))
  }
  class(out) <- "lcc.bootstrap"
  out
}

#######################################################################
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

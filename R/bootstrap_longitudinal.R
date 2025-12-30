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
##' @importFrom doRNG %dorng% registerDoRNG
##' @importFrom doSNOW registerDoSNOW
##' @importFrom parallel makeCluster stopCluster
##' @importFrom stats coef
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
#' @param rng_seed optional integer seed used to initialise the parallel
#'   random number stream via \code{doRNG}. When \code{NULL}, the current
#'   global RNG state is used.
#' @importFrom MASS mvrnorm
bootstrapSamples <- function(nboot, model, q_f, q_r, interaction, covar,
                             var.class, pdmat, weights.form, show.warnings,
                             tk, diffbeta, ldb, components,
                             lme.control, method.init, numCore,
                             boot.scheme = "np_case",
                             keep_models = FALSE,
                             rng_seed = NULL) {
  pkgload_available <- FALSE
  if (numCore > 1L) {
    has_lcc <- requireNamespace("lcc", quietly = TRUE)
    pkgload_available <- requireNamespace("pkgload", quietly = TRUE)
    if (!has_lcc && !pkgload_available) {
      warn_general("Parallel bootstrap requires either an installed 'lcc' package or the 'pkgload' helper; falling back to serial execution.")
      numCore <- 1L
    }
  }

  if (!is.null(rng_seed)) {
    rng_seed <- check_scalar_integer(rng_seed, arg = "rng_seed")
  }

  lcc_model_env <- new.env(parent = emptyenv())

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
  boot_state <- vector("list", nboot)
  
  Data           <- model$data
  subj_idx_list  <- split(seq_len(nrow(Data)), Data$subject)
  subj_names     <- names(subj_idx_list)
  n_subj         <- length(subj_idx_list)

  ## Helper: extract random-effects covariance G_hat from VarCorr.lme
  extract_G_hat <- function(model) {
    vc <- nlme::VarCorr(model)
    vc_mat <- as.matrix(vc)

    n_row <- nrow(vc_mat)
    if (n_row < 2L) {
      abort_internal("VarCorr(model) has no random-effects rows; cannot build G_hat.")
    }

    re_rows <- seq_len(n_row - 1L)  ## drop the residual row

    ## Prefer StdDev column; fall back to sqrt(Variance)
    sd_col <- intersect(colnames(vc_mat), c("StdDev", "Std.Dev", "Std.Dev."))
    if (length(sd_col) == 0L && "Variance" %in% colnames(vc_mat)) {
      sd_re <- sqrt(as.numeric(vc_mat[re_rows, "Variance"]))
    } else if (length(sd_col)) {
      sd_re <- as.numeric(vc_mat[re_rows, sd_col[1L]])
    } else {
      abort_internal("Could not locate StdDev/Variance columns in VarCorr(model).")
    }

    names(sd_re) <- rownames(vc_mat)[re_rows]
    if (!length(sd_re) || anyNA(sd_re)) {
      abort_internal("Could not extract finite random-effects standard deviations.")
    }

    corr_re <- attr(vc, "correlation")
    q <- length(sd_re)
    if (is.null(corr_re)) {
      Cmat <- diag(1, nrow = q, ncol = q)
    } else {
      corr_re <- tryCatch(as.matrix(corr_re), error = function(e) NULL)
      if (is.null(corr_re) || any(dim(corr_re) != c(q, q))) {
        Cmat <- diag(1, nrow = q, ncol = q)
      } else {
        Cmat <- corr_re
      }
    }
    G_hat <- diag(sd_re, nrow = q, ncol = q) %*% Cmat %*% diag(sd_re, nrow = q, ncol = q)
    rownames(G_hat) <- colnames(G_hat) <- names(sd_re)
    G_hat
  }

  G_hat <- extract_G_hat(model)

  safe_chol <- function(mat, jitter = 1e-8, max_iter = 5L) {
    if (!is.matrix(mat) || !nrow(mat) || !ncol(mat)) {
      return(NULL)
    }
    attempt <- tryCatch(chol(mat), error = function(e) NULL)
    if (!is.null(attempt)) {
      return(attempt)
    }
    dim_mat <- nrow(mat)
    for (k in seq_len(max_iter)) {
      step <- jitter * (10^(k - 1L))
      attempt <- tryCatch(
        chol(mat + diag(step, dim_mat)),
        error = function(e) NULL
      )
      if (!is.null(attempt)) {
        return(attempt)
      }
    }
    NULL
  }
  
  ## Subject-level bootstrap via indices to reduce copying
  split_by_subject <- function(return_idx = FALSE) {
    draw <- sample.int(n_subj, n_subj, replace = TRUE)
    out  <- vector("list", length(draw))
    idx_list <- vector("list", length(draw))

    for (k in seq_along(draw)) {
      rows <- subj_idx_list[[draw[k]]]
      tmp  <- Data[rows, , drop = FALSE]

      subj_chr <- as.character(tmp$subject)
      tmp$subject <- factor(
        paste0(subj_chr, "__boot", sprintf("%05d", k)),
        levels = NULL
      )

      out[[k]] <- tmp
      idx_list[[k]] <- rows
    }

    db <- do.call(rbind, out)
    rownames(db) <- NULL
    db$subject <- factor(as.character(db$subject))

    if (!return_idx) {
      return(db)
    }

    idx <- unlist(idx_list, use.names = FALSE)
    list(data = db, idx = idx)
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

  basis_cache <- new.env(parent = emptyenv())
  assign(as.character(q_r), basis_template, envir = basis_cache)

  get_basis <- function(q_r_val) {
    key <- as.character(q_r_val)
    if (exists(key, envir = basis_cache, inherits = FALSE)) {
      return(get(key, envir = basis_cache, inherits = FALSE))
    }
    basis_new <- list(
      tk   = tk,
      q_f  = q_f,
      q_r  = q_r_val,
      Tk_f = basis_template$Tk_f,
      Tk_r = outer(tk, 0:q_r_val, `^`)
    )
    assign(key, basis_new, envir = basis_cache)
    basis_new
  }

  fitted_orig <- stats::fitted(model)
  eps_hat     <- residuals(model)
  eps_by_subj <- split(eps_hat, Data$subject)
  ranef_mat   <- as.matrix(nlme::ranef(model))
  pred_level1 <- predict(model, level = 1)
  pred_level0 <- predict(model, level = 0)

  ## Align ranef rows with subject order used in subj_idx_list
  if (!is.null(rownames(ranef_mat))) {
    ranef_mat <- ranef_mat[subj_names, , drop = FALSE]
  }

  ## Fixed and random effects from the model
  beta_hat <- nlme::fixef(model)

  ## Precompute subject-specific design matrices and residual covariances R_i
  ## This will be used for the fully parametric bootstrap (p_re_pr)
  subj_RE_mats <- vector("list", n_subj)
  names(subj_RE_mats) <- subj_names

  for (s in seq_len(n_subj)) {
    rows     <- subj_idx_list[[s]]
    subj_dat <- Data[rows, , drop = FALSE]

    ## lccModel stored the fixed and random design matrices
    X_i <- as.matrix(subj_dat$fixed)
    Z_i <- if (!is.null(subj_dat$fmla.rand)) {
      as.matrix(subj_dat$fmla.rand)
    } else {
      ## Fallback: random intercept only
      matrix(1, nrow = nrow(subj_dat), ncol = ncol(G_hat))
    }

    ## Marginal covariance of y_i from nlme (includes random + residual)
    V_i_obj <- nlme::getVarCov(
      model,
      individual = subj_names[s],
      type       = "marginal"
    )

    ## nlme::getVarCov(..., type = "marginal") returns a 1x1 container whose
    ## first element is the n_i x n_i covariance matrix.
    if (is.list(V_i_obj) || inherits(V_i_obj, "VarCov")) {
      if (length(V_i_obj) != 1L) {
        abort_internal("Unexpected structure of getVarCov(model, type = 'marginal')")
      }
      V_i <- as.matrix(V_i_obj[[1L]])
    } else {
      V_i <- as.matrix(V_i_obj)
    }

    if (!is.numeric(V_i) || !is.matrix(V_i)) {
      abort_internal("V_i is not a numeric matrix in bootstrapSamples()")
    }
    if (nrow(V_i) != nrow(Z_i)) {
      abort_internal(
        "Dimension mismatch between V_i and Z_i in bootstrapSamples(): nrow(V_i) = {.val {nrow(V_i)}}, nrow(Z_i) = {.val {nrow(Z_i)}}"
      )
    }

    ## Residual covariance: R_i = V_i - Z_i G_hat Z_i'
    R_i <- V_i - Z_i %*% G_hat %*% t(Z_i)
    ## Enforce symmetry and small positive-definite correction
    R_i <- (R_i + t(R_i)) / 2
    eig <- eigen(R_i, symmetric = TRUE)
    eig$values[eig$values < 0] <- 1e-8
    L_Ri <- eig$vectors %*% diag(sqrt(eig$values))

    subj_RE_mats[[s]] <- list(X = X_i, Z = Z_i, L_R = L_Ri)
  }

  ## ------------------------------------------------------------
  ## Shrinkage corrections (Thai et al. style)
  ## ------------------------------------------------------------

  ## 3.1 Random effects: transform empirical BLUP covariance to G_hat
  ##    (multivariate extension via linear map A)
  ranef_corr <- ranef_mat
  if (nrow(ranef_mat) > 0L && ncol(ranef_mat) > 0L) {
    G_emp   <- stats::cov(ranef_mat)
    finite  <- is.matrix(G_emp) && all(is.finite(G_emp))
    U_emp   <- if (finite) safe_chol(G_emp) else NULL
    U_hat   <- safe_chol(G_hat)

    if (!is.null(U_emp) && !is.null(U_hat)) {
      A <- t(U_hat) %*% solve(t(U_emp))
      ranef_c <- scale(ranef_mat, center = TRUE, scale = FALSE)
      ranef_corr <- ranef_c %*% t(A)
    }
  }

  ## Build corrected level-1 predictions for the NP RE schemes
  pred_level1_corr <- numeric(length(pred_level1))
  for (s in seq_len(n_subj)) {
    rows <- subj_idx_list[[s]]

    X_i <- subj_RE_mats[[s]]$X
    Z_i <- subj_RE_mats[[s]]$Z
    u_i_corr <- as.numeric(ranef_corr[s, ])

    pred_level1_corr[rows] <-
      as.numeric(X_i %*% beta_hat + Z_i %*% u_i_corr)
  }

  ## 3.2 Residuals: global variance rescaling to match model sigma^2
  ##     (Thai et al.'s homoscedastic case; varStruct handled by g(delta) in LCC)
  sig2_hat <- model$sigma^2
  eps_c    <- eps_hat - mean(eps_hat)
  var_eps_emp <- stats::var(eps_c)

  if (is.finite(var_eps_emp) && var_eps_emp > 0) {
    A_eps   <- sqrt(sig2_hat / var_eps_emp)
    eps_corr <- eps_c * A_eps
  } else {
    eps_corr <- eps_hat
  }
  eps_by_subj_corr <- split(eps_corr, Data$subject)

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
      res_star <- sample(eps_corr, length(idx), replace = TRUE)
    } else {
      subj_boot  <- as.character(db$subject)
      subj_split <- split(seq_along(subj_boot), subj_boot)
      subj_orig  <- as.character(Data$subject[idx])
      res_star   <- numeric(length(idx))

      for (boot_id in names(subj_split)) {
        pos  <- subj_split[[boot_id]]
        orig <- subj_orig[pos][1L]
        pool <- eps_by_subj_corr[[orig]]
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
    subj_levels <- names(subj_idx_list)
    for (k in seq_len(n_subj)) {
      sid  <- sample(subj_levels, 1L, replace = TRUE)
      rows <- subj_idx_list[[sid]]

      fitted_sid <- pred_level1_corr[rows]

      if (global) {
        res_sid <- sample(eps_corr, length(rows), replace = TRUE)
      } else {
        pool    <- eps_by_subj_corr[[sid]]
        res_sid <- sample(pool, length(rows), replace = TRUE)
      }
      
      tmp        <- Data[rows, , drop = FALSE]
      subj_new   <- paste0(as.character(tmp$subject[1L]), "__boot", sprintf("%05d", k))
      tmp$subject <- factor(subj_new, levels = subj_new)
      tmp$resp   <- fitted_sid + res_sid
      out_list[[k]] <- tmp
    }
    db <- do.call(rbind, out_list)
    rownames(db) <- NULL
    db$subject <- factor(as.character(db$subject))
    db
  }
  
  # Semi-parametric: case resampling + Gaussian residual noise around level-1 fits
  # (assumes homoscedastic, independent residuals with variance sigma^2)
  generate_sp_case_pr <- function() {
    res <- split_by_subject(return_idx = TRUE)
    db  <- res$data
    idx <- res$idx

    fitted_db <- pred_level1[idx]

    res_star <- stats::rnorm(nrow(db), mean = 0, sd = model$sigma)
    db$resp  <- fitted_db + res_star
    db
  }
  
  # Fully parametric: simulate random effects and residuals from full fitted covariances
  generate_p_re_pr <- function() {
    out_list <- vector("list", n_subj)
    beta_hat <- nlme::fixef(model)

    for (s in seq_len(n_subj)) {
      rows     <- subj_idx_list[[s]]
      subj_dat <- Data[rows, , drop = FALSE]
      mats     <- subj_RE_mats[[s]]

      X_i <- mats$X
      Z_i <- mats$Z
      L_R <- mats$L_R

      ## u_i* ~ N(0, G_hat)
      u_i <- as.numeric(MASS::mvrnorm(
        n     = 1L,
        mu    = rep(0, ncol(G_hat)),
        Sigma = G_hat
      ))

      ## eps_i* ~ N(0, R_i_hat) via L_R %*% z
      z      <- stats::rnorm(nrow(subj_dat))
      eps_i  <- as.numeric(L_R %*% z)

      subj_dat$resp <- as.numeric(X_i %*% beta_hat + Z_i %*% u_i + eps_i)
      out_list[[s]] <- subj_dat
    }

    db <- do.call(rbind, out_list)
    rownames(db) <- NULL
    db
  }
  
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
      abort_input("Unknown 'boot.scheme' in bootstrapSamples()")
    )
    
    if (!exists("fn", envir = lcc_model_env, inherits = FALSE)) {
      assign(
        "fn",
        get("lccModel", envir = asNamespace("lcc"), inherits = FALSE),
        envir = lcc_model_env
      )
    }

    fit <- get("fn", envir = lcc_model_env)(
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
    summary_boot <- summary(boot_mod)
    G_i_info <- extract_random_effects_cov(boot_mod)
    q_r_i    <- max(0L, G_i_info$n_re - 1L)
    basis_i  <- get_basis(q_r_i)
    
    pre <- .precompute_longitudinal(
      model = boot_mod,
      tk    = tk,
      q_f   = q_f,
      q_r   = q_r_i,
      basis = basis_i,
      G_info = G_i_info,
      summary_obj = summary_boot
    )
    
    fx_boot       <- nlme::fixef(boot_mod)
    diffbeta_boot <- compute_betas(fx_boot)
    var_struct_boot <- summary_boot$modelStruct$varStruct
    state_record <- if (keep_models) {
      boot_mod
    } else {
      list(
        fixef       = fx_boot,
        ranef_cov   = G_i_info$G,
        sigma       = boot_mod$sigma,
        var_struct  = if (!is.null(var_struct_boot)) {
          list(
            class        = class(var_struct_boot)[1L],
            coefficients = stats::coef(var_struct_boot, unconstrained = FALSE)
          )
        } else NULL,
        wcount      = fit$wcount,
        message     = if (fit$wcount > 0L) fit$message else NULL
      )
    }
    
    ## ----- Metrics via helper loop -----
    rho_boot <- .bootstrap_metric_loop(
      pre                  = pre,
      ldb                  = ldb,
      compute_metric_from_pre = function(pre_obj, diffbeta_vec) {
        .compute_LCC(pre_obj, diffbeta = diffbeta_vec)
      },
      diffbeta             = diffbeta_boot,
      use_delta_by_level   = use_delta_by_level,
      needs_diff           = TRUE,
      selector             = .select_metric_lcc,
      label                = "LCC"
    )

    if (components) {
      rho.pearson_boot <- .bootstrap_metric_loop(
        pre                  = pre,
        ldb                  = ldb,
        compute_metric_from_pre = function(pre_obj, diffbeta_unused) {
          .compute_LPC(pre_obj)
        },
        use_delta_by_level   = use_delta_by_level,
        needs_diff           = FALSE,
        selector             = .select_metric_default,
        label                = "LPC"
      )

      Cb_boot <- .bootstrap_metric_loop(
        pre                  = pre,
        ldb                  = ldb,
        compute_metric_from_pre = function(pre_obj, diffbeta_vec) {
          .compute_LA(pre_obj, diffbeta = diffbeta_vec)
        },
        diffbeta             = diffbeta_boot,
        use_delta_by_level   = use_delta_by_level,
        needs_diff           = TRUE,
        selector             = .select_metric_default,
        label                = "LA"
      )
    } else {
      rho.pearson_boot <- NULL
      Cb_boot <- NULL
    }

    list(
      LCC    = rho_boot,
      LPC    = rho.pearson_boot,
      Cb     = Cb_boot,
      wcount = fit$wcount,
      state  = state_record
    )
  }
  
  ## -------------------------------------------------------------------
  ## Serial vs parallel execution + progress bar
  ## -------------------------------------------------------------------
  if (numCore <= 1L) {
    ## Serial
    pb <- utils::txtProgressBar(min = 0, max = nboot, style = 3)
    on.exit(close(pb), add = TRUE)
    
    if (!is.null(rng_seed)) {
      set.seed(rng_seed)
    }

    for (i in seq_len(nboot)) {
      res <- one_bootstrap(i)
      boot_state[[i]] <- res$state
      res$state <- NULL
      
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
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doSNOW::registerDoSNOW(cl)
    if (is.null(rng_seed)) {
      doRNG::registerDoRNG()
    } else {
      doRNG::registerDoRNG(rng_seed)
    }
    
    pb <- utils::txtProgressBar(min = 0, max = nboot, style = 3)
    on.exit(close(pb), add = TRUE)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    lib_paths <- .libPaths()
    parallel::clusterCall(
      cl,
      function(paths) {
        .libPaths(paths)
        NULL
      },
      lib_paths
    )

    pkg_path <- getNamespaceInfo(asNamespace("lcc"), "path")
    parallel::clusterCall(
      cl,
      function(path, can_use_pkgload) {
        load_dev <- FALSE
        if (can_use_pkgload && requireNamespace("pkgload", quietly = TRUE)) {
          loader <- get("load_all", envir = asNamespace("pkgload"))
          load_dev <- tryCatch({
            loader(
              path,
              compile = FALSE,
              helpers = FALSE,
              attach_testthat = FALSE,
              reset = FALSE,
              quiet = TRUE
            )
            TRUE
          }, error = function(e) FALSE)
        }
        if (!load_dev) {
          if (!requireNamespace("lcc", quietly = TRUE)) {
            stop("Failed to load the 'lcc' package on a worker.")
          }
        }
        NULL
      },
      pkg_path,
      pkgload_available
    )

    results <- foreach::foreach(
      i = seq_len(nboot),
      .options.snow = opts,
      .packages = c("nlme", "MASS", "lcc"),
      .export = c("abort_internal", "abort_input", "warn_general", "inform_general")
    ) %dorng% {
      one_bootstrap(i)
    }
    
    for (i in seq_along(results)) {
      res <- results[[i]]
      boot_state[[i]] <- res$state
      res$state <- NULL
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
    inform_general(
      "Convergence error in {.val {warnings}} out of {.val {nboot}} bootstrap samples."
    )
    if (length(fail_indices)) {
      inform_general(
        "Failed sample indices: {paste(fail_indices, collapse = ', ')}"
      )
    }
  }
  
  out <- list(
    LCC_Boot = LCC_Boot,
    LPC_Boot = LPC_Boot,
    Cb_Boot  = Cb_Boot,
    state    = boot_state
  )
  .check_range <- function(x, lower, upper, name) {
    if (!all(x >= lower & x <= upper, na.rm = TRUE)) {
      abort_internal("Bootstrap {name} values fell outside [{lower}, {upper}].")
    }
  }
  if (ldb == 1L) {
    if (!is.null(LCC_Boot)) .check_range(LCC_Boot, -1, 1, "LCC")
    if (!is.null(LPC_Boot)) .check_range(LPC_Boot, -1, 1, "LPC")
    if (!is.null(Cb_Boot))  .check_range(Cb_Boot,  0, 1, "Cb")
  }
  class(out) <- "lcc.bootstrap"
  out
}

.select_metric_default <- function(metric_stack, n_delta) {
  if (is.list(metric_stack)) {
    idx <- min(n_delta, length(metric_stack))
    return(metric_stack[[idx]])
  }
  metric_stack
}

.select_metric_lcc <- function(metric_stack, n_delta) {
  if (!is.list(metric_stack) || !length(metric_stack)) {
    abort_internal("LCC metric stack must be a non-empty list.")
  }
  fallback <- metric_stack[[1L]]
  if (length(metric_stack) < 2L) {
    return(fallback)
  }
  candidate_idx <- min(n_delta, length(metric_stack))
  candidate <- metric_stack[[candidate_idx]]
  second <- metric_stack[[2L]]
  if (anyNA(second)) {
    return(fallback)
  }
  candidate
}

.bootstrap_metric_loop <- function(pre, ldb, compute_metric_from_pre,
                                   diffbeta = NULL,
                                   use_delta_by_level = FALSE,
                                   needs_diff = FALSE,
                                   selector = .select_metric_default,
                                   label = "metric") {
  if (!is.function(compute_metric_from_pre)) {
    abort_internal(".bootstrap_metric_loop requires a function for {.arg compute_metric_from_pre}.")
  }
  if (!is.function(selector)) {
    abort_internal(".bootstrap_metric_loop requires a selector function.")
  }

  to_numeric <- function(x) {
    if (is.null(x)) return(NULL)
    as.numeric(x)
  }

  eval_metric <- function(diff_arg) {
    compute_metric_from_pre(pre, diff_arg)
  }

  if (ldb <= 1L) {
    diff_arg <- if (needs_diff) {
      if (is.null(diffbeta) || !length(diffbeta)) {
        abort_internal("{label}: missing diffbeta for single comparison.")
      }
      to_numeric(diffbeta[[1L]])
    } else NULL
    metric_stack <- eval_metric(diff_arg)
    return(selector(metric_stack, 1L))
  }

  out <- vector("list", ldb)

  if (!needs_diff) {
    metric_stack <- eval_metric(NULL)
    for (j in seq_len(ldb)) {
      n_delta <- if (use_delta_by_level) j else 1L
      out[[j]] <- selector(metric_stack, n_delta)
    }
    return(out)
  }

  if (is.null(diffbeta) || length(diffbeta) < ldb) {
    abort_internal("{label}: expected diffbeta list of length {.val {ldb}}.")
  }

  for (j in seq_len(ldb)) {
    diff_arg <- to_numeric(diffbeta[[j]])
    metric_stack <- eval_metric(diff_arg)
    n_delta <- if (use_delta_by_level) j else 1L
    out[[j]] <- selector(metric_stack, n_delta)
  }

  out
}

#######################################################################
# lccBootstrap / lpcBootstrap / laBootstrap                           #
#######################################################################

##' @keywords internal
lccBootstrap <- function(model_boot, diff_boot, ldb, nboot, tk, q_f) {
  CCC_Boot <- vector("list", nboot)

  for (i in seq_len(nboot)) {
    model_i <- model_boot[[i]]
    if (is.null(model_i)) {
      CCC_Boot[[i]] <- NULL
      next
    }

    summary_i <- summary(model_i)
    G_i_info  <- extract_random_effects_cov(model_i)
    q_r_i     <- max(0L, G_i_info$n_re - 1L)

    pre_i <- .precompute_longitudinal(model_i, tk, q_f = q_f, q_r = q_r_i,
                                      summary_obj = summary_i, G_info = G_i_info)

    CCC_Boot[[i]] <- .bootstrap_metric_loop(
      pre                  = pre_i,
      ldb                  = ldb,
      compute_metric_from_pre = function(pre_obj, diffbeta_vec) {
        .compute_LCC(pre_obj, diffbeta = diffbeta_vec)
      },
      diffbeta             = diff_boot[[i]],
      use_delta_by_level   = length(summary_i$modelStruct$varStruct) > 1L,
      needs_diff           = TRUE,
      selector             = .select_metric_lcc,
      label                = "LCC"
    )
  }

  CCC_Boot
}

##' @keywords internal
lpcBootstrap <- function(model_boot, ldb, nboot, tk, q_f) {
  LPC_Boot <- vector("list", nboot)

  for (i in seq_len(nboot)) {
    model_i <- model_boot[[i]]
    if (is.null(model_i)) {
      LPC_Boot[[i]] <- NULL
      next
    }

    summary_i <- summary(model_i)
    G_i_info  <- extract_random_effects_cov(model_i)
    q_r_i     <- max(0L, G_i_info$n_re - 1L)

    pre_i <- .precompute_longitudinal(model_i, tk, q_f = q_f, q_r = q_r_i,
                                      summary_obj = summary_i, G_info = G_i_info)

    LPC_Boot[[i]] <- .bootstrap_metric_loop(
      pre                  = pre_i,
      ldb                  = ldb,
      compute_metric_from_pre = function(pre_obj, diffbeta_unused) {
        .compute_LPC(pre_obj)
      },
      use_delta_by_level   = length(summary_i$modelStruct$varStruct) > 1L,
      needs_diff           = FALSE,
      selector             = .select_metric_default,
      label                = "LPC"
    )
  }

  LPC_Boot
}

##' @keywords internal
laBootstrap <- function(model_boot, diff_boot, ldb, nboot, tk, q_f) {
  Cb_Boot <- vector("list", nboot)

  for (i in seq_len(nboot)) {
    model_i <- model_boot[[i]]
    if (is.null(model_i)) {
      Cb_Boot[[i]] <- NULL
      next
    }

    summary_i <- summary(model_i)
    G_i_info  <- extract_random_effects_cov(model_i)
    q_r_i     <- max(0L, G_i_info$n_re - 1L)

    pre_i <- .precompute_longitudinal(model_i, tk, q_f = q_f, q_r = q_r_i,
                                      summary_obj = summary_i, G_info = G_i_info)

    Cb_Boot[[i]] <- .bootstrap_metric_loop(
      pre                  = pre_i,
      ldb                  = ldb,
      compute_metric_from_pre = function(pre_obj, diffbeta_vec) {
        .compute_LA(pre_obj, diffbeta = diffbeta_vec)
      },
      diffbeta             = diff_boot[[i]],
      use_delta_by_level   = length(summary_i$modelStruct$varStruct) > 1L,
      needs_diff           = TRUE,
      selector             = .select_metric_default,
      label                = "LA"
    )
  }

  Cb_Boot
}

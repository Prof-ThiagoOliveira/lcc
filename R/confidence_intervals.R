
#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: confidence_intervals.R                                        #
# Contains: lcc_intervals, ciBuilder, ciCompute and CI helpers        #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 29/07/2019                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# --------------------------------------------------------------------
# Shared CI transforms + helper used by both lcc_intervals() and ciCompute()
# --------------------------------------------------------------------

ZFisher     <- function(x) 0.5 * log((1 + x) / (1 - x))
ZFisher_inv <- function(x) (exp(2 * x) - 1) / (exp(2 * x) + 1)

Arcsin     <- function(x) asin(sqrt(x))
Arcsin_inv <- function(x) sin(x)^2

.clamp_r <- function(x, eps = 1e-8) {
  pmax(-1 + eps, pmin(1 - eps, x))
}

build_ci_metric <- function(boot_list, alpha,
                            transform, inv_transform,
                            percentile) {
  if (is.null(boot_list) || (is.list(boot_list) && !length(boot_list))) {
    return(matrix(NA_real_, nrow = 2L, ncol = 0L))
  }
  if (is.list(boot_list)) {
    boot_list <- boot_list[ vapply(boot_list, function(x) length(x) > 0, logical(1L)) ]
    if (!length(boot_list)) {
      return(matrix(NA_real_, nrow = 2L, ncol = 0L))
    }
  }
  res <- try(
    .build_ci_from_boot(
      boot_list     = boot_list,
      alpha         = alpha,
      transform     = if (!percentile) transform else NULL,
      inv_transform = if (!percentile) inv_transform else NULL,
      percentile    = percentile
    ),
    silent = TRUE
  )
  if (inherits(res, "try-error")) {
    # Fallback: return NA matrix with columns matching tk length if available
    cols <- 0L
    if (is.list(boot_list) && length(boot_list)) {
      cols <- length(boot_list[[1]])
    } else if (is.matrix(boot_list)) {
      cols <- ncol(boot_list)
    }
    return(matrix(NA_real_, nrow = 2L, ncol = cols))
  }
  res
}

##' @keywords internal
.align_ci_to_grid <- function(ci_mat, len_time, name = "CI") {
  if (is.null(ci_mat)) return(ci_mat)
  if (!is.matrix(ci_mat)) {
    abort_internal("{name} must be a matrix with 2 rows.")
  }
  if (nrow(ci_mat) != 2L && ncol(ci_mat) == 2L) {
    ci_mat <- t(ci_mat)
  }
  if (nrow(ci_mat) != 2L) {
    abort_internal(
      "{name} must be a 2-row matrix (lower/upper). Got {nrow(ci_mat)} rows."
    )
  }
  T_dim <- ncol(ci_mat)
  if (T_dim == len_time) {
    return(ci_mat)
  }
  if (T_dim > len_time) {
    return(ci_mat[, seq_len(len_time), drop = FALSE])
  }
  extra <- matrix(NA_real_, nrow = 2L, ncol = len_time - T_dim)
  cbind(ci_mat, extra)
}

##' @keywords internal
.align_ci_list_to_grid <- function(ci_list, len_time, name = "CI") {
  if (is.null(ci_list)) return(ci_list)
  if (!is.list(ci_list)) {
    abort_internal("{name} must be NULL or a list for ldb > 1.")
  }
  lapply(seq_along(ci_list), function(i) {
    .align_ci_to_grid(ci_list[[i]], len_time, sprintf("%s[%d]", name, i))
  })
}

##' @title Internal Functions to Compute the Non-Parametric Confidence
##'   Intervals for LCC.
##'
##' @description This is an internally called functions used to compute
##'   the non-parametric confidence intervals for LCC.
##'
##' @param boot.scheme character. Bootstrap scheme; defaults to
##'   \code{"np_case"} (non-parametric case resampling). Other options
##'   include residual-based and parametric variants.
##' @param ci.method character. Confidence interval method:
##'   \code{"normal"}, \code{"percentile"}, or \code{"bca"}.
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @importFrom stats quantile sd qnorm
##'
##' @keywords internal
##' @keywords internal
lcc_intervals <- function(rho, tk.plot, tk.plot2, ldb, model, ci,
                          LCC_Boot, alpha,
                          ci.method = "normal") {
  percentile <- identical(ci.method, "percentile")
  len_time <- if (is.null(dim(rho))) length(rho) else nrow(rho)
  
  if (ldb == 1L) {
    ## LCC_Boot is list over bootstrap samples; each element numeric over time
    ENV.LCC <- build_ci_metric(
      boot_list     = LCC_Boot,
      alpha         = alpha,
      transform     = ZFisher,
      inv_transform = ZFisher_inv,
      percentile    = percentile
    )
    ENV.LCC <- .align_ci_to_grid(ENV.LCC, len_time, name = "ENV.LCC")
    # Fallback: if CI is still all NA, recompute directly ignoring any NA replicates
    if (all(is.na(ENV.LCC)) && is.matrix(LCC_Boot) && ncol(LCC_Boot) > 0L) {
      lower <- apply(
        LCC_Boot, 1L, stats::quantile, probs = alpha / 2, na.rm = TRUE
      )
      upper <- apply(
        LCC_Boot, 1L, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE
      )
      ENV.LCC <- rbind(lower, upper)
      ENV.LCC <- .align_ci_to_grid(ENV.LCC, len_time, name = "ENV.LCC")
    }
  } else {
    ## LCC_Boot is list over bootstrap samples; each [[b]] is list over methods
    ENV.LCC <- vector("list", ldb)
    for (i in seq_len(ldb)) {
      boot_i <- lapply(LCC_Boot, function(x) if (!is.null(x)) x[[i]] else NULL)
      ENV.LCC[[i]] <- build_ci_metric(
        boot_list     = boot_i,
        alpha         = alpha,
        transform     = ZFisher,
        inv_transform = ZFisher_inv,
        percentile    = percentile
      )
    }
    ENV.LCC <- .align_ci_list_to_grid(ENV.LCC, len_time, name = "ENV.LCC")
  }
  
  CI.LCC <- list("rho" = rho, "ENV.LCC" = ENV.LCC)

  # Ensure CI matrices have columns matching the time grid length
  fill_if_empty <- function(env_obj, len) {
    if (is.null(env_obj)) return(env_obj)
    if (is.matrix(env_obj)) {
      if (!ncol(env_obj)) {
        return(matrix(NA_real_, nrow = 2L, ncol = len))
      }
      return(env_obj)
    }
    # list case
    lapply(env_obj, function(mat) {
      if (!ncol(mat)) {
        matrix(NA_real_, nrow = 2L, ncol = len)
      } else {
        mat
      }
    })
  }
  ENV.LCC <- fill_if_empty(ENV.LCC, length(tk.plot))
  CI.LCC$ENV.LCC <- ENV.LCC

  CI.LCC
}

##' @keywords Internal
.build_ci_from_boot <- function(boot_list, alpha,
                                transform = NULL, inv_transform = NULL,
                                percentile = FALSE) {
  ## boot_list: list over bootstrap replicates (vectors), or matrix with
  ## rows = time points and cols = replicates.

  if (is.matrix(boot_list)) {
    boot_mat <- boot_list
    if (!ncol(boot_mat)) {
      return(matrix(NA_real_, nrow = 2L, ncol = 0L))
    }
  } else {
    if (!is.list(boot_list)) {
      boot_list <- list(boot_list)
    }
    ## Drop NULL or all-NA replicates
    boot_list <- boot_list[!vapply(boot_list, is.null, logical(1L))]
    if (!length(boot_list)) {
      return(matrix(NA_real_, nrow = 2L, ncol = 0L))
    }
    boot_mat <- do.call(cbind, boot_list)
  }

  if (is.null(dim(boot_mat)) || !length(boot_mat) || any(dim(boot_mat) == 0)) {
    return(matrix(NA_real_, nrow = 2L, ncol = 0L))
  }
  
  if (percentile) {
    lower <- apply(
      boot_mat, 1L, stats::quantile, probs = alpha / 2, na.rm = TRUE
    )
    upper <- apply(
      boot_mat, 1L, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE
    )
    ci    <- rbind(lower, upper)
    return(ci)
  }
  
  ## Normal-approximation on a transformed scale
  if (!is.null(transform)) {
    boot_mat <- transform(.clamp_r(boot_mat))
  }
  
  se <- apply(boot_mat, 1L, sd, na.rm = TRUE)
  mu <- apply(boot_mat, 1L, mean, na.rm = TRUE)
  z  <- stats::qnorm(1 - alpha / 2)
  
  ci <- rbind(mu - z * se, mu + z * se)
  
  if (!is.null(inv_transform)) {
    ci <- inv_transform(ci)
  }
  
  ci
}
# --------------------------------------------------------------------
# ciBuilder / ciCompute
# --------------------------------------------------------------------

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
                      diffbeta, ldb, tk.plot, tk.plot2, ci,
                      boot.scheme = "np_case",
                      ci.method   = "normal",
                      alpha, components, lme.control, method.init,
                      numCore) {
  ci.method <- check_choice(ci.method, c("normal", "percentile", "bca"), arg = "ci.method")
  boot.scheme <- check_choice(
    boot.scheme,
    c(
      "np_case", "np_case_resid_gr", "np_case_resid_ir",
      "np_re_resid_gr", "np_re_resid_ir", "sp_case_pr", "p_re_pr"
    ),
    arg = "boot.scheme"
  )
  alpha <- check_scalar_numeric(alpha, arg = "alpha", lower = 0, upper = 1)
  nboot <- check_scalar_integer(nboot, arg = "nboot", lower = 1L)
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
    diffbeta     = diffbeta,
    ldb          = ldb,
    components   = components,
    boot.scheme  = boot.scheme
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
      LCC_Boot     = LCC_Boot,
      alpha        = alpha,
      ci.method    = ci.method
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
    LCC_Boot     = LCC_Boot,
    LPC_Boot     = LPC_Boot,
    Cb_Boot      = Cb_Boot,
    alpha        = alpha,
    ci.method    = ci.method,
    q_f          = q_f,
    q_r          = q_r,
    interaction  = interaction,
    covar        = covar,
    pdmat        = pdmat,
    var.class    = var.class,
    weights.form = weights.form,
    diffbeta     = diffbeta,
    components   = components,
    lme.control  = lme.control,
    method.init  = method.init
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
                      ci, LCC_Boot, LPC_Boot, Cb_Boot, alpha,
                      ci.method = "normal",
                      q_f = NULL, q_r = NULL, interaction = NULL, covar = NULL,
                      pdmat = NULL, var.class = NULL, weights.form = NULL,
                      diffbeta = NULL, components = FALSE,
                      lme.control = NULL, method.init = NULL) {
  percentile <- identical(ci.method, "percentile")
  
  .bca_interval <- function(theta_hat, boot, jack, alpha) {
    boot <- boot[is.finite(boot)]
    jack <- jack[is.finite(jack)]
    B <- length(boot)
    N <- length(jack)
    if (B < 10L || N < 3L) return(c(NA_real_, NA_real_))
    p0 <- mean(boot < theta_hat)
    z0 <- stats::qnorm(p0)
    jack_bar <- mean(jack)
    u <- jack_bar - jack
    a <- sum(u^3) / (6 * (sum(u^2)^(3/2)))
    z_alpha <- stats::qnorm(c(alpha / 2, 1 - alpha / 2))
    adj <- function(z) stats::pnorm(z0 + (z0 + z) / (1 - a * (z0 + z)))
    probs <- adj(z_alpha)
    q <- stats::quantile(boot, probs = probs, names = FALSE, type = 6)
    c(q[1L], q[2L])
  }
  
  jackknife_agreement <- function() {
    Data <- model$data
    subj_levels <- levels(Data$subject)
    N <- length(subj_levels)
    out <- vector("list", N)
    for (i in seq_len(N)) {
      keep <- Data$subject != subj_levels[i]
      Data_i <- Data[keep, , drop = FALSE]
      fit_i <- try(
        lccModel(
          dataset      = Data_i,
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
        ),
        silent = TRUE
      )
      if (inherits(fit_i, "try-error") || fit_i$wcount == 1L) {
        out[[i]] <- NULL
        next
      }
      mod_i <- fit_i$model
      G_i   <- nlme::getVarCov(mod_i)
      q_r_i <- nrow(G_i) - 1L
      pre_i <- .precompute_longitudinal(mod_i, tk.plot, q_f = q_f, q_r = q_r_i)
      nd    <- length(summary(mod_i)$modelStruct$varStruct)
      
      rho_i <- vector("list", length(diffbeta))
      for (j in seq_along(diffbeta)) {
        rho_all <- .compute_LCC(pre_i, diffbeta = as.numeric(diffbeta[[j]]))
        if (length(rho_all) == 1L || sum(is.na(rho_all[[2L]])) != 0) {
          rho_i[[j]] <- rho_all[[1L]]
        } else {
          n_delta <- if (nd <= 1L) 1L else j
          rho_i[[j]] <- rho_all[[n_delta]]
        }
      }
      if (length(rho_i) == 1L) rho_i <- rho_i[[1L]]
      
      rhoP_i <- NULL
      Cb_i   <- NULL
      if (components) {
        rhoP_all <- .compute_LPC(pre_i)
        Cb_all <- lapply(seq_along(diffbeta), function(j) {
          .compute_LA(pre_i, diffbeta = as.numeric(diffbeta[[j]]))
        })
        if (length(diffbeta) == 1L) {
          rhoP_i <- rhoP_all[[1L]]
          Cb_i   <- Cb_all[[1L]][[1L]]
        } else {
          rhoP_i <- vector("list", length(diffbeta))
          Cb_i   <- vector("list", length(diffbeta))
          for (j in seq_along(diffbeta)) {
            n_delta <- if (nd <= 1L) 1L else j
            rhoP_i[[j]] <- rhoP_all[[n_delta]]
            Cb_i[[j]]   <- Cb_all[[j]][[n_delta]]
          }
        }
      }
      out[[i]] <- list(rho = rho_i, rho.pearson = rhoP_i, Cb = Cb_i)
    }
    out
  }
  
  if (ci.method %in% c("normal", "percentile")) {
    if (ldb == 1L) {
      ENV.LCC <- build_ci_metric(LCC_Boot, alpha, ZFisher, ZFisher_inv, percentile)
      ENV.LCC <- .align_ci_to_grid(ENV.LCC, length(tk.plot), name = "ENV.LCC")
      if (all(is.na(ENV.LCC)) && is.matrix(LCC_Boot) && ncol(LCC_Boot) > 0L) {
        lower <- apply(
          LCC_Boot, 1L, stats::quantile, probs = alpha / 2, na.rm = TRUE
        )
        upper <- apply(
          LCC_Boot, 1L, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE
        )
        ENV.LCC <- rbind(lower, upper)
        ENV.LCC <- .align_ci_to_grid(ENV.LCC, length(tk.plot), name = "ENV.LCC")
      }
      ENV.LPC <- if (!is.null(LPC_Boot)) build_ci_metric(LPC_Boot, alpha, ZFisher, ZFisher_inv, percentile) else NULL
      if (!is.null(ENV.LPC)) {
        ENV.LPC <- .align_ci_to_grid(ENV.LPC, length(tk.plot), name = "ENV.LPC")
      }
      ENV.Cb  <- if (!is.null(Cb_Boot))  build_ci_metric(Cb_Boot,  alpha, Arcsin,  Arcsin_inv, percentile) else NULL
      if (!is.null(ENV.Cb)) {
        ENV.Cb <- .align_ci_to_grid(ENV.Cb, length(tk.plot), name = "ENV.Cb")
      }
    } else {
      ENV.LCC <- vector("list", ldb)
      ENV.LPC <- if (!is.null(LPC_Boot)) vector("list", ldb) else NULL
      ENV.Cb  <- if (!is.null(Cb_Boot))  vector("list", ldb) else NULL
      
      for (i in seq_len(ldb)) {
        LCC_i <- lapply(LCC_Boot, function(x) if (!is.null(x)) x[[i]] else NULL)
        ENV.LCC[[i]] <- build_ci_metric(LCC_i, alpha, ZFisher, ZFisher_inv, percentile)
        if (!is.null(LPC_Boot)) {
          LPC_i <- lapply(LPC_Boot, function(x) if (!is.null(x)) x[[i]] else NULL)
          ENV.LPC[[i]] <- build_ci_metric(LPC_i, alpha, ZFisher, ZFisher_inv, percentile)
        }
        if (!is.null(Cb_Boot)) {
          Cb_i  <- lapply(Cb_Boot,  function(x) if (!is.null(x)) x[[i]] else NULL)
          ENV.Cb[[i]]  <- build_ci_metric(Cb_i,  alpha, Arcsin,  Arcsin_inv, percentile)
        }
      }
      ENV.LCC <- .align_ci_list_to_grid(ENV.LCC, length(tk.plot), name = "ENV.LCC")
      if (!is.null(ENV.LPC)) {
        ENV.LPC <- .align_ci_list_to_grid(ENV.LPC, length(tk.plot), name = "ENV.LPC")
      }
      if (!is.null(ENV.Cb)) {
        ENV.Cb <- .align_ci_list_to_grid(ENV.Cb, length(tk.plot), name = "ENV.Cb")
      }
    }
  } else {  # BCa
    jack <- jackknife_agreement()
    if (ldb == 1L) {
      boot_mat <- if (is.matrix(LCC_Boot)) LCC_Boot else do.call(cbind, LCC_Boot)
      ENV.LCC <- matrix(NA_real_, nrow = 2L, ncol = nrow(boot_mat))
      for (t in seq_len(nrow(boot_mat))) {
        theta_hat <- rho[t]
        boot_t    <- boot_mat[t, ]
        jack_t    <- vapply(jack, function(j) if (!is.null(j)) j$rho[t] else NA_real_, numeric(1L))
        ENV.LCC[, t] <- .bca_interval(theta_hat, boot_t, jack_t, alpha)
      }
      ENV.LCC <- .align_ci_to_grid(ENV.LCC, length(tk.plot), name = "ENV.LCC")
      if (!is.null(LPC_Boot)) {
        boot_mat <- if (is.matrix(LPC_Boot)) LPC_Boot else do.call(cbind, LPC_Boot)
        ENV.LPC <- matrix(NA_real_, nrow = 2L, ncol = nrow(boot_mat))
        for (t in seq_len(nrow(boot_mat))) {
          theta_hat <- rho.pearson[t]
          boot_t    <- boot_mat[t, ]
          jack_t    <- vapply(jack, function(j) if (!is.null(j)) j$rho.pearson[t] else NA_real_, numeric(1L))
          ENV.LPC[, t] <- .bca_interval(theta_hat, boot_t, jack_t, alpha)
        }
      } else {
        ENV.LPC <- NULL
      }
      if (!is.null(Cb_Boot)) {
        boot_mat <- if (is.matrix(Cb_Boot)) Cb_Boot else do.call(cbind, Cb_Boot)
        ENV.Cb <- matrix(NA_real_, nrow = 2L, ncol = nrow(boot_mat))
        for (t in seq_len(nrow(boot_mat))) {
          theta_hat <- Cb[t]
          boot_t    <- boot_mat[t, ]
          jack_t    <- vapply(jack, function(j) if (!is.null(j)) j$Cb[t] else NA_real_, numeric(1L))
          ENV.Cb[, t] <- .bca_interval(theta_hat, boot_t, jack_t, alpha)
        }
      } else {
        ENV.Cb <- NULL
      }
    } else {
      ENV.LCC <- vector("list", ldb)
      ENV.LPC <- if (!is.null(LPC_Boot)) vector("list", ldb) else NULL
      ENV.Cb  <- if (!is.null(Cb_Boot))  vector("list", ldb) else NULL
      for (i in seq_len(ldb)) {
        boot_i <- lapply(LCC_Boot, function(x) if (!is.null(x)) x[[i]] else NULL)
        boot_mat <- do.call(cbind, boot_i)
        if (nrow(boot_mat) == 0L || ncol(boot_mat) == 0L) {
          ENV.LCC[[i]] <- matrix(NA_real_, nrow = 2L, ncol = length(tk.plot))
        } else {
          env_i <- matrix(NA_real_, nrow = 2L, ncol = nrow(boot_mat))
          for (t in seq_len(nrow(boot_mat))) {
            theta_hat <- rho[[i]][t]
            boot_t    <- boot_mat[t, ]
            jack_t    <- vapply(jack, function(j) if (!is.null(j)) j$rho[[i]][t] else NA_real_, numeric(1L))
            env_i[, t] <- .bca_interval(theta_hat, boot_t, jack_t, alpha)
          }
          ENV.LCC[[i]] <- env_i
        }
        
        if (!is.null(LPC_Boot)) {
          boot_i <- lapply(LPC_Boot, function(x) if (!is.null(x)) x[[i]] else NULL)
          boot_mat <- do.call(cbind, boot_i)
          if (nrow(boot_mat) == 0L || ncol(boot_mat) == 0L) {
            ENV.LPC[[i]] <- matrix(NA_real_, nrow = 2L, ncol = length(tk.plot))
          } else {
            env_lp <- matrix(NA_real_, nrow = 2L, ncol = nrow(boot_mat))
            for (t in seq_len(nrow(boot_mat))) {
              theta_hat <- rho.pearson[[i]][t]
              boot_t    <- boot_mat[t, ]
              jack_t    <- vapply(jack, function(j) if (!is.null(j)) j$rho.pearson[[i]][t] else NA_real_, numeric(1L))
              env_lp[, t] <- .bca_interval(theta_hat, boot_t, jack_t, alpha)
            }
            ENV.LPC[[i]] <- env_lp
          }
        }
        if (!is.null(Cb_Boot)) {
          boot_i <- lapply(Cb_Boot, function(x) if (!is.null(x)) x[[i]] else NULL)
          boot_mat <- do.call(cbind, boot_i)
          if (nrow(boot_mat) == 0L || ncol(boot_mat) == 0L) {
            ENV.Cb[[i]] <- matrix(NA_real_, nrow = 2L, ncol = length(tk.plot))
          } else {
            env_cb <- matrix(NA_real_, nrow = 2L, ncol = nrow(boot_mat))
            for (t in seq_len(nrow(boot_mat))) {
              theta_hat <- Cb[[i]][t]
              boot_t    <- boot_mat[t, ]
              jack_t    <- vapply(jack, function(j) if (!is.null(j)) j$Cb[[i]][t] else NA_real_, numeric(1L))
              env_cb[, t] <- .bca_interval(theta_hat, boot_t, jack_t, alpha)
            }
            ENV.Cb[[i]] <- env_cb
          }
        }
      }
      ENV.LCC <- .align_ci_list_to_grid(ENV.LCC, length(tk.plot), name = "ENV.LCC")
      if (!is.null(ENV.LPC)) {
        ENV.LPC <- .align_ci_list_to_grid(ENV.LPC, length(tk.plot), name = "ENV.LPC")
      }
      if (!is.null(ENV.Cb)) {
        ENV.Cb <- .align_ci_list_to_grid(ENV.Cb, length(tk.plot), name = "ENV.Cb")
      }
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

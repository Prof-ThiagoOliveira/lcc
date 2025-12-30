
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
##' @description Internally constructs the metric estimates, bootstrap
##'   summaries, and plotting auxiliaries consumed by downstream
##'   helpers. This routine is not exported.
##'
##' @param model fitted \code{nlme} object returned by \code{lccModel}.
##' @param q_f polynomial degree for fixed effects.
##' @param q_r polynomial degree requested for random effects.
##' @param interaction logical indicating whether a method-by-time
##'   interaction was estimated.
##' @param tk sorted vector of observed time points.
##' @param covar optional character vector of additional covariates in the
##'   mixed model.
##' @param pdmat variance-covariance structure identifier passed from
##'   \code{lcc}.
##' @param diffbeta list of fixed-effect contrasts between methods.
##' @param time_lcc optional list controlling the prediction grid.
##' @param ci logical indicating whether confidence intervals are requested.
##' @param boot.scheme bootstrap resampling strategy.
##' @param ci.method bootstrap confidence-interval construction method.
##' @param alpha significance level used for confidence intervals.
##' @param nboot number of bootstrap replicates.
##' @param labels internal labeling data produced in \code{lcc}.
##' @param var.class variance-function specification.
##' @param weights.form character describing the variance-function formula.
##' @param show.warnings logical; propagate convergence warnings from
##'   bootstrap fits.
##' @param components logical; include longitudinal Pearson correlation
##'   and accuracy components.
##' @param lme.control optional control list for \code{nlme::lme}.
##' @param method.init fitting method used for the initial model.
##' @param numCore number of cores allocated for bootstrap computation.
##' @param keep_models logical; retain bootstrap fit objects.
##' @param boot_seed optional integer seed for bootstrap reproducibility.
##'
##' @return Invisibly returns a list containing \code{Summary.lcc},
##'   plotting auxiliaries, cached metrics, and bootstrap state.
##'
##' @keywords internal
lccInternal <- function(model, q_f, q_r, interaction, tk,
                        covar, pdmat, diffbeta, time_lcc,
                        ci, boot.scheme, ci.method, alpha, nboot,
                        labels, var.class, weights.form, show.warnings,
                        components, lme.control, method.init, numCore,
                        keep_models, boot_seed) {

  Data <- model$data
  if (is.null(Data)) {
    abort_internal("Fitted model does not contain the original data.")
  }

  method_levels <- levels(Data$method)
  if (length(method_levels) < 2L) {
    abort_internal("At least two methods are required to compute longitudinal concordance.")
  }

  ldb <- length(diffbeta)
  if (ldb < 1L) {
    abort_internal("No method contrasts supplied in 'diffbeta'.")
  }

  tk_observed <- sort(unique(Data$time))
  if (!length(tk_observed)) {
    abort_internal("No time information available to build longitudinal summaries.")
  }

  if (missing(tk) || !length(tk)) {
    tk <- tk_observed
  }

  build_time_grid <- function(spec, observed) {
    if (is.null(spec)) {
      return(observed)
    }
    if (!is.list(spec)) {
      abort_input("'time_lcc' must be NULL or a list with elements 'time', 'from', 'to', and 'n'.")
    }
    candidate <- spec$time
    from <- if (!is.null(spec$from)) spec$from else min(observed)
    to   <- if (!is.null(spec$to))   spec$to   else max(observed)
    if (!is.finite(from) || !is.finite(to) || from > to) {
      abort_input("Invalid 'time_lcc' interval; ensure 'from' <= 'to' and both finite.")
    }
    n_raw <- spec$n
    if (is.null(n_raw) || !is.numeric(n_raw) || n_raw < 2) {
      n_raw <- max(length(observed), 2L)
    }
    n <- max(1L, as.integer(round(n_raw)))
    seq_time <- if (n > 1L || from == to) seq.int(from, to, length.out = n) else from
    base <- numeric(0)
    if (!is.null(candidate)) {
      base <- as.numeric(candidate)
      base <- base[is.finite(base)]
    }
    grid <- unique(c(base, seq_time, observed))
    sort(grid)
  }

  tk.plot  <- build_time_grid(time_lcc, tk_observed)
  tk.plot2 <- tk_observed

  comp_names <- names(diffbeta)
  if (is.null(comp_names) || length(comp_names) != ldb || anyNA(comp_names) || any(!nzchar(comp_names))) {
    comp_names <- NULL
  } else {
    comp_names <- as.character(comp_names)
  }

  if (is.null(comp_names)) {
    ref_method <- method_levels[1L]
    other_methods <- method_levels[seq_len(min(length(method_levels) - 1L, ldb)) + 1L]
    if (ldb == 1L) {
      if (length(other_methods) >= 1L) {
        comp_names <- paste(other_methods[1L], "vs", ref_method)
      } else {
        comp_names <- "Comparison_1"
      }
    } else {
      if (length(other_methods) < ldb) {
        comp_names <- sprintf("Comparison_%d", seq_len(ldb))
      } else {
        comp_names <- paste(other_methods, "vs", ref_method)
      }
    }
  }

  if (length(comp_names) != ldb) {
    comp_names <- sprintf("Comparison_%d", seq_len(ldb))
  }

  model_summary <- summary(model)
  nd <- length(model_summary$modelStruct$varStruct)
  resp_var <- stats::var(model$data$resp)
  degenerate_resp <- !is.finite(resp_var) || resp_var <= 0

  mask_metric_bundle <- function(metric_obj) {
    if (is.null(metric_obj) || is.null(metric_obj$estimate)) {
      return(metric_obj)
    }
    fill_len <- length(tk.plot)
    if (!is.null(metric_obj$grid) && !is.null(metric_obj$grid$primary)) {
      fill_len <- length(metric_obj$grid$primary)
    }
    est <- metric_obj$estimate
    if (!is.list(est)) {
      est <- list(as.numeric(est))
    }
    nm <- names(est)
    est <- lapply(est, function(x) rep(NA_real_, fill_len))
    if (!is.null(nm) && length(nm) == length(est)) {
      names(est) <- nm
    }
    metric_obj$estimate <- est
    if (!is.null(metric_obj$ci)) {
      metric_obj$ci <- lapply(est, function(x) {
        matrix(
          NA_real_,
          nrow = 2L,
          ncol = fill_len,
          dimnames = list(c("lower", "upper"), NULL)
        )
      })
    }
    metric_obj$bootstrap <- NULL
    metric_obj
  }

  as_comp_estimate <- function(values) {
    if (is.null(values)) return(NULL)
    if (is.data.frame(values) || is.matrix(values)) {
      values <- lapply(seq_len(ncol(values)), function(j) values[, j])
    } else if (!is.list(values) || is.numeric(values)) {
      values <- list(values)
    }
    out <- lapply(values, function(val) as.numeric(val))
    if (length(out) == ldb) {
      names(out) <- comp_names
    }
    out
  }

  if (!ci) {
    G_info  <- extract_random_effects_cov(model)
    q_r_eff <- max(0L, G_info$n_re - 1L)
    pre     <- .precompute_longitudinal(
      model = model,
      tk    = tk.plot,
      q_f   = q_f,
      q_r   = q_r_eff,
      G_info = G_info,
      summary_obj = model_summary
    )

    metrics <- list()

    if (ldb == 1L) {
      beta1 <- as.numeric(diffbeta[[1L]])
      rho_list <- .compute_LCC(pre, diffbeta = beta1)
      rho_est  <- rho_list[[1L]]

      metrics$lcc <- list(
        estimate  = as_comp_estimate(rho_est),
        bootstrap = NULL,
        ci        = NULL,
        grid      = list(primary = tk.plot, secondary = tk.plot2)
      )

      if (components) {
        rho_pearson_list <- .compute_LPC(pre)
        LA_list          <- .compute_LA(pre, diffbeta = beta1)

        metrics$lpc <- list(
          estimate  = as_comp_estimate(rho_pearson_list[[1L]]),
          bootstrap = NULL,
          ci        = NULL,
          grid      = list(primary = tk.plot, secondary = tk.plot2)
        )
        metrics$la <- list(
          estimate  = as_comp_estimate(LA_list[[1L]]),
          bootstrap = NULL,
          ci        = NULL,
          grid      = list(primary = tk.plot, secondary = tk.plot2)
        )
      }

    } else {
      rho_entries <- vector("list", ldb)
      for (i in seq_len(ldb)) {
        beta_i  <- as.numeric(diffbeta[[i]])
        rho_all <- .compute_LCC(pre, diffbeta = beta_i)
        if (length(rho_all) == 1L || sum(is.na(rho_all[[2L]])) != 0) {
          rho_entries[[i]] <- rho_all[[1L]]
        } else {
          n_delta <- if (nd <= 1L) 1L else i
          rho_entries[[i]] <- rho_all[[n_delta]]
        }
      }

      metrics$lcc <- list(
        estimate  = as_comp_estimate(rho_entries),
        bootstrap = NULL,
        ci        = NULL,
        grid      = list(primary = tk.plot, secondary = tk.plot2)
      )

      if (components) {
        rho_pearson_all <- .compute_LPC(pre)
        rhoP_entries    <- vector("list", ldb)
        Cb_entries      <- vector("list", ldb)

        for (i in seq_len(ldb)) {
          n_delta <- if (nd <= 1L) 1L else i
          rhoP_entries[[i]] <- rho_pearson_all[[n_delta]]
          beta_i            <- as.numeric(diffbeta[[i]])
          LA_list           <- .compute_LA(pre, diffbeta = beta_i)
          Cb_entries[[i]]   <- LA_list[[n_delta]]
        }

        metrics$lpc <- list(
          estimate  = as_comp_estimate(rhoP_entries),
          bootstrap = NULL,
          ci        = NULL,
          grid      = list(primary = tk.plot, secondary = tk.plot2)
        )
        metrics$la <- list(
          estimate  = as_comp_estimate(Cb_entries),
          bootstrap = NULL,
          ci        = NULL,
          grid      = list(primary = tk.plot, secondary = tk.plot2)
        )
      }
    }

    if (degenerate_resp) {
      metrics$lcc <- mask_metric_bundle(metrics$lcc)
      if (components) {
        metrics$lpc <- mask_metric_bundle(metrics$lpc)
        metrics$la  <- mask_metric_bundle(metrics$la)
      }
    }

    summary.lcc <- lccSummary(
      model      = model,
      tk         = tk,
      tk.plot    = tk.plot,
      tk.plot2   = tk.plot2,
      metrics    = metrics,
      ldb        = ldb,
      ci         = FALSE,
      components = components,
      degenerate_resp = degenerate_resp
    )

    CI <- NULL

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
      comp_names   = comp_names,
      boot.scheme  = boot.scheme,
      ci.method    = ci.method,
      alpha        = alpha,
      components   = components,
      lme.control  = lme.control,
      method.init  = method.init,
      numCore      = numCore,
      keep_models  = keep_models,
      rng_seed     = boot_seed
    )

    metrics <- CI$metrics
    if (!is.null(CI$comp_names)) {
      comp_names <- CI$comp_names
    }
    
    if (degenerate_resp) {
      metrics$lcc <- mask_metric_bundle(metrics$lcc)
      if (components) {
        metrics$lpc <- mask_metric_bundle(metrics$lpc)
        metrics$la  <- mask_metric_bundle(metrics$la)
      }
    }

    summary.lcc <- lccSummary(
      model      = model,
      tk         = tk,
      tk.plot    = tk.plot,
      tk.plot2   = tk.plot2,
      metrics    = metrics,
      ldb        = ldb,
      ci         = TRUE,
      components = components,
      degenerate_resp = degenerate_resp
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
    "covar"       = covar,
    "metrics"     = metrics,
    "comp_names"  = comp_names
  )
  if (isTRUE(ci)) {
    internal_lcc$bootstrap_state <- CI$bootstrap_state
    internal_lcc$alpha <- alpha
    internal_lcc$nboot <- nboot
  }
  
  if (ldb == 1L) {
    internal_lcc$rho <- metrics$lcc$estimate[[1L]]
    if (components) {
      internal_lcc$rho.pearson <- metrics$lpc$estimate[[1L]]
      internal_lcc$Cb          <- metrics$la$estimate[[1L]]
    }
  } else if (ldb > 1L) {
    extract_matrix <- function(est_list) {
      out <- as.data.frame(do.call(cbind, est_list))
      colnames(out) <- names(est_list)
      out
    }
    internal_lcc$rho <- extract_matrix(metrics$lcc$estimate)
    if (components) {
      internal_lcc$rho.pearson <- extract_matrix(metrics$lpc$estimate)
      internal_lcc$Cb          <- extract_matrix(metrics$la$estimate)
    }
  }

  if (isTRUE(ci)) {
    if (ldb == 1L) {
      internal_lcc$ENV.LCC <- metrics$lcc$ci[[1L]]
      if (components) {
        internal_lcc$ENV.LPC <- metrics$lpc$ci[[1L]]
        internal_lcc$ENV.LA  <- metrics$la$ci[[1L]]
      }
    } else {
      internal_lcc$ENV.LCC <- metrics$lcc$ci
      if (components) {
        internal_lcc$ENV.LPC <- metrics$lpc$ci
        internal_lcc$ENV.LA  <- metrics$la$ci
      }
    }
  }
  
  invisible(internal_lcc)
}

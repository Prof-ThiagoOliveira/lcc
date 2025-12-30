#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: summary_outputs.R                                             #
# Contains: lccSummary and fittedBuilder functions                    #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 22/11/2025                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal Function to Summarize Fitted and Sampled Values for \code{lcc} Objects
##'
##' @description Internally called helper that converts the unified metric bundles into the
##'   structures consumed by downstream plotting and summary routines.
##'
##' @param model the fitted \code{nlme} model stored within the \code{lcc} object.
##' @param tk original time grid from the data.
##' @param tk.plot time grid used for fitted trajectories.
##' @param tk.plot2 time grid used for sampled trajectories.
##' @param metrics list of metric bundles (LCC, LPC, LA) produced by \code{ciBuilder} or
##'   the no-bootstrap path.
##' @param ldb integer giving the number of method comparisons.
##' @param ci logical flag indicating whether confidence intervals were requested.
##' @param components logical flag indicating whether LPC/LA components were requested.
##'
##' @importFrom stats predict
##'
##' @keywords internal
lccSummary <- function(model, tk, tk.plot, tk.plot2, metrics,
                       ldb, ci, components,
                       degenerate_resp = FALSE) {

  method_levels <- levels(model$data$method)

  if (ldb == 1L) {
    comp <- paste0(method_levels[2L], " vs. ", method_levels[1L])
  } else {
    comp <- vector("list", ldb)
    for (i in seq_len(ldb)) {
      comp[[i]] <- paste0(method_levels[i + 1L], " vs. ", method_levels[1L])
    }
  }

  if (isTRUE(degenerate_resp)) {
    GF <- NA_real_
    CCC_vals <- lapply(
      seq_len(ldb),
      function(i) data.frame(V1 = rep(NA_real_, length(tk.plot2)))
    )
  } else {
    GF <- CCC(stats::predict(model), model$data$resp)

    CCC_vals <- CCC_lin(
      dataset = model$data,
      resp    = "resp",
      subject = "subject",
      method  = "method",
      time    = "time"
    )
  }

  to_comp_list <- function(values) {
    extract_numeric <- function(x) {
      if (is.null(x)) {
        return(NA_real_)
      }
      if (is.list(x) && length(x) == 1L &&
          (is.data.frame(x[[1L]]) || is.matrix(x[[1L]]) || is.list(x[[1L]]) || is.numeric(x[[1L]]))) {
        return(extract_numeric(x[[1L]]))
      }
      if (is.data.frame(x)) {
        return(as.numeric(x[[1L]]))
      }
      if (is.matrix(x)) {
        return(as.numeric(x[, 1L]))
      }
      as.numeric(x)
    }

    if (ldb == 1L) {
      return(list(extract_numeric(values)))
    }

    out <- NULL
    if (is.list(values) && length(values) == ldb) {
      out <- lapply(values, extract_numeric)
    } else if (is.data.frame(values) || is.matrix(values)) {
      out <- lapply(seq_len(ldb), function(j) extract_numeric(values[, j]))
    } else {
      out <- lapply(seq_len(ldb), function(j) extract_numeric(values[[j]]))
    }

    out
  }

  extract_ci_cols <- function(ci_mat) {
    if (is.null(ci_mat) || !is.matrix(ci_mat) || nrow(ci_mat) < 2L) {
      return(NULL)
    }
    data.frame(
      Lower = as.numeric(ci_mat["lower", , drop = TRUE]),
      Upper = as.numeric(ci_mat["upper", , drop = TRUE]),
      check.names = FALSE
    )
  }

  build_metric_df <- function(est_vec, ci_mat, value_name) {
    df <- data.frame(
      Time  = tk.plot,
      value = as.numeric(est_vec),
      check.names = FALSE
    )
    names(df)[2L] <- value_name
    ci_cols <- extract_ci_cols(ci_mat)
    if (!is.null(ci_cols)) {
      df <- cbind(df, ci_cols)
    }
    df
  }

  build_metric_list <- function(est_list, ci_list, value_name) {
    out <- vector("list", length(est_list))
    for (i in seq_along(est_list)) {
      ci_mat <- if (!is.null(ci_list)) ci_list[[i]] else NULL
      out[[i]] <- build_metric_df(est_list[[i]], ci_mat, value_name)
    }
    out
  }

  lcc_est <- metrics$lcc$estimate
  lcc_ci  <- metrics$lcc$ci

  if (!components) {
    if (ldb == 1L) {
      fitted <- build_metric_df(lcc_est[[1L]], if (!is.null(lcc_ci)) lcc_ci[[1L]] else NULL, "LCC")
    } else {
      fitted <- build_metric_list(lcc_est, lcc_ci, "LCC")
    }

    sampled <- data.frame(
      Time = tk.plot2,
      CCC  = CCC_vals
    )

    plot.data <- list(
      fitted  = fitted,
      sampled = sampled,
      gof     = GF,
      comp    = comp
    )

    return(invisible(plot.data))
  }

  if (isTRUE(degenerate_resp)) {
    Pearson_vals <- lapply(
      seq_len(ldb),
      function(i) data.frame(V1 = rep(NA_real_, length(tk.plot2)))
    )
  } else {
    Pearson_vals <- Pearson(
      dataset = model$data,
      resp    = "resp",
      subject = "subject",
      method  = "method",
      time    = "time"
    )
  }

  lpc_est <- metrics$lpc$estimate
  lpc_ci  <- metrics$lpc$ci
  la_est  <- metrics$la$estimate
  la_ci   <- metrics$la$ci

  CCC_list     <- to_comp_list(CCC_vals)
  Pearson_list <- to_comp_list(Pearson_vals)
  LA_sample    <- mapply(
    function(ccc, pear) ccc / pear,
    CCC_list,
    Pearson_list,
    SIMPLIFY = FALSE
  )

  if (!ci) {
    if (ldb == 1L) {
      fitted <- data.frame(
        Time = tk.plot,
        LCC  = as.numeric(lcc_est[[1L]]),
        LPC  = as.numeric(lpc_est[[1L]]),
        LA   = as.numeric(la_est[[1L]]),
        check.names = FALSE
      )
    } else {
      fitted <- vector("list", ldb)
      for (i in seq_len(ldb)) {
        fitted[[i]] <- data.frame(
          Time = tk.plot,
          LCC  = as.numeric(lcc_est[[i]]),
          LPC  = as.numeric(lpc_est[[i]]),
          LA   = as.numeric(la_est[[i]]),
          check.names = FALSE
        )
      }
    }

    if (ldb == 1L) {
      sampled <- data.frame(
        Time    = tk.plot2,
        CCC     = as.numeric(CCC_list[[1L]]),
        Pearson = as.numeric(Pearson_list[[1L]]),
        Cb      = as.numeric(LA_sample[[1L]]),
        check.names = FALSE
      )
    } else {
      sampled <- vector("list", ldb)
      for (i in seq_len(ldb)) {
        sampled[[i]] <- data.frame(
          Time    = tk.plot2,
          CCC     = as.numeric(CCC_list[[i]]),
          Pearson = as.numeric(Pearson_list[[i]]),
          Cb      = as.numeric(LA_sample[[i]]),
          check.names = FALSE
        )
      }
    }

    plot.data <- list(
      fitted  = fitted,
      sampled = sampled,
      gof     = GF,
      comp    = comp
    )

    return(invisible(plot.data))
  }

  if (ldb == 1L) {
    fitted <- list(
      LCC = build_metric_df(lcc_est[[1L]], if (!is.null(lcc_ci)) lcc_ci[[1L]] else NULL, "LCC"),
      LPC = build_metric_df(lpc_est[[1L]], if (!is.null(lpc_ci)) lpc_ci[[1L]] else NULL, "LPC"),
      LA  = build_metric_df(la_est[[1L]],  if (!is.null(la_ci))  la_ci[[1L]]  else NULL, "LA")
    )
    sampled <- data.frame(
      Time    = tk.plot2,
      CCC     = as.numeric(CCC_list[[1L]]),
      Pearson = as.numeric(Pearson_list[[1L]]),
      Cb      = as.numeric(LA_sample[[1L]]),
      check.names = FALSE
    )
  } else {
    fitted <- list(
      LCC = build_metric_list(lcc_est, lcc_ci, "LCC"),
      LPC = build_metric_list(lpc_est, lpc_ci, "LPC"),
      LA  = build_metric_list(la_est,  la_ci,  "LA")
    )
    sampled <- vector("list", ldb)
    for (i in seq_len(ldb)) {
      sampled[[i]] <- data.frame(
        Time    = tk.plot2,
        CCC     = as.numeric(CCC_list[[i]]),
        Pearson = as.numeric(Pearson_list[[i]]),
        Cb      = as.numeric(LA_sample[[i]]),
        check.names = FALSE
      )
    }
  }

  plot.data <- list(
    fitted  = fitted,
    sampled = sampled,
    gof     = GF,
    comp    = comp
  )

  invisible(plot.data)
}

#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# Section: fittedBuilder                                              #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 22/11/2025                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal Function to Build Fitted Values for
##'   \code{lcc} Objects
##'
##' @description This is an internally called function used to build
##'   fitted values.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @keywords internal
fittedBuilder <- function(object, type) {
  # Map type -> index and column name
  .form <- switch(type,
                  "lcc" = 1L,
                  "lpc" = 2L,
                  "la"  = 3L)
  .name <- switch(type,
                  "lcc" = "LCC",
                  "lpc" = "LPC",
                  "la"  = "LA")
  
  comp   <- object$Summary.lcc$comp
  fitted <- object$Summary.lcc$fitted
  
  #------------------------------------------------------------------
  # Single comparison (ldb == 1): comp is a character vector
  #------------------------------------------------------------------
  if (inherits(comp, "character")) {
    # fitted is either:
    #  - data.frame (no components or ci==FALSE)
    #  - list(LCC = df, LPC = df, LA = df) when components=TRUE & ci=TRUE
    df <- if (inherits(fitted, "data.frame")) {
      fitted
    } else {
      fitted[[.form]]
    }
    
    ret <- data.frame(
      Methods = comp,
      Time    = df[, "Time"],
      LCC     = df[, .name]
    )
    
  } else {
    #----------------------------------------------------------------
    # Multiple comparisons (ldb > 1): comp is a list of labels
    #----------------------------------------------------------------
    # 'fitted' is either:
    #  - list of data.frames (no components or ci==FALSE)
    #  - list(LCC = list(df_i), LPC = list(df_i), LA = list(df_i))
    #    when components=TRUE & ci=TRUE
    if (is.null(fitted$LCC)) {
      # No named components: list of data.frames for this type
      df_list <- fitted
    } else {
      # Named components: pick the list corresponding to this type
      df_list <- fitted[[.form]]
    }
    
    ncomp <- length(comp)
    fit   <- vector("list", ncomp)
    rn    <- vector("list", ncomp)
    met   <- vector("list", ncomp)
    
    for (i in seq_len(ncomp)) {
      df_i <- df_list[[i]]
      
      fit[[i]] <- data.frame(LCC = df_i[, .name])
      rn[[i]]  <- data.frame(Time = df_i[, "Time"])
      met[[i]] <- data.frame(
        Methods = rep(comp[[i]], nrow(df_i))
      )
    }
    
    ret <- data.frame(
      do.call(rbind.data.frame, met),
      do.call(rbind.data.frame, rn),
      do.call(rbind.data.frame, fit)
    )
  }
  
  # Keep the final column name convention identical to original
  colnames(ret)[colnames(ret) == "LCC"] <- paste0("fitted.", .name)
  
  ret
}


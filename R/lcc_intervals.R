#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: lcc_intervals.R                                               #
# Contains: lcc_intervals  function                                   #
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

##' @title Internal Functions to Compute the Non-Parametric Confidence
##'   Intervals for LCC.
##'
##' @description This is an internally called functions used to compute
##'   the non-parametric confidence intervals for LCC.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @importFrom stats quantile sd qnorm
##'
##' @keywords internal
##' @keywords internal
lcc_intervals <- function(rho, tk.plot, tk.plot2, ldb, model, ci,
                          percentileMet, LCC_Boot, alpha) {
  ## Fisher z transform and its inverse
  ZFisher     <- function(x) 0.5 * log((1 + x) / (1 - x))
  ZFisher_inv <- function(x) (exp(2 * x) - 1) / (exp(2 * x) + 1)
  
  ## percentileMet may come as "TRUE"/"FALSE" or logical
  percentile <- isTRUE(percentileMet) || identical(percentileMet, "TRUE")
  
  if (ldb == 1L) {
    ## LCC_Boot is list over bootstrap samples; each element numeric over time
    ENV.LCC <- .build_ci_from_boot(
      boot_list     = LCC_Boot,
      alpha         = alpha,
      transform     = if (!percentile) ZFisher else NULL,
      inv_transform = if (!percentile) ZFisher_inv else NULL,
      percentile    = percentile
    )
  } else {
    ## LCC_Boot is list over bootstrap samples; each [[b]] is list over methods
    ENV.LCC <- vector("list", ldb)
    for (i in seq_len(ldb)) {
      boot_i <- lapply(LCC_Boot, function(x) if (!is.null(x)) x[[i]] else NULL)
      ENV.LCC[[i]] <- .build_ci_from_boot(
        boot_list     = boot_i,
        alpha         = alpha,
        transform     = if (!percentile) ZFisher else NULL,
        inv_transform = if (!percentile) ZFisher_inv else NULL,
        percentile    = percentile
      )
    }
  }
  
  CI.LCC <- list("rho" = rho, "ENV.LCC" = ENV.LCC)
  CI.LCC
}

##' @keywords Internal
.build_ci_from_boot <- function(boot_list, alpha,
                                transform = NULL, inv_transform = NULL,
                                percentile = FALSE) {
  ## boot_list: list over bootstrap replicates; each element is a
  ## numeric vector over time (or NULL for failed fits)
  
  ## Drop NULL or all-NA replicates
  boot_list <- boot_list[!vapply(boot_list, is.null, logical(1L))]
  if (!length(boot_list)) {
    return(matrix(NA_real_, nrow = 2L, ncol = 0L))
  }
  
  ## Coerce to matrix: rows = time points; cols = bootstrap replicates
  boot_mat <- do.call(cbind, boot_list)
  
  if (percentile) {
    lower <- apply(boot_mat, 1L, stats::quantile, probs = alpha / 2)
    upper <- apply(boot_mat, 1L, stats::quantile, probs = 1 - alpha / 2)
    ci    <- rbind(lower, upper)
    return(ci)
  }
  
  ## Normal-approximation on a transformed scale
  if (!is.null(transform)) {
    boot_mat <- transform(boot_mat)
  }
  
  se <- apply(boot_mat, 1L, sd)
  mu <- apply(boot_mat, 1L, mean)
  z  <- stats::qnorm(1 - alpha / 2)
  
  ci <- rbind(mu - z * se, mu + z * se)
  
  if (!is.null(inv_transform)) {
    ci <- inv_transform(ci)
  }
  
  ci
}

ZFisher <- function(x) 0.5 * log((1 + x) / (1 - x))
ZFisher_inv <- function(x) (exp(2 * x) - 1) / (exp(2 * x) + 1)

Arcsin <- function(x) asin(sqrt(x))
Arcsin_inv <- function(x) sign(x) * sin(x)^2

build_ci_metric <- function(boot_list, alpha,
                            transform, inv_transform,
                            percentileMet) {
  # percentileMet may be logical or "TRUE"/"FALSE"
  percentile <- isTRUE(percentileMet) || identical(percentileMet, "TRUE")
  
  .build_ci_from_boot(
    boot_list     = boot_list,
    alpha         = alpha,
    transform     = if (!percentile) transform else NULL,
    inv_transform = if (!percentile) inv_transform else NULL,
    percentile    = percentile
  )
}

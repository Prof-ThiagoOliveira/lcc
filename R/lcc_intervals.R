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
lcc_intervals <- function(rho, tk.plot, tk.plot2, ldb, model, ci,
                          percentileMet, LCC_Boot, alpha) {
  # Fisher z transform and its inverse
  ZFisher      <- function(x) 0.5 * log((1 + x) / (1 - x))
  ZFisher_inv  <- function(x) (exp(2 * x) - 1) / (exp(2 * x) + 1)
  
  # percentileMet may come as "TRUE"/"FALSE" or logical
  percentile <- isTRUE(percentileMet) || identical(percentileMet, "TRUE")
  
  if (ldb == 1L) {
    # LCC_Boot is a list over bootstrap samples, each element is a numeric vector over time
    ENV.LCC <- .build_ci_from_boot(
      boot_list     = LCC_Boot,
      alpha         = alpha,
      transform     = if (!percentile) ZFisher else NULL,
      inv_transform = if (!percentile) ZFisher_inv else NULL,
      percentile    = percentile
    )
    
  } else {
    # LCC_Boot is a list over bootstrap samples; each [[b]] is a list over methods (ldb)
    ENV.LCC <- vector("list", ldb)
    
    for (i in seq_len(ldb)) {
      # Collect i-th method from each bootstrap replicate, keeping NULLs if whole replicate is NULL
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
  
  # Keep return structure unchanged
  CI.LCC <- list("rho" = rho, "ENV.LCC" = ENV.LCC)
  CI.LCC
}


#' @keywords Internal
.build_ci_from_boot <- function(boot_list, alpha,
                                transform = NULL, inv_transform = NULL,
                                percentile = FALSE) {
  # boot_list: list of numeric vectors, one per bootstrap replicate
  # (elements may be NULL; they are dropped)
  
  if (length(boot_list) == 0L) {
    return(NULL)
  }
  
  not_null <- !vapply(boot_list, is.null, logical(1L))
  
  # if all bootstrap replicates are NULL, we can't build CI
  if (!any(not_null)) {
    return(NULL)
  }
  
  boot_list <- boot_list[not_null]
  
  # rows: time, cols: bootstrap replicates
  boot_mat <- do.call(cbind, boot_list)
  
  if (percentile) {
    ci <- apply(
      boot_mat, 1L, quantile,
      probs = c(alpha / 2, 1 - alpha / 2)
    )
    return(ci)
  }
  
  # transform-based normal approximation
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
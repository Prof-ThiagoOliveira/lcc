#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: validators.R                                                  #
# Purpose: Shared validation and numerical safety helpers             #
#                                                                     #
#######################################################################

#' Internal validation helpers
#'
#' Utility functions for consistent input checking and numerical safety used
#' across the package.
#'
#' @keywords internal
#' @name lcc_validators
#' @importFrom stats var
NULL

# --------------------------------------------------------------------
# Length / type / missingness checks
# --------------------------------------------------------------------

validate_equal_length <- function(x, y, name_x = "x", name_y = "y") {
  len_x <- length(x)
  len_y <- length(y)
  if (!identical(len_x, len_y)) {
    abort_input(
      "Arguments {.arg {name_x}} and {.arg {name_y}} must have the same length (got {.val {len_x}} vs {.val {len_y}})."
    )
  }
  invisible(TRUE)
}

validate_numeric_no_na <- function(x, name_x = "x") {
  if (!is.numeric(x)) {
    abort_input("Argument {.arg {name_x}} must be numeric.")
  }
  if (anyNA(x)) {
    abort_input("Argument {.arg {name_x}} contains missing values; remove or impute them before calling this function.")
  }
  if (!all(is.finite(x))) {
    abort_input("Argument {.arg {name_x}} contains non-finite values.")
  }
  invisible(TRUE)
}

validate_non_degenerate_var <- function(x, name_x = "x") {
  vx <- stats::var(x)
  if (!is.finite(vx) || vx <= 0) {
    warn_general(
      sprintf("Variance of '%s' is non-positive; returning NA_real_.", name_x)
    )
    return(FALSE)
  }
  TRUE
}

# --------------------------------------------------------------------
# Fisher transform safety wrappers
# --------------------------------------------------------------------

safe_clamp_r <- function(x, eps = 1e-8) {
  dims <- dim(x)
  clamped <- base::pmax(-1 + eps, base::pmin(1 - eps, x))
  if (!is.null(dims)) {
    dim(clamped) <- dims
  }
  clamped
}

safe_fisher <- function(r, eps = 1e-8) {
  dims <- dim(r)
  transformed <- ZFisher(safe_clamp_r(r, eps = eps))
  if (!is.null(dims)) {
    dim(transformed) <- dims
  }
  transformed
}

safe_fisher_inv <- function(z, eps = 1e-8) {
  dims <- dim(z)
  inv <- safe_clamp_r(ZFisher_inv(z), eps = eps)
  if (!is.null(dims)) {
    dim(inv) <- dims
  }
  inv
}
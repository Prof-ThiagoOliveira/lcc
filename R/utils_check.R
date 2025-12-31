#' Internal helpers for argument checking and diagnostics
#'
#' These helpers centralise all user-facing errors, warnings and
#' informational messages using the \pkg{cli} and \pkg{rlang} packages.
#'
#' @name utils_check
#' @keywords internal
#' @importFrom cli cli_abort cli_warn cli_inform format_inline
#' @importFrom rlang is_bool is_scalar_integerish is_string arg_match0 caller_arg caller_env abort
NULL

# -------------------------------------------------------------------
# Core abort / warn / inform wrappers
# -------------------------------------------------------------------

# Generic package error
abort_lcc <- function(message, ..., .subclass = NULL, .call = caller_env()) {
  msg <- cli::format_inline(message, ..., .envir = parent.frame())
  rlang::abort(msg, class = c(.subclass, "lcc_error"), call = .call)
}

# Errors specifically due to invalid user input
abort_input <- function(message, ..., .subclass = NULL, .call = caller_env()) {
  msg <- cli::format_inline(message, ..., .envir = parent.frame())
  rlang::abort(msg, class = c(.subclass, "lcc_error_input"), call = .call)
}

# Errors that should indicate internal bugs / unreachable states
abort_internal <- function(message, ..., .call = caller_env()) {
  header <- cli::format_inline("Internal error in {.pkg lcc}.")
  info   <- cli::format_inline(message, ..., .envir = parent.frame())
  msg    <- paste0(header, " ", info)
  rlang::abort(msg, class = "lcc_error_internal", call = .call)
}

# Generic warning
warn_general <- function(message, ..., .subclass = NULL) {
  if (!isTRUE(getOption("lcc.show.warnings", TRUE))) {
    return(invisible(NULL))
  }
  old <- options(cli.width = 1000L, width = 1000L)
  on.exit(options(old), add = TRUE)
  msg <- cli::format_inline(message, ..., .envir = parent.frame())
  cli::cli_warn(
    msg,
    class = c(.subclass, "lcc_warning")
  )
}

# Generic informational message (progress, convergence summaries, etc.)
inform_general <- function(message, ...) {
  old <- options(cli.width = 1000L, width = 1000L)
  on.exit(options(old), add = TRUE)
  cli::cli_inform(cli::format_inline(message, ..., .envir = parent.frame()))
}

# -------------------------------------------------------------------
# Argument / value check helpers
# -------------------------------------------------------------------

# Logical flag
check_flag <- function(x, arg = caller_arg(x)) {
  if (!rlang::is_bool(x)) {
    abort_input("Argument {.arg {arg}} must be a single TRUE or FALSE.")
  }
  invisible(x)
}

# Scalar integer (optionally bounded)
check_scalar_integer <- function(x,
                                 arg   = caller_arg(x),
                                 lower = -Inf,
                                 upper = Inf,
                                 allow_na = FALSE) {
  if (!rlang::is_scalar_integerish(x) || (!allow_na && is.na(x))) {
    abort_input("Argument {.arg {arg}} must be a single integer.")
  }
  if (!is.na(lower) && x < lower) {
    abort_input(
      "Argument {.arg {arg}} must be >= {.val {lower}} (got {.val {x}})."
    )
  }
  if (!is.na(upper) && x > upper) {
    abort_input(
      "Argument {.arg {arg}} must be <= {.val {upper}} (got {.val {x}})."
    )
  }
  invisible(x)
}

# Scalar numeric (optionally bounded)
check_scalar_numeric <- function(x,
                                 arg   = caller_arg(x),
                                 lower = -Inf,
                                 upper = Inf,
                                 allow_na = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || (!allow_na && is.na(x))) {
    abort_input("Argument {.arg {arg}} must be a single numeric value.")
  }
  if (!is.na(lower) && x < lower) {
    abort_input(
      "Argument {.arg {arg}} must be >= {.val {lower}} (got {.val {x}})."
    )
  }
  if (!is.na(upper) && x > upper) {
    abort_input(
      "Argument {.arg {arg}} must be <= {.val {upper}} (got {.val {x}})."
    )
  }
  invisible(x)
}

# Scalar character
check_character_scalar <- function(x, arg = caller_arg(x)) {
  if (!rlang::is_string(x)) {
    abort_input("Argument {.arg {arg}} must be a single character string.")
  }
  invisible(x)
}

# Choice among allowed values (character)
check_choice <- function(x, choices, arg = caller_arg(x)) {
  x <- rlang::arg_match0(x, choices)
  invisible(x)
}

# -------------------------------------------------------------------
# Data / column helpers
# -------------------------------------------------------------------

# Ensure columns exist in a data frame
check_has_columns <- function(data, cols, data_name = caller_arg(data)) {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols)) {
    abort_input(
      "Object {.arg {data_name}} is missing required column{?s} {.field {missing_cols}}."
    )
  }
  invisible(data)
}

# Check that a column is factor
check_is_factor_col <- function(data, col, label = col) {
  if (!is.factor(data[[col]])) {
    abort_input("Please, '{label}' variable should be factor")
  }
  invisible(data)
}

# Check that a column is numeric
check_is_numeric_col <- function(data, col, label = col) {
  if (!is.numeric(data[[col]])) {
    abort_input("Please, '{label}' variable should be numeric")
  }
  invisible(data)
}

# -------------------------------------------------------------------
# High-level helpers tailored to existing code patterns
# -------------------------------------------------------------------

# Validate polynomial degrees qf and qr used in LCC models
check_polynomial_degrees <- function(qf, qr) {
  check_scalar_integer(qf, arg = "qf", lower = 1L)
  check_scalar_integer(qr, arg = "qr", lower = 0L)
  if (qr > qf) {
    abort_input("'qr' should be less or equal 'qf'")
  }
  invisible(list(qf = qf, qr = qr))
}

# Validate method for REML/ML
check_REML_flag <- function(REML) {
  check_flag(REML, arg = "REML")
  invisible(REML)
}

# Validate number of cores
check_num_core <- function(numCore) {
  check_scalar_integer(numCore, arg = "numCore", lower = 1L)
  invisible(numCore)
}

# Validate gold-standard method (gs)
check_gs <- function(gs, dataset, method) {
  if (is.null(gs)) {
    return(invisible(gs))
  }
  check_character_scalar(gs, arg = "gs")
  levels_method <- levels(dataset[[method]])
  if (!gs %in% levels_method) {
    abort_input("There is no method level called: '{gs}'")
  }
  invisible(gs)
}

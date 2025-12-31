build_fixture_dataset <- function() {
  df <- expand.grid(
    Fruit = factor(sprintf("F%d", 1:4)),
    Method = factor(c("A", "B"), levels = c("A", "B")),
    Time = 0:2
  )
  df <- df[order(df$Fruit, df$Method, df$Time), ]

  df$H_mean <- with(
    df,
    as.numeric(Fruit) * 0.5 +
      ifelse(Method == "B", 0.4, 0) +
      0.2 * Time +
      ifelse(Method == "B", 0.1, 0) * Time
  )

  df
}

build_test_lcc <- local({
  cache <- new.env(parent = emptyenv())

  function(components = TRUE) {
    key <- if (isTRUE(components)) "with_comp" else "no_comp"
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }

    data <- build_fixture_dataset()

    fit <- lcc(
      data       = data,
      subject    = "Fruit",
      resp       = "H_mean",
      method     = "Method",
      time       = "Time",
      qf         = 1,
      qr         = 0,
      components = components,
      ci         = FALSE,
      show.warnings = FALSE,
      keep.boot.models = FALSE,
      numCore    = 1
    )

    assign(key, fit, envir = cache)
    fit
  }
})

use_namespace_stub <- function(name, replacement, ns = "lcc") {
  ns_env <- asNamespace(ns)
  was_locked <- bindingIsLocked(name, ns_env)
  if (was_locked) {
    unlockBinding(name, ns_env)
  }
  original <- get(name, envir = ns_env)
  assign(name, replacement, envir = ns_env)
  list(name = name, original = original, ns_env = ns_env, was_locked = was_locked)
}

restore_namespace_stub <- function(stub) {
  assign(stub$name, stub$original, envir = stub$ns_env)
  if (stub$was_locked) {
    lockBinding(stub$name, stub$ns_env)
  }
}

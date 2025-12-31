test_that("lccWrapper selects comparison based on available metrics", {
  assign("summary.mock_lme", function(object, ...) list(modelStruct = list(varStruct = NULL)), envir = .GlobalEnv)
  on.exit(rm("summary.mock_lme", envir = .GlobalEnv), add = TRUE)

  stub_cov <- use_namespace_stub(
    "extract_random_effects_cov",
    function(model) list(G = matrix(1), n_re = 2L)
  )
  on.exit(restore_namespace_stub(stub_cov), add = TRUE)

  pre_calls <- list()
  stub_pre <- use_namespace_stub(
    ".precompute_longitudinal",
    function(model, tk, q_f, q_r, ...) {
      pre_calls <<- append(pre_calls, list(list(tk = tk, q_f = q_f, q_r = q_r)))
      list(tag = paste0("pre-", model$label))
    }
  )
  on.exit(restore_namespace_stub(stub_pre), add = TRUE)

  rho_queue <- list(
    list(c(0.4, 0.5), c(0.6, 0.7)),
    list(c(0.2, 0.25)),
    list(c(0.8, 0.85), c(NA, 0.9))
  )
  stub_lcc <- use_namespace_stub(
    ".compute_LCC",
    function(pre, diffbeta) {
      res <- rho_queue[[1]]
      rho_queue <<- rho_queue[-1]
      res
    }
  )
  on.exit(restore_namespace_stub(stub_lcc), add = TRUE)

  model <- structure(list(label = "fit"), class = "mock_lme")

  val_first <- lcc:::lccWrapper(model, q_f = 1, tk = 0:2, diffbeta = list(c(0.1, 0.2)), n.delta = 2)
  expect_equal(val_first, c(0.6, 0.7))
  expect_equal(pre_calls[[1]]$q_r, 1)

  val_second <- lcc:::lccWrapper(model, q_f = 1, tk = 0:1, diffbeta = list(c(0.1)), n.delta = 1)
  expect_equal(val_second, c(0.2, 0.25))

  val_third <- lcc:::lccWrapper(model, q_f = 1, tk = 0:1, diffbeta = list(c(0.1)), n.delta = 2)
  expect_equal(val_third, c(0.8, 0.85))
})

test_that("lpcWrapper and laWrapper return requested levels", {
  assign("summary.mock_lme", function(object, ...) list(modelStruct = list(varStruct = NULL)), envir = .GlobalEnv)
  on.exit(rm("summary.mock_lme", envir = .GlobalEnv), add = TRUE)

  stub_cov <- use_namespace_stub(
    "extract_random_effects_cov",
    function(model) list(G = matrix(1), n_re = 3L)
  )
  on.exit(restore_namespace_stub(stub_cov), add = TRUE)

  stub_pre <- use_namespace_stub(
    ".precompute_longitudinal",
    function(model, tk, q_f, q_r, ...) list(tag = paste0("pre-", model$label))
  )
  on.exit(restore_namespace_stub(stub_pre), add = TRUE)

  stub_lpc <- use_namespace_stub(
    ".compute_LPC",
    function(pre) list(c(0.5, 0.55), c(0.6, 0.65))
  )
  on.exit(restore_namespace_stub(stub_lpc), add = TRUE)

  stub_la <- use_namespace_stub(
    ".compute_LA",
    function(pre, diffbeta) list(c(0.7, 0.72), c(0.75, 0.78))
  )
  on.exit(restore_namespace_stub(stub_la), add = TRUE)

  model <- structure(list(label = "fit"), class = "mock_lme")

  lpc_val <- lcc:::lpcWrapper(model, q_f = 1, tk = 0:2, n.delta = 2)
  expect_equal(lpc_val, c(0.6, 0.65))

  la_val <- lcc:::laWrapper(model, q_f = 1, tk = 0:2, diffbeta = list(c(0.1)), n.delta = 1)
  expect_equal(la_val, c(0.7, 0.72))
})

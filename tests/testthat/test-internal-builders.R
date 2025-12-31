test_that("laBuilder forwards arguments to precompute helpers", {
  inputs <- list()

  stub_pre_lbuilder <- use_namespace_stub(
    ".precompute_longitudinal",
    function(model, tk, q_f, q_r, ...) {
      inputs <<- append(inputs, list(list(model = model, tk = tk, q_f = q_f, q_r = q_r)))
      list(pre_id = "mock_pre")
    }
  )
  on.exit(restore_namespace_stub(stub_pre_lbuilder), add = TRUE)

  stub_la_lbuilder <- use_namespace_stub(
    ".compute_LA",
    function(pre, diffbeta) {
      list(pre = pre$pre_id, diff = diffbeta)
    }
  )
  on.exit(restore_namespace_stub(stub_la_lbuilder), add = TRUE)

  out <- lcc:::laBuilder(
    G = NULL,
    diffbeta = c(0.1, 0.2),
    tk = 0:1,
    q_r = 0,
    q_f = 1,
    g = NULL,
    sig2_epsilon = 1,
    delta = NULL,
    deltal = NULL,
    model = "model-tag"
  )

  expect_equal(inputs[[1]]$model, "model-tag")
  expect_equal(inputs[[1]]$tk, 0:1)
  expect_equal(out$pre, "mock_pre")
  expect_equal(out$diff, c(0.1, 0.2))
})

test_that("lccBuilder and lpcBuilder reuse precomputation", {
  lcc_calls <- list()

  stub_pre_lcc <- use_namespace_stub(
    ".precompute_longitudinal",
    function(model, tk, q_f, q_r, ...) {
      lcc_calls <<- append(lcc_calls, list(list(model = model, tk = tk, q_f = q_f, q_r = q_r)))
      list(pre_id = "pre")
    }
  )
  on.exit(restore_namespace_stub(stub_pre_lcc), add = TRUE)

  stub_compute_lcc <- use_namespace_stub(
    ".compute_LCC",
    function(pre, diffbeta) {
      list(seq_along(diffbeta))
    }
  )
  on.exit(restore_namespace_stub(stub_compute_lcc), add = TRUE)

  stub_compute_lpc <- use_namespace_stub(
    ".compute_LPC",
    function(pre) {
      list(0.5)
    }
  )
  on.exit(restore_namespace_stub(stub_compute_lpc), add = TRUE)

  res_lcc <- lcc:::lccBuilder(
    G = NULL,
    diffbeta = c(0.2, 0.1),
    tk = 0:1,
    q_r = 0,
    q_f = 1,
    g = NULL,
    sig2_epsilon = 1,
    delta = NULL,
    deltal = NULL,
    model = "m"
  )
  res_lpc <- lcc:::lpcBuilder(
    G = NULL,
    tk = 0:1,
    q_r = 0,
    q_f = 1,
    g = NULL,
    sig2_epsilon = 1,
    delta = NULL,
    deltal = NULL,
    model = "m"
  )

  expect_equal(lcc_calls[[1]]$model, "m")
  expect_identical(res_lcc, list(seq_len(2)))
  expect_identical(res_lpc, list(0.5))
})

test_that("bootstrap helpers skip NULL models and honour selectors", {
  loop_calls <- list()

  assign("summary.mock_lme", function(object, ...) object$summary_obj, envir = .GlobalEnv)
  on.exit(rm("summary.mock_lme", envir = .GlobalEnv), add = TRUE)
  stub_cov_boot <- use_namespace_stub(
    "extract_random_effects_cov",
    function(model) list(G = matrix(1), n_re = 1L)
  )
  on.exit(restore_namespace_stub(stub_cov_boot), add = TRUE)

  stub_pre_boot <- use_namespace_stub(
    ".precompute_longitudinal",
    function(model, tk, q_f, q_r, summary_obj, G_info) {
      list(pre_id = paste0("pre-", model$id))
    }
  )
  on.exit(restore_namespace_stub(stub_pre_boot), add = TRUE)

  stub_boot_lcc <- use_namespace_stub(
    ".compute_LCC",
    function(pre, diffbeta) sum(diffbeta)
  )
  on.exit(restore_namespace_stub(stub_boot_lcc), add = TRUE)

  stub_boot_lpc <- use_namespace_stub(
    ".compute_LPC",
    function(pre) 0.2
  )
  on.exit(restore_namespace_stub(stub_boot_lpc), add = TRUE)

  stub_boot_la <- use_namespace_stub(
    ".compute_LA",
    function(pre, diffbeta) sum(diffbeta) / 2
  )
  on.exit(restore_namespace_stub(stub_boot_la), add = TRUE)

  stub_boot_loop <- use_namespace_stub(
    ".bootstrap_metric_loop",
    function(pre, ldb, compute_metric_from_pre, diffbeta = NULL,
             use_delta_by_level, needs_diff, selector, label) {
      loop_calls <<- append(loop_calls, list(list(
        label = label,
        pre   = pre$pre_id,
        diff  = diffbeta,
        use_delta = use_delta_by_level,
        needs_diff = needs_diff
      )))
      compute_metric_from_pre(pre, if (needs_diff) diffbeta[[1L]] else NULL)
    }
  )
  on.exit(restore_namespace_stub(stub_boot_loop), add = TRUE)

  make_model <- function(id, var_len) {
    structure(
      list(
        id = id,
        summary_obj = list(modelStruct = list(varStruct = seq_len(var_len)))
      ),
      class = "mock_lme"
    )
  }

  model_boot <- list(
    make_model("m1", 1L),
    NULL,
    make_model("m3", 2L)
  )
  diff_boot <- list(list(c(0.1, 0.2)), list(c(0.1, 0.2)), list(c(0.3, 0.4)))

  lcc_out <- lcc:::lccBootstrap(model_boot, diff_boot, ldb = 1, nboot = 3, tk = 0:1, q_f = 1)
  lpc_out <- lcc:::lpcBootstrap(model_boot, ldb = 1, nboot = 3, tk = 0:1, q_f = 1)
  la_out  <- lcc:::laBootstrap(model_boot, diff_boot, ldb = 1, nboot = 3, tk = 0:1, q_f = 1)

  expect_equal(length(lcc_out), 3)
  expect_null(lcc_out[[2]])
  expect_equal(loop_calls[[1]]$label, "LCC")
  expect_false(loop_calls[[1]]$use_delta)
  expect_true(loop_calls[[1]]$needs_diff)
  lpc_m3 <- Filter(function(x) identical(x$label, "LPC") && identical(x$pre, "pre-m3"), loop_calls)
  expect_equal(length(lpc_m3), 1)
  expect_true(lpc_m3[[1]]$use_delta)
  expect_equal(lpc_out[[1]], 0.2)
  expect_equal(la_out[[1]], sum(diff_boot[[1]][[1]]) / 2)
})

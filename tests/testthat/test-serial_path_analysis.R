test_that("serial_path_analysis runs and selects expected model (m=3)", {
  skip_if_not_installed("data.table")

  set.seed(123)
  n <- 300

  X <- rnorm(n)
  M1 <- 0.6 * X + rnorm(n, sd = 0.8)
  M2 <- 0.4 * X + 0.5 * M1 + rnorm(n, sd = 0.8)
  M3 <- 0.3 * X + 0.4 * M1 + 0.6 * M2 + rnorm(n, sd = 0.8)
  # cp = 0 (Semifull is the data-generating model)
  Y <- 0.7 * M1 + 0.2 * M2 + 0.5 * M3 + rnorm(n, sd = 0.8)

  df <- data.frame(X = X, M1 = M1, M2 = M2, M3 = M3, Y = Y)

  tmp <- tempfile("serial_path_")
  out_dir <- dirname(tmp)
  out_prefix <- basename(tmp)

  serial_path_analysis(
    data = df,
    columns = list(
      x = "X",
      mediators = list(m1 = "M1", m2 = "M2", m3 = "M3"),
      y = "Y"
    ),
    output_dir = out_dir,
    output_prefix = out_prefix,
    nrep = 20,
    num_threads = 1,
    seed = 1,
    alpha = 0.05,
    match = "exact",
    chunk_size = 128,
    grain_size = 1,
    overwrite = TRUE,
    output = c("fit", "params", "effects")
  )

  fit <- data.table::fread(file.path(out_dir, paste0(out_prefix, "_fit.csv")))
  params <- data.table::fread(file.path(out_dir, paste0(out_prefix, "_params.csv")))
  effects <- data.table::fread(file.path(out_dir, paste0(out_prefix, "_effects.csv")))

  expect_equal(nrow(fit), 1)
  expect_equal(fit$m[[1]], 3)
  expect_equal(fit$semifull_df[[1]], 1)
  expect_equal(fit$full_df[[1]], 6)

  # With this fixed seed/data, Semifull should be selected.
  expect_equal(fit$selected_model[[1]], "Semifull")

  cp <- params[label == "cp"]
  expect_equal(nrow(cp), 1)
  expect_lt(abs(cp$est[[1]]), 0.2)

  total_ind <- effects[effect == "total_indirect"]
  expect_equal(nrow(total_ind), 1)
  expect_gt(total_ind$est[[1]], 0)
})

test_that("df=0 fit indices are fixed by contract (Partial)", {
  skip_if_not_installed("data.table")

  set.seed(1)
  n <- 200
  X <- rnorm(n)
  M1 <- 0.5 * X + rnorm(n)
  M2 <- 0.3 * X + 0.6 * M1 + rnorm(n)
  Y <- 0.4 * X + 0.7 * M1 + 0.2 * M2 + rnorm(n)
  df <- data.frame(X = X, M1 = M1, M2 = M2, Y = Y)

  tmp <- tempfile("serial_path_df0_")
  out_dir <- dirname(tmp)
  out_prefix <- basename(tmp)

  serial_path_analysis(
    data = df,
    columns = list(x = "X", mediators = list(m1 = "M1", m2 = "M2"), y = "Y"),
    output_dir = out_dir,
    output_prefix = out_prefix,
    nrep = 0,
    num_threads = 1,
    seed = 1,
    match = "exact",
    overwrite = TRUE,
    output = "fit"
  )

  fit <- data.table::fread(file.path(out_dir, paste0(out_prefix, "_fit.csv")))
  expect_equal(nrow(fit), 1)
  expect_equal(fit$partial_df[[1]], 0)
  expect_equal(fit$partial_chisq[[1]], 0)
  expect_equal(fit$partial_p[[1]], 1)
  expect_equal(fit$partial_cfi[[1]], 1)
  expect_equal(fit$partial_tli[[1]], 1)
  expect_equal(fit$partial_rmsea[[1]], 0)
})

test_that("serial_path_analysis runs for m=1 and writes expected df/effects", {
  skip_if_not_installed("data.table")

  set.seed(11)
  n <- 250

  X <- rnorm(n)
  M1 <- 0.6 * X + rnorm(n, sd = 0.8)
  Y <- 0.7 * M1 + rnorm(n, sd = 0.8)

  df <- data.frame(X = X, M1 = M1, Y = Y)

  tmp <- tempfile("serial_path_m1_")
  out_dir <- dirname(tmp)
  out_prefix <- basename(tmp)

  serial_path_analysis(
    data = df,
    columns = list(x = "X", mediators = list(m1 = "M1"), y = "Y"),
    output_dir = out_dir,
    output_prefix = out_prefix,
    nrep = 0,
    num_threads = 1,
    seed = 1,
    match = "exact",
    overwrite = TRUE,
    output = c("fit", "effects")
  )

  fit <- data.table::fread(file.path(out_dir, paste0(out_prefix, "_fit.csv")))
  effects <- data.table::fread(file.path(out_dir, paste0(out_prefix, "_effects.csv")))

  expect_equal(nrow(fit), 1)
  expect_equal(fit$m[[1]], 1)
  expect_equal(fit$full_df[[1]], 1)

  expect_true(any(effects$effect == "x_to_m1_total"))
  expect_true(any(effects$effect == "adjacent_chain"))
})

test_that("m missingness semantics: infer from stage list when m omitted", {
  skip_if_not_installed("data.table")

  set.seed(2)
  n <- 250
  X <- rnorm(n)
  M1 <- 0.5 * X + rnorm(n)
  M2 <- 0.4 * X + 0.6 * M1 + rnorm(n)
  M3 <- 0.3 * X + 0.2 * M1 + 0.7 * M2 + rnorm(n)
  Y <- 0.2 * X + 0.3 * M1 + 0.4 * M2 + 0.5 * M3 + rnorm(n)
  df <- data.frame(X = X, M1 = M1, M2 = M2, M3 = M3, Y = Y)

  tmp <- tempfile("serial_path_m_missing_")
  out_dir <- dirname(tmp)
  out_prefix <- basename(tmp)

  # m is omitted; infer m=3 from list length.
  serial_path_analysis(
    data = df,
    columns = list(x = "X", mediators = list(m1 = "M1", m2 = "M2", m3 = "M3"), y = "Y"),
    output_dir = out_dir,
    output_prefix = out_prefix,
    nrep = 0,
    num_threads = 1,
    seed = 1,
    match = "exact",
    overwrite = TRUE,
    output = "fit"
  )

  fit <- data.table::fread(file.path(out_dir, paste0(out_prefix, "_fit.csv")))
  expect_equal(fit$m[[1]], 3)

  # If user explicitly supplies conflicting m, error.
  expect_error(
    serial_path_analysis(
      data = df,
      columns = list(x = "X", mediators = list(m1 = "M1", m2 = "M2", m3 = "M3"), y = "Y"),
      m = 2,
      output_dir = out_dir,
      output_prefix = paste0(out_prefix, "_conflict"),
      nrep = 0,
      num_threads = 1,
      seed = 1,
      match = "exact",
      overwrite = TRUE,
      output = "fit"
    ),
    "does not match"
  )
})

test_that("serial_path_analysis runs for m=4", {
  skip_if_not_installed("data.table")

  set.seed(12)
  n <- 300

  X <- rnorm(n)
  M1 <- 0.6 * X + rnorm(n, sd = 0.8)
  M2 <- 0.5 * X + 0.4 * M1 + rnorm(n, sd = 0.8)
  M3 <- 0.4 * X + 0.3 * M1 + 0.5 * M2 + rnorm(n, sd = 0.8)
  M4 <- 0.3 * X + 0.2 * M2 + 0.6 * M3 + rnorm(n, sd = 0.8)
  Y <- 0.1 * X + 0.3 * M1 + 0.2 * M2 + 0.4 * M3 + 0.5 * M4 + rnorm(n, sd = 0.8)

  df <- data.frame(X = X, M1 = M1, M2 = M2, M3 = M3, M4 = M4, Y = Y)

  tmp <- tempfile("serial_path_m4_")
  out_dir <- dirname(tmp)
  out_prefix <- basename(tmp)

  serial_path_analysis(
    data = df,
    columns = list(x = "X", mediators = list(m1 = "M1", m2 = "M2", m3 = "M3", m4 = "M4"), y = "Y"),
    output_dir = out_dir,
    output_prefix = out_prefix,
    nrep = 0,
    num_threads = 1,
    seed = 1,
    match = "exact",
    overwrite = TRUE,
    output = "fit"
  )

  fit <- data.table::fread(file.path(out_dir, paste0(out_prefix, "_fit.csv")))
  expect_equal(nrow(fit), 1)
  expect_equal(fit$m[[1]], 4)
  expect_equal(fit$full_df[[1]], 10)
})

test_that("default m=2 when mediators is a shared pool and m omitted", {
  skip_if_not_installed("data.table")

  set.seed(3)
  n <- 200
  X <- rnorm(n)
  M1 <- 0.4 * X + rnorm(n)
  M2 <- 0.3 * X + 0.5 * M1 + rnorm(n)
  M3 <- 0.2 * X + 0.2 * M1 + 0.4 * M2 + rnorm(n)
  Y <- 0.3 * X + 0.3 * M1 + 0.3 * M2 + 0.3 * M3 + rnorm(n)
  df <- data.frame(X = X, M1 = M1, M2 = M2, M3 = M3, Y = Y)

  tmp <- tempfile("serial_path_default_m2_")
  out_dir <- dirname(tmp)
  out_prefix <- basename(tmp)

  serial_path_analysis(
    data = df,
    columns = list(x = "X", mediators = c("M1", "M2", "M3"), y = "Y"),
    output_dir = out_dir,
    output_prefix = out_prefix,
    nrep = 0,
    num_threads = 1,
    seed = 1,
    match = "exact",
    overwrite = TRUE,
    output = "fit"
  )

  fit <- data.table::fread(file.path(out_dir, paste0(out_prefix, "_fit.csv")))
  expect_true(all(fit$m == 2))
  expect_equal(nrow(fit), 3 * 3) # 3 mediators in each of 2 stages
})

test_that("combination guard stops early when n_to_process exceeds max_combinations", {
  skip_if_not_installed("data.table")

  set.seed(4)
  n <- 50
  df <- as.data.frame(matrix(rnorm(n * 12), nrow = n))
  names(df) <- c(
    "X1", "X2",
    "M1", "M2", "M3", "M4",
    "Y1", "Y2",
    "Z1", "Z2", "Z3", "Z4"
  )

  tmp <- tempfile("serial_path_guard_")
  out_dir <- dirname(tmp)
  out_prefix <- basename(tmp)

  # nx=2; shared mediator pool size=4 expanded to m=2 => 4*4; ny=2 => total=64
  expect_error(
    serial_path_analysis(
      data = df,
      columns = list(x = c("X1", "X2"), mediators = c("M1", "M2", "M3", "M4"), y = c("Y1", "Y2")),
      output_dir = out_dir,
      output_prefix = out_prefix,
      nrep = 0,
      num_threads = 1,
      seed = 1,
      match = "exact",
      overwrite = TRUE,
      max_combinations = 10,
      output = "fit"
    ),
    "exceeds max_combinations"
  )
})

test_that("output selection controls file creation", {
  skip_if_not_installed("data.table")

  set.seed(5)
  n <- 200
  X <- rnorm(n)
  M1 <- 0.5 * X + rnorm(n)
  M2 <- 0.4 * X + 0.6 * M1 + rnorm(n)
  Y <- 0.2 * X + 0.3 * M1 + 0.4 * M2 + rnorm(n)
  df <- data.frame(X = X, M1 = M1, M2 = M2, Y = Y)

  tmp <- tempfile("serial_path_output_select_")
  out_dir <- dirname(tmp)
  out_prefix <- basename(tmp)

  serial_path_analysis(
    data = df,
    columns = list(x = "X", mediators = list(m1 = "M1", m2 = "M2"), y = "Y"),
    output_dir = out_dir,
    output_prefix = out_prefix,
    nrep = 0,
    num_threads = 1,
    seed = 1,
    match = "exact",
    overwrite = TRUE,
    output = "fit"
  )

  expect_true(file.exists(file.path(out_dir, paste0(out_prefix, "_fit.csv"))))
  expect_false(file.exists(file.path(out_dir, paste0(out_prefix, "_params.csv"))))
  expect_false(file.exists(file.path(out_dir, paste0(out_prefix, "_effects.csv"))))
})

test_that("overwrite=FALSE stops before running when files exist", {
  skip_if_not_installed("data.table")

  set.seed(6)
  n <- 50
  X <- rnorm(n)
  M1 <- 0.5 * X + rnorm(n)
  M2 <- 0.4 * X + 0.6 * M1 + rnorm(n)
  Y <- 0.2 * X + 0.3 * M1 + 0.4 * M2 + rnorm(n)
  df <- data.frame(X = X, M1 = M1, M2 = M2, Y = Y)

  tmp <- tempfile("serial_path_overwrite_false_")
  out_dir <- dirname(tmp)
  out_prefix <- basename(tmp)
  fit_path <- file.path(out_dir, paste0(out_prefix, "_fit.csv"))

  writeLines("sentinel\n", fit_path)
  info_before <- file.info(fit_path)

  expect_error(
    serial_path_analysis(
      data = df,
      columns = list(x = "X", mediators = list(m1 = "M1", m2 = "M2"), y = "Y"),
      output_dir = out_dir,
      output_prefix = out_prefix,
      nrep = 0,
      num_threads = 1,
      seed = 1,
      match = "exact",
      overwrite = FALSE,
      output = "fit"
    ),
    "already exist"
  )

  info_after <- file.info(fit_path)
  expect_equal(info_after$size, info_before$size)
})

test_that("combination id uses canonical format and escapes %,|,=", {
  skip_if_not_installed("data.table")

  set.seed(7)
  n <- 200
  df <- data.frame(
    `X|a` = rnorm(n),
    `M=1` = rnorm(n),
    `M%2` = rnorm(n),
    `Y%1` = rnorm(n),
    check.names = FALSE
  )

  tmp <- tempfile("serial_path_combo_id_")
  out_dir <- dirname(tmp)
  out_prefix <- basename(tmp)

  serial_path_analysis(
    data = df,
    columns = list(x = "X|a", mediators = list(m1 = "M=1", m2 = "M%2"), y = "Y%1"),
    output_dir = out_dir,
    output_prefix = out_prefix,
    nrep = 0,
    num_threads = 1,
    seed = 1,
    match = "exact",
    overwrite = TRUE,
    output = "fit"
  )

  fit <- data.table::fread(file.path(out_dir, paste0(out_prefix, "_fit.csv")))
  expect_equal(nrow(fit), 1)
  expect_equal(
    fit$Combination[[1]],
    "X=X%7Ca|M1=M%3D1|M2=M%252|Y=Y%251"
  )
})

test_that("duplicate columns in a combination yield a failed (NA) fit row", {
  skip_if_not_installed("data.table")

  set.seed(8)
  n <- 100
  X <- rnorm(n)
  M2 <- 0.5 * X + rnorm(n)
  Y <- 0.3 * X + 0.4 * M2 + rnorm(n)
  df <- data.frame(X = X, M2 = M2, Y = Y)

  tmp <- tempfile("serial_path_dups_")
  out_dir <- dirname(tmp)
  out_prefix <- basename(tmp)

  serial_path_analysis(
    data = df,
    columns = list(x = "X", mediators = list(m1 = "X", m2 = "M2"), y = "Y"),
    output_dir = out_dir,
    output_prefix = out_prefix,
    nrep = 0,
    num_threads = 1,
    seed = 1,
    match = "exact",
    overwrite = TRUE,
    output = "fit"
  )

  fit <- data.table::fread(file.path(out_dir, paste0(out_prefix, "_fit.csv")))
  expect_equal(nrow(fit), 1)
  expect_true(is.na(fit$selected_model[[1]]))
})

test_that("include_failure_reason annotates failed and successful fit rows", {
  skip_if_not_installed("data.table")

  set.seed(81)
  n <- 120
  X <- rnorm(n)
  M1 <- 0.5 * X + rnorm(n)
  M2 <- 0.4 * X + 0.6 * M1 + rnorm(n)
  Y <- 0.2 * X + 0.3 * M1 + 0.4 * M2 + rnorm(n)
  df <- data.frame(X = X, M1 = M1, M2 = M2, Y = Y)

  tmp <- tempfile("serial_path_failure_reason_")
  out_dir <- dirname(tmp)
  out_prefix <- basename(tmp)

  serial_path_analysis(
    data = df,
    columns = list(x = "X", mediators = list(m1 = c("X", "M1"), m2 = "M2"), y = "Y"),
    output_dir = out_dir,
    output_prefix = out_prefix,
    nrep = 0,
    num_threads = 1,
    seed = 1,
    match = "exact",
    overwrite = TRUE,
    include_failure_reason = TRUE,
    output = "fit"
  )

  fit <- data.table::fread(file.path(out_dir, paste0(out_prefix, "_fit.csv")))
  expect_true("failure_reason" %in% names(fit))
  expect_equal(nrow(fit), 2)

  failed <- fit[is.na(selected_model)]
  succeeded <- fit[!is.na(selected_model)]
  expect_equal(nrow(failed), 1)
  expect_equal(failed$failure_reason[[1]], "duplicate_columns")
  expect_true(all(is.na(succeeded$failure_reason)))
})

test_that("sharding partitions global indices deterministically", {
  skip_if_not_installed("data.table")

  set.seed(9)
  n <- 250
  df <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n),
    M1 = rnorm(n),
    M2 = rnorm(n),
    Y = rnorm(n)
  )

  columns <- list(x = c("X1", "X2"), mediators = list(m1 = c("M1", "M2"), m2 = c("M1", "M2")), y = "Y")

  tmp_full <- tempfile("serial_path_shard_full_")
  out_dir <- dirname(tmp_full)
  prefix_full <- basename(tmp_full)

  serial_path_analysis(
    data = df,
    columns = columns,
    output_dir = out_dir,
    output_prefix = prefix_full,
    nrep = 5,
    num_threads = 2,
    seed = 42,
    match = "exact",
    overwrite = TRUE,
    output = "fit"
  )
  fit_full <- data.table::fread(file.path(out_dir, paste0(prefix_full, "_fit.csv")))

  tmp0 <- tempfile("serial_path_shard0_")
  tmp1 <- tempfile("serial_path_shard1_")
  tmp2 <- tempfile("serial_path_shard2_")

  serial_path_analysis(
    data = df,
    columns = columns,
    output_dir = dirname(tmp0),
    output_prefix = basename(tmp0),
    nrep = 5,
    num_threads = 1,
    seed = 42,
    match = "exact",
    shard_id = 0,
    shard_count = 3,
    overwrite = TRUE,
    output = "fit"
  )
  serial_path_analysis(
    data = df,
    columns = columns,
    output_dir = dirname(tmp1),
    output_prefix = basename(tmp1),
    nrep = 5,
    num_threads = 1,
    seed = 42,
    match = "exact",
    shard_id = 1,
    shard_count = 3,
    overwrite = TRUE,
    output = "fit"
  )
  serial_path_analysis(
    data = df,
    columns = columns,
    output_dir = dirname(tmp2),
    output_prefix = basename(tmp2),
    nrep = 5,
    num_threads = 1,
    seed = 42,
    match = "exact",
    shard_id = 2,
    shard_count = 3,
    overwrite = TRUE,
    output = "fit"
  )

  fit0 <- data.table::fread(file.path(dirname(tmp0), paste0(basename(tmp0), "_fit.csv")))
  fit1 <- data.table::fread(file.path(dirname(tmp1), paste0(basename(tmp1), "_fit.csv")))
  fit2 <- data.table::fread(file.path(dirname(tmp2), paste0(basename(tmp2), "_fit.csv")))

  fit_shards <- data.table::rbindlist(list(fit0, fit1, fit2), use.names = TRUE, fill = TRUE)

  data.table::setorder(fit_full, Combination)
  data.table::setorder(fit_shards, Combination)
  expect_equal(fit_shards$Combination, fit_full$Combination)
  expect_equal(fit_shards$selected_model, fit_full$selected_model)
  expect_equal(fit_shards$full_chisq, fit_full$full_chisq)
})

test_that("results are deterministic across num_threads/grain_size/chunk_size", {
  skip_if_not_installed("data.table")

  set.seed(10)
  n <- 300
  df <- data.frame(
    X = rnorm(n),
    M1 = rnorm(n),
    M2 = rnorm(n),
    Y = rnorm(n)
  )

  columns <- list(x = "X", mediators = list(m1 = "M1", m2 = "M2"), y = "Y")

  tmp_a <- tempfile("serial_path_det_a_")
  out_dir <- dirname(tmp_a)
  prefix_a <- basename(tmp_a)

  serial_path_analysis(
    data = df,
    columns = columns,
    output_dir = out_dir,
    output_prefix = prefix_a,
    nrep = 10,
    num_threads = 1,
    seed = 99,
    match = "exact",
    chunk_size = 2,
    grain_size = 1,
    overwrite = TRUE,
    output = c("fit", "params", "effects")
  )

  tmp_b <- tempfile("serial_path_det_b_")
  prefix_b <- basename(tmp_b)
  serial_path_analysis(
    data = df,
    columns = columns,
    output_dir = out_dir,
    output_prefix = prefix_b,
    nrep = 10,
    num_threads = 2,
    seed = 99,
    match = "exact",
    chunk_size = 3,
    grain_size = 2,
    overwrite = TRUE,
    output = c("fit", "params", "effects")
  )

  fit_a <- data.table::fread(file.path(out_dir, paste0(prefix_a, "_fit.csv")))
  fit_b <- data.table::fread(file.path(out_dir, paste0(prefix_b, "_fit.csv")))
  data.table::setorder(fit_a, Combination)
  data.table::setorder(fit_b, Combination)
  expect_equal(fit_a, fit_b)

  params_a <- data.table::fread(file.path(out_dir, paste0(prefix_a, "_params.csv")))
  params_b <- data.table::fread(file.path(out_dir, paste0(prefix_b, "_params.csv")))
  data.table::setorder(params_a, Combination, label, lhs, rhs)
  data.table::setorder(params_b, Combination, label, lhs, rhs)
  expect_equal(params_a, params_b)

  effects_a <- data.table::fread(file.path(out_dir, paste0(prefix_a, "_effects.csv")))
  effects_b <- data.table::fread(file.path(out_dir, paste0(prefix_b, "_effects.csv")))
  data.table::setorder(effects_a, Combination, effect)
  data.table::setorder(effects_b, Combination, effect)
  expect_equal(effects_a, effects_b)
})

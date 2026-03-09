test_that("overwrite=FALSE prevents clobbering output_file", {
  set.seed(42)
  n <- 60

  test_data <- data.table::data.table(
    EXP1 = stats::rbinom(n, 1, 0.5),
    MED1 = stats::rnorm(n),
    OUT1 = stats::rnorm(n)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")
  writeLines("sentinel", output_csv)

  expect_error(
    mediation_analysis(
      data = test_data,
      columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
      nrep = 10,
      output_file = output_csv,
      num_threads = 1,
      seed = 1,
      overwrite = FALSE
    ),
    regexp = "overwrite"
  )

  expect_equal(readLines(output_csv, warn = FALSE), "sentinel")
})

test_that("excel_safe_csv prefixes formula-like Combination fields", {
  set.seed(42)
  n <- 80

  dt <- data.table::data.table(
    x = stats::rbinom(n, 1, 0.5),
    m = stats::rnorm(n),
    y = stats::rnorm(n)
  )
  data.table::setnames(dt, c("=EXP1", "MED1", "OUT1"))

  output_csv <- withr::local_tempfile(fileext = ".csv")
  mediation_analysis(
    data = dt,
    columns = list(exposure = c("=EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 10,
    output_file = output_csv,
    num_threads = 1,
    seed = 1,
    excel_safe_csv = TRUE
  )

  results <- data.table::fread(output_csv)
  expect_equal(nrow(results), 1)
  expect_equal(results$Combination, "'=EXP1_MED1_OUT1")
})

test_that("Poisson mediator with extreme treat.value does not hang", {
  set.seed(42)
  n <- 80

  x <- stats::rbinom(n, 1, 0.5)
  m <- stats::rpois(n, lambda = exp(0.1 + 0.4 * x))
  y <- 0.2 + 0.1 * m + 0.3 * x + stats::rnorm(n)

  dt <- data.table::data.table(EXP1 = x, MED1 = m, OUT1 = y)
  output_csv <- withr::local_tempfile(fileext = ".csv")

  mediation_analysis(
    data = dt,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 5,
    output_file = output_csv,
    num_threads = 1,
    pert = "asymptotic",
    seed = 1,
    mediator.family = "poisson",
    outcome.family = "gaussian",
    treat.value = 1e6,
    control.value = 0
  )

  results <- data.table::fread(output_csv)
  expect_equal(nrow(results), 1)
})

test_that("absurdly large nrep errors early", {
  set.seed(42)
  n <- 60

  test_data <- data.table::data.table(
    EXP1 = stats::rbinom(n, 1, 0.5),
    MED1 = stats::rnorm(n),
    OUT1 = stats::rnorm(n)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")

  expect_error(
    mediation_analysis(
      data = test_data,
      columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
      nrep = .Machine$integer.max,
      output_file = output_csv,
      num_threads = 1,
      seed = 1
    ),
    regexp = "nrep is too large"
  )
})

test_that("failure_mode='error' stops on failing combinations", {
  set.seed(11)
  n <- 60
  dt <- data.table::data.table(
    EXP1 = rep(0, n),
    MED1 = stats::rnorm(n),
    OUT1 = stats::rnorm(n)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")
  expect_error(
    mediation_analysis(
      data = dt,
      columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
      nrep = 10,
      output_file = output_csv,
      num_threads = 1,
      seed = 1,
      failure_mode = "error"
    ),
    regexp = "failed"
  )
})

test_that("include_failure_reason annotates failed mediation rows", {
  set.seed(12)
  n <- 100
  dt <- data.table::data.table(
    EXP1 = stats::rbinom(n, 1, 0.5),
    MED_good = stats::rnorm(n),
    MED_bad = rep(0, n),
    OUT1 = stats::rnorm(n)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")
  mediation_analysis(
    data = dt,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 20,
    output_file = output_csv,
    num_threads = 1,
    seed = 2,
    include_failure_reason = TRUE
  )

  results <- data.table::fread(output_csv)
  expect_true("failure_reason" %in% names(results))
  expect_equal(nrow(results), 2)

  failed <- results[is.na(d0_estimate)]
  succeeded <- results[!is.na(d0_estimate)]

  expect_true(nrow(failed) >= 1)
  expect_true(all(!is.na(failed$failure_reason)))
  expect_true(all(is.na(succeeded$failure_reason)))
})

test_that("default rng_mode with seed is equivalent to explicit rng_mode='fast'", {
  set.seed(13)
  n <- 80
  dt <- data.table::data.table(
    EXP1 = stats::rbinom(n, 1, 0.5),
    MED1 = stats::rnorm(n),
    OUT1 = stats::rnorm(n)
  )

  out_default <- withr::local_tempfile(fileext = ".csv")
  out_fast <- withr::local_tempfile(fileext = ".csv")

  mediation_analysis(
    data = dt,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 30,
    output_file = out_default,
    num_threads = 1,
    seed = 7
  )

  mediation_analysis(
    data = dt,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 30,
    output_file = out_fast,
    num_threads = 1,
    seed = 7,
    rng_mode = "fast"
  )

  res_default <- data.table::fread(out_default)
  res_fast <- data.table::fread(out_fast)
  expect_equal(res_default, res_fast, tolerance = 1e-12)
})

test_that("match_mediation is deprecated but overrides rng_mode", {
  set.seed(14)
  n <- 120
  dt <- data.table::data.table(
    EXP1 = stats::rbinom(n, 1, 0.5),
    MED1 = stats::rnorm(n),
    OUT1 = stats::rnorm(n)
  )

  out_alias <- withr::local_tempfile(fileext = ".csv")
  out_mode <- withr::local_tempfile(fileext = ".csv")

  expect_warning(
    mediation_analysis(
      data = dt,
      columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
      nrep = 25,
      output_file = out_alias,
      num_threads = 1,
      seed = 9,
      rng_mode = "fast",
      match_mediation = TRUE
    ),
    regexp = "deprecated"
  )

  mediation_analysis(
    data = dt,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 25,
    output_file = out_mode,
    num_threads = 1,
    seed = 9,
    rng_mode = "mediate"
  )

  res_alias <- data.table::fread(out_alias)
  res_mode <- data.table::fread(out_mode)
  expect_equal(res_alias, res_mode, tolerance = 1e-12)
})

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


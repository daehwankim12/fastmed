test_that("rank-deficient combinations produce NA rows (not abort)", {
  set.seed(42)
  n <- 120

  test_data <- data.table::data.table(
    EXP1 = stats::rbinom(n, 1, 0.5),
    MED_good = stats::rnorm(n),
    MED_bad = rep(0, n),
    OUT1 = stats::rnorm(n)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")

  expect_error(
    mediation_analysis(
      data = test_data,
      columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
      nrep = 20,
      output_file = output_csv,
      num_threads = 2,
      seed = 1,
      chunk_size = 1
    ),
    NA
  )

  results <- data.table::fread(output_csv)
  expect_equal(nrow(results), 2)

  expect_true(any(is.na(results$d0_estimate)))
  expect_true(any(!is.na(results$d0_estimate)))
})
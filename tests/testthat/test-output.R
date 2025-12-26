test_that("legacy output format writes the legacy schema", {
  set.seed(42)
  n <- 60
  test_data <- data.table::data.table(
    EXP1 = stats::rbinom(n, 1, 0.5),
    MED1 = stats::rnorm(n),
    OUT1 = stats::rnorm(n)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")

  mediation_analysis(
    data = test_data,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 10,
    output_file = output_csv,
    num_threads = 1,
    seed = 1,
    output.format = "legacy"
  )

  results <- data.table::fread(output_csv)

  expected_names <- c(
    "Combination",
    "ACME_Mean", "ACME_2.5%", "ACME_97.5%", "ACME_p-value",
    "ADE_Mean", "ADE_2.5%", "ADE_97.5%", "ADE_p-value",
    "Total_Effect_Mean", "Total_Effect_2.5%", "Total_Effect_97.5%", "Total_Effect_p-value"
  )
  expect_equal(names(results), expected_names)
})

test_that("NA rows have NA for all effect columns", {
  set.seed(42)
  n <- 120

  test_data <- data.table::data.table(
    EXP1 = stats::rbinom(n, 1, 0.5),
    MED_good = stats::rnorm(n),
    MED_bad = rep(0, n),
    OUT1 = stats::rnorm(n)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")

  mediation_analysis(
    data = test_data,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 20,
    output_file = output_csv,
    num_threads = 2,
    seed = 1,
    chunk_size = 1
  )

  results <- data.table::fread(output_csv)
  expect_equal(nrow(results), 2)

  na_rows <- results[is.na(d0_estimate)]
  expect_true(nrow(na_rows) >= 1)

  effect_cols <- grep("(_estimate|_ci_lower|_ci_upper|_p)$", names(results), value = TRUE)
  expect_true(all(is.na(na_rows[, ..effect_cols])))
})


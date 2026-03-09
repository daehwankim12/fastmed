test_that("loop_order changes output order without changing estimates (rng_mode = 'fast')", {
  set.seed(123)
  n <- 80
  test_data <- data.table::data.table(
    EXP1 = stats::rnorm(n),
    EXP2 = stats::rnorm(n),
    MED1 = stats::rnorm(n),
    MED2 = stats::rnorm(n),
    OUT1 = stats::rnorm(n),
    OUT2 = stats::rnorm(n)
  )

  output_emo <- withr::local_tempfile(fileext = ".csv")
  output_meo <- withr::local_tempfile(fileext = ".csv")

  mediation_analysis(
    data = test_data,
    columns = list(exposure = "EXP", mediator = "MED", outcome = "OUT"),
    nrep = 25,
    output_file = output_emo,
    num_threads = 1,
    pert = "asymptotic",
    seed = 1,
    rng_mode = "fast",
    loop_order = "EMO"
  )

  mediation_analysis(
    data = test_data,
    columns = list(exposure = "EXP", mediator = "MED", outcome = "OUT"),
    nrep = 25,
    output_file = output_meo,
    num_threads = 1,
    pert = "asymptotic",
    seed = 1,
    rng_mode = "fast",
    loop_order = "MEO"
  )

  results_emo <- data.table::fread(output_emo)
  results_meo <- data.table::fread(output_meo)

  expect_equal(nrow(results_emo), 8)
  expect_equal(nrow(results_meo), 8)

  expect_false(identical(results_emo$Combination, results_meo$Combination))

  setkey(results_emo, Combination)
  setkey(results_meo, Combination)

  expect_equal(results_emo$Combination, results_meo$Combination)

  effect_cols <- grep("(_estimate|_ci_lower|_ci_upper|_p)$", names(results_emo), value = TRUE)
  for (col in effect_cols) {
    expect_equal(results_emo[[col]], results_meo[[col]], tolerance = 1e-12)
  }
})

test_that("match_mediation bootstrap tolerates failed replicates", {
  n <- 20
  x <- c(rep(0, n - 1), 1)

  set.seed(123)
  test_data <- data.table::data.table(
    EXP1 = x,
    MED1 = rnorm(n),
    OUT1 = rnorm(n)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")
  nrep <- 100

  mediation_analysis(
    data = test_data,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = nrep,
    output_file = output_csv,
    num_threads = 2,
    pert = "bootstrap",
    seed = 1,
    mediator.family = "gaussian",
    outcome.family = "gaussian"
  )

  results <- data.table::fread(output_csv)
  expect_equal(nrow(results), 1)

  for (col in c("d0_estimate", "d0_ci_lower", "d0_ci_upper", "d0_p")) {
    expect_false(is.na(results[[col]]))
  }
})

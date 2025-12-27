test_that("bootstrap Estimate is invariant to nrep (t0 semantics)", {
  data <- generate_mediation_data(
    n = 250,
    treat_family = "binary",
    mediator_family = "gaussian",
    outcome_family = "gaussian",
    seed = 101
  )

  run_boot <- function(nrep) {
    output_csv <- withr::local_tempfile(fileext = ".csv")
    mediation_analysis(
      data = data,
      columns = list(exposure = "T", mediator = "M", outcome = "Y"),
      nrep = nrep,
      output_file = output_csv,
      num_threads = 1,
      pert = "bootstrap",
      seed = 99,
      mediator.family = "gaussian",
      outcome.family = "gaussian"
    )
    data.table::fread(output_csv)
  }

  res_small <- run_boot(50)
  res_large <- run_boot(250)

  expect_equal(nrow(res_small), 1)
  expect_equal(nrow(res_large), 1)

  for (col in c("d0_estimate", "d1_estimate", "z0_estimate", "z1_estimate", "tau_estimate")) {
    expect_equal(res_small[[col]], res_large[[col]], tolerance = 1e-6)
  }
})


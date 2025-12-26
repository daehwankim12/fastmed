run_mediation_once <- function(data,
                               nrep = 40,
                               seed = 123,
                               mediator.family = "auto",
                               outcome.family = "auto",
                               replace.outcome = FALSE) {
  output_csv <- withr::local_tempfile(fileext = ".csv")
  mediation_analysis(
    data = data,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = nrep,
    output_file = output_csv,
    num_threads = 1,
    pert = "asymptotic",
    seed = seed,
    mediator.family = mediator.family,
    outcome.family = outcome.family,
    replace.outcome = replace.outcome
  )
  data.table::fread(output_csv)
}

test_that("replace.outcome changes results", {
  set.seed(42)
  n <- 200
  x <- stats::rbinom(n, 1, 0.5)
  m <- 0.2 + 0.6 * x + stats::rnorm(n)
  y <- 1.0 + 0.7 * m + 0.5 * x + 2.0 * m * x + stats::rnorm(n)

  test_data <- data.table::data.table(EXP1 = x, MED1 = m, OUT1 = y)

  res_no_replace <- run_mediation_once(test_data, replace.outcome = FALSE)
  res_replace <- run_mediation_once(test_data, replace.outcome = TRUE)

  expect_equal(nrow(res_no_replace), 1)
  expect_equal(nrow(res_replace), 1)

  expect_gt(
    abs(res_no_replace[["Total_Effect_Mean"]] - res_replace[["Total_Effect_Mean"]]),
    1e-6
  )
})

test_that("family auto-detection matches explicit family selection", {
  set.seed(1)
  n <- 180

  x <- stats::rbinom(n, 1, 0.5)
  m1 <- stats::rpois(n, lambda = exp(0.3 + 0.4 * x))
  p1 <- stats::plogis(-0.2 + 0.25 * m1 + 0.15 * x)
  y1 <- stats::rbinom(n, 1, p1)
  data_pois_binom <- data.table::data.table(EXP1 = x, MED1 = m1, OUT1 = y1)

  auto_1 <- run_mediation_once(data_pois_binom, mediator.family = "auto", outcome.family = "auto")
  exp_1 <- run_mediation_once(data_pois_binom, mediator.family = "poisson", outcome.family = "binomial")
  expect_equal(auto_1[["ACME_Mean"]], exp_1[["ACME_Mean"]])
  expect_equal(auto_1[["ADE_Mean"]], exp_1[["ADE_Mean"]])
  expect_equal(auto_1[["Total_Effect_Mean"]], exp_1[["Total_Effect_Mean"]])

  m2 <- stats::rbinom(n, 1, stats::plogis(0.2 + 0.7 * x))
  y2 <- stats::rpois(n, lambda = exp(0.1 + 0.2 * m2 + 0.3 * x))
  data_binom_pois <- data.table::data.table(EXP1 = x, MED1 = m2, OUT1 = y2)

  auto_2 <- run_mediation_once(data_binom_pois, mediator.family = "auto", outcome.family = "auto")
  exp_2 <- run_mediation_once(data_binom_pois, mediator.family = "binomial", outcome.family = "poisson")
  expect_equal(auto_2[["ACME_Mean"]], exp_2[["ACME_Mean"]])
  expect_equal(auto_2[["ADE_Mean"]], exp_2[["ADE_Mean"]])
  expect_equal(auto_2[["Total_Effect_Mean"]], exp_2[["Total_Effect_Mean"]])

  m3 <- 0.1 + 0.5 * x + stats::rnorm(n)
  y3 <- 1.0 + 0.3 * m3 + 0.2 * x + stats::rnorm(n)
  data_gauss_gauss <- data.table::data.table(EXP1 = x, MED1 = m3, OUT1 = y3)

  auto_3 <- run_mediation_once(data_gauss_gauss, mediator.family = "auto", outcome.family = "auto")
  exp_3 <- run_mediation_once(data_gauss_gauss, mediator.family = "gaussian", outcome.family = "gaussian")
  expect_equal(auto_3[["ACME_Mean"]], exp_3[["ACME_Mean"]])
  expect_equal(auto_3[["ADE_Mean"]], exp_3[["ADE_Mean"]])
  expect_equal(auto_3[["Total_Effect_Mean"]], exp_3[["Total_Effect_Mean"]])
})

test_that("explicit family selection overrides auto-detection", {
  set.seed(123)
  n <- 100
  x <- stats::rbinom(n, 1, 0.5)
  m <- stats::runif(n, 0, 3) + 0.1
  y <- 1 + 0.3 * m + 0.4 * x + stats::rnorm(n)
  test_data <- data.table::data.table(EXP1 = x, MED1 = m, OUT1 = y)

  expect_error(
    run_mediation_once(test_data, mediator.family = "poisson", outcome.family = "gaussian"),
    regexp = "Poisson requires integer y"
  )
})

test_that("invalid family strings error in R wrapper", {
  set.seed(123)
  n <- 50
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
      nrep = 10,
      output_file = output_csv,
      num_threads = 1,
      mediator.family = "not-a-family"
    ),
    regexp = "mediator\\.family must be one of:"
  )
})

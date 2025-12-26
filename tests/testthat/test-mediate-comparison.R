test_that("fastmed roughly matches mediation::mediate() for Gaussian/Gaussian (asymptotic)", {
  skip_if_not_installed("mediation")

  data <- generate_mediation_data(
    n = 250,
    treat_family = "binary",
    mediator_family = "gaussian",
    outcome_family = "gaussian",
    true_acme = 0.3,
    true_ade = 0.4,
    seed = 42
  )

  sims <- 400

  med_model <- stats::lm(M ~ T, data = data)
  out_model <- stats::lm(Y ~ M + T, data = data)
  mediate_args <- list(
    model.m = med_model,
    model.y = out_model,
    treat = "T",
    mediator = "M",
    sims = sims,
    boot = FALSE
  )
  mediate_formals <- names(formals(mediation::mediate))
  if ("treat.value" %in% mediate_formals) mediate_args$treat.value <- 1
  if ("control.value" %in% mediate_formals) mediate_args$control.value <- 0
  med_result <- do.call(mediation::mediate, mediate_args)

  output_csv <- withr::local_tempfile(fileext = ".csv")
  mediation_analysis(
    data = data,
    columns = list(exposure = "T", mediator = "M", outcome = "Y"),
    nrep = sims,
    output_file = output_csv,
    num_threads = 1,
    pert = "asymptotic",
    seed = 42,
    mediator.family = "gaussian",
    outcome.family = "gaussian"
  )
  fast <- data.table::fread(output_csv)
  expect_equal(nrow(fast), 1)

  z975 <- stats::qnorm(0.975)
  se_or_ci <- function(se, ci) {
    if (is.numeric(se) && length(se) == 1L && is.finite(se) && se >= 0) {
      return(se)
    }
    if (is.numeric(ci) && length(ci) == 2L && all(is.finite(ci))) {
      return((ci[2] - ci[1]) / (2 * z975))
    }
    NA_real_
  }

  tol_from <- function(se, ci, extra = 0.05, mult = 3) {
    se_val <- se_or_ci(se, ci)
    if (!is.finite(se_val)) {
      return(0.5)
    }
    mult * se_val + extra
  }

  expect_equal(fast$d0_estimate, med_result$d0, tolerance = tol_from(med_result$d0.se, med_result$d0.ci))
  expect_equal(fast$d1_estimate, med_result$d1, tolerance = tol_from(med_result$d1.se, med_result$d1.ci))
  expect_equal(fast$z0_estimate, med_result$z0, tolerance = tol_from(med_result$z0.se, med_result$z0.ci))
  expect_equal(fast$z1_estimate, med_result$z1, tolerance = tol_from(med_result$z1.se, med_result$z1.ci))
  expect_equal(fast$tau_estimate, med_result$tau.coef, tolerance = tol_from(med_result$tau.se, med_result$tau.ci))

  fast_d0_width <- fast$d0_ci_upper - fast$d0_ci_lower
  med_d0_width <- med_result$d0.ci[2] - med_result$d0.ci[1]
  expect_equal(fast_d0_width, med_d0_width, tolerance = 0.35 * med_d0_width)
})

test_that("fastmed roughly matches mediation::mediate() for Gaussian/Gaussian (bootstrap)", {
  skip_if_not_installed("mediation")

  data <- generate_mediation_data(
    n = 200,
    treat_family = "binary",
    mediator_family = "gaussian",
    outcome_family = "gaussian",
    true_acme = 0.3,
    true_ade = 0.4,
    seed = 123
  )

  sims <- 200

  med_model <- stats::lm(M ~ T, data = data)
  out_model <- stats::lm(Y ~ M + T, data = data)
  mediate_args <- list(
    model.m = med_model,
    model.y = out_model,
    treat = "T",
    mediator = "M",
    sims = sims,
    boot = TRUE
  )
  mediate_formals <- names(formals(mediation::mediate))
  if ("treat.value" %in% mediate_formals) mediate_args$treat.value <- 1
  if ("control.value" %in% mediate_formals) mediate_args$control.value <- 0
  med_result <- do.call(mediation::mediate, mediate_args)

  output_csv <- withr::local_tempfile(fileext = ".csv")
  mediation_analysis(
    data = data,
    columns = list(exposure = "T", mediator = "M", outcome = "Y"),
    nrep = sims,
    output_file = output_csv,
    num_threads = 1,
    pert = "bootstrap",
    seed = 123,
    mediator.family = "gaussian",
    outcome.family = "gaussian"
  )
  fast <- data.table::fread(output_csv)
  expect_equal(nrow(fast), 1)

  z975 <- stats::qnorm(0.975)
  se_or_ci <- function(se, ci) {
    if (is.numeric(se) && length(se) == 1L && is.finite(se) && se >= 0) {
      return(se)
    }
    if (is.numeric(ci) && length(ci) == 2L && all(is.finite(ci))) {
      return((ci[2] - ci[1]) / (2 * z975))
    }
    NA_real_
  }

  tol_from <- function(se, ci, extra = 0.1, mult = 4) {
    se_val <- se_or_ci(se, ci)
    if (!is.finite(se_val)) {
      return(0.75)
    }
    mult * se_val + extra
  }

  expect_equal(fast$d0_estimate, med_result$d0, tolerance = tol_from(med_result$d0.se, med_result$d0.ci))
  expect_equal(fast$tau_estimate, med_result$tau.coef, tolerance = tol_from(med_result$tau.se, med_result$tau.ci))
})

test_that("fastmed completes across all family combinations (smoke)", {
  combos <- expand.grid(
    mediator_family = c("gaussian", "binomial", "poisson"),
    outcome_family = c("gaussian", "binomial", "poisson"),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(combos))) {
    med_fam <- combos$mediator_family[i]
    out_fam <- combos$outcome_family[i]
    data <- generate_mediation_data(
      n = 200,
      treat_family = "binary",
      mediator_family = med_fam,
      outcome_family = out_fam,
      true_acme = 0.2,
      true_ade = 0.2,
      seed = 100 + i
    )

    output_csv <- withr::local_tempfile(fileext = ".csv")
    mediation_analysis(
      data = data,
      columns = list(exposure = "T", mediator = "M", outcome = "Y"),
      nrep = 30,
      output_file = output_csv,
      num_threads = 1,
      pert = "asymptotic",
      seed = 200 + i,
      mediator.family = med_fam,
      outcome.family = out_fam
    )

    res <- data.table::fread(output_csv)
    expect_equal(nrow(res), 1)
    expect_true(is.finite(res$d0_estimate))
    expect_true(is.finite(res$z0_estimate))
    expect_true(is.finite(res$tau_estimate))
    expect_true(res$d0_ci_lower <= res$d0_estimate)
    expect_true(res$d0_ci_upper >= res$d0_estimate)
  }
})

test_that("Monte Carlo variability decreases as nrep increases", {
  data <- generate_mediation_data(n = 300, seed = 999)

  run_d0 <- function(nrep, seed) {
    output_csv <- withr::local_tempfile(fileext = ".csv")
    mediation_analysis(
      data = data,
      columns = list(exposure = "T", mediator = "M", outcome = "Y"),
      nrep = nrep,
      output_file = output_csv,
      num_threads = 1,
      pert = "asymptotic",
      seed = seed,
      mediator.family = "gaussian",
      outcome.family = "gaussian"
    )
    res <- data.table::fread(output_csv)
    res$d0_estimate[1]
  }

  seeds <- 1:12
  d0_50 <- vapply(seeds, function(s) run_d0(50, s), numeric(1))
  d0_200 <- vapply(seeds, function(s) run_d0(200, s), numeric(1))

  expect_gt(stats::var(d0_50), 0)
  expect_lt(stats::var(d0_200), stats::var(d0_50))
})

test_that("Gaussian/Gaussian CI contains the true d0 effect reasonably often", {
  true_acme <- 0.3
  true_d0 <- 0.5 * true_acme

  n_datasets <- 12
  covered <- logical(n_datasets)

  for (i in seq_len(n_datasets)) {
    data <- generate_mediation_data(
      n = 250,
      treat_family = "binary",
      mediator_family = "gaussian",
      outcome_family = "gaussian",
      true_acme = true_acme,
      true_ade = 0.4,
      seed = 500 + i
    )

    output_csv <- withr::local_tempfile(fileext = ".csv")
    mediation_analysis(
      data = data,
      columns = list(exposure = "T", mediator = "M", outcome = "Y"),
      nrep = 250,
      output_file = output_csv,
      num_threads = 1,
      pert = "asymptotic",
      seed = 800 + i,
      mediator.family = "gaussian",
      outcome.family = "gaussian"
    )

    res <- data.table::fread(output_csv)
    covered[i] <- isTRUE(res$d0_ci_lower[1] <= true_d0 && res$d0_ci_upper[1] >= true_d0)
  }

  expect_gte(mean(covered), 0.7)
})
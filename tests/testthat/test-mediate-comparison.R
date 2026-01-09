tol_from_ci <- function(ci, frac = 0.05, min_abs = 0.002) {
  ci <- unname(ci)
  width <- ci[2] - ci[1]
  max(min_abs, frac * abs(width))
}

check_effect <- function(effect, fast, med_est, med_ci, med_p, tol = 1e-6) {
  med_ci <- unname(med_ci)

  fast_est <- fast[[paste0(effect, "_estimate")]][1]
  fast_ci_lower <- fast[[paste0(effect, "_ci_lower")]][1]
  fast_ci_upper <- fast[[paste0(effect, "_ci_upper")]][1]
  fast_p <- fast[[paste0(effect, "_p")]][1]

  expect_true(abs(fast_est - med_est) <= tol, info = paste(effect, "estimate mismatch"))
  expect_true(abs(fast_ci_lower - med_ci[1]) <= tol, info = paste(effect, "CI lower mismatch"))
  expect_true(abs(fast_ci_upper - med_ci[2]) <= tol, info = paste(effect, "CI upper mismatch"))

  fast_width <- fast_ci_upper - fast_ci_lower
  med_width <- med_ci[2] - med_ci[1]
  expect_true(abs(fast_width - med_width) <= 2 * tol, info = paste(effect, "CI width mismatch"))

  expect_true(abs(fast_p - med_p) <= tol, info = paste(effect, "p-value mismatch"))
}

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

  sims <- 1000

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
  set.seed(42)
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

  check_effect("d0", fast, med_result$d0, med_result$d0.ci, med_result$d0.p)
  check_effect("d1", fast, med_result$d1, med_result$d1.ci, med_result$d1.p)
  check_effect("z0", fast, med_result$z0, med_result$z0.ci, med_result$z0.p)
  check_effect("z1", fast, med_result$z1, med_result$z1.ci, med_result$z1.p)
  check_effect("tau", fast, med_result$tau.coef, med_result$tau.ci, med_result$tau.p)
})

test_that("fastmed matches mediation::mediate() with NA omitted per combination", {
  skip_if_not_installed("mediation")

  data <- generate_mediation_data(
    n = 220,
    treat_family = "binary",
    mediator_family = "gaussian",
    outcome_family = "gaussian",
    true_acme = 0.3,
    true_ade = 0.4,
    seed = 2024
  )

  set.seed(2024)
  miss_idx <- sample.int(nrow(data), size = 12)
  data$T[miss_idx] <- NA_real_
  data$M[miss_idx] <- NA_real_

  sims <- 300

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
  set.seed(2024)
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
    seed = 2024,
    mediator.family = "gaussian",
    outcome.family = "gaussian"
  )
  fast <- data.table::fread(output_csv)
  expect_equal(nrow(fast), 1)

  check_effect("d0", fast, med_result$d0, med_result$d0.ci, med_result$d0.p)
  check_effect("d1", fast, med_result$d1, med_result$d1.ci, med_result$d1.p)
  check_effect("z0", fast, med_result$z0, med_result$z0.ci, med_result$z0.p)
  check_effect("z1", fast, med_result$z1, med_result$z1.ci, med_result$z1.p)
  check_effect("tau", fast, med_result$tau.coef, med_result$tau.ci, med_result$tau.p)
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
  set.seed(123)
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

  check_effect("d0", fast, med_result$d0, med_result$d0.ci, med_result$d0.p)
  check_effect("d1", fast, med_result$d1, med_result$d1.ci, med_result$d1.p)
  check_effect("z0", fast, med_result$z0, med_result$z0.ci, med_result$z0.p)
  check_effect("z1", fast, med_result$z1, med_result$z1.ci, med_result$z1.p)
  check_effect("tau", fast, med_result$tau.coef, med_result$tau.ci, med_result$tau.p)
})

test_that("fastmed roughly matches mediation::mediate() for Gaussian/Gaussian with non-0/1 treatment values", {
  skip_if_not_installed("mediation")

  set.seed(321)
  n <- 250
  T <- sample(c(-1, 1), n, replace = TRUE)
  M <- 1 + 0.5 * T + stats::rnorm(n)
  Y <- 2 + 0.3 * M + 0.4 * T + stats::rnorm(n)
  data <- data.frame(T = T, M = M, Y = Y)

  sims <- 400

  med_model <- stats::lm(M ~ T, data = data)
  out_model <- stats::lm(Y ~ M + T, data = data)
  mediate_args <- list(
    model.m = med_model,
    model.y = out_model,
    treat = "T",
    mediator = "M",
    sims = sims,
    boot = FALSE,
    treat.value = 1,
    control.value = -1
  )
  set.seed(321)
  med_result <- do.call(mediation::mediate, mediate_args)

  output_csv <- withr::local_tempfile(fileext = ".csv")
  mediation_analysis(
    data = data,
    columns = list(exposure = "T", mediator = "M", outcome = "Y"),
    nrep = sims,
    output_file = output_csv,
    num_threads = 1,
    pert = "asymptotic",
    seed = 321,
    mediator.family = "gaussian",
    outcome.family = "gaussian",
    treat.value = 1,
    control.value = -1
  )
  fast <- data.table::fread(output_csv)
  expect_equal(nrow(fast), 1)

  check_effect("d0", fast, med_result$d0, med_result$d0.ci, med_result$d0.p)
  check_effect("d1", fast, med_result$d1, med_result$d1.ci, med_result$d1.p)
  check_effect("z0", fast, med_result$z0, med_result$z0.ci, med_result$z0.p)
  check_effect("z1", fast, med_result$z1, med_result$z1.ci, med_result$z1.p)
  check_effect("tau", fast, med_result$tau.coef, med_result$tau.ci, med_result$tau.p)
})

test_that("fastmed roughly matches mediation::mediate() for Gaussian/Gaussian with weights (bootstrap)", {
  skip_if_not_installed("mediation")

  data <- generate_mediation_data(
    n = 250,
    treat_family = "binary",
    mediator_family = "gaussian",
    outcome_family = "gaussian",
    true_acme = 0.3,
    true_ade = 0.4,
    seed = 111
  )

  set.seed(111)
  w <- stats::runif(nrow(data), min = 0, max = 3)
  w <- w / mean(w)

  sims <- 200

  med_model <- stats::lm(M ~ T, data = data, weights = w)
  out_model <- stats::lm(Y ~ M + T, data = data, weights = w)
  mediate_args <- list(
    model.m = med_model,
    model.y = out_model,
    treat = "T",
    mediator = "M",
    sims = sims,
    boot = TRUE
  )
  set.seed(111)
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
    seed = 111,
    mediator.family = "gaussian",
    outcome.family = "gaussian",
    weights = w
  )
  fast <- data.table::fread(output_csv)
  expect_equal(nrow(fast), 1)

  check_effect("d0", fast, med_result$d0, med_result$d0.ci, med_result$d0.p)
  check_effect("d1", fast, med_result$d1, med_result$d1.ci, med_result$d1.p)
  check_effect("z0", fast, med_result$z0, med_result$z0.ci, med_result$z0.p)
  check_effect("z1", fast, med_result$z1, med_result$z1.ci, med_result$z1.p)
  check_effect("tau", fast, med_result$tau.coef, med_result$tau.ci, med_result$tau.p)
})

test_that("fastmed roughly matches mediation::mediate() for Gaussian/Binomial (asymptotic)", {
  skip_if_not_installed("mediation")

  data <- generate_mediation_data(
    n = 300,
    treat_family = "binary",
    mediator_family = "gaussian",
    outcome_family = "binomial",
    true_acme = 0.3,
    true_ade = 0.4,
    seed = 10
  )

  sims <- 1000

  med_model <- stats::lm(M ~ T, data = data)
  out_model <- stats::glm(Y ~ M + T, family = stats::binomial(), data = data)

  set.seed(42)
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
    outcome.family = "binomial"
  )
  fast <- data.table::fread(output_csv)
  expect_equal(nrow(fast), 1)

  check_effect("d0", fast, med_result$d0, med_result$d0.ci, med_result$d0.p)
  check_effect("d1", fast, med_result$d1, med_result$d1.ci, med_result$d1.p)
  check_effect("z0", fast, med_result$z0, med_result$z0.ci, med_result$z0.p)
  check_effect("z1", fast, med_result$z1, med_result$z1.ci, med_result$z1.p)
  check_effect("tau", fast, med_result$tau.coef, med_result$tau.ci, med_result$tau.p)
})

test_that("fastmed roughly matches mediation::mediate() for Binomial/Gaussian (asymptotic)", {
  skip_if_not_installed("mediation")

  data <- generate_mediation_data(
    n = 300,
    treat_family = "binary",
    mediator_family = "binomial",
    outcome_family = "gaussian",
    true_acme = 0.3,
    true_ade = 0.4,
    seed = 11
  )

  sims <- 1000

  med_model <- stats::glm(M ~ T, family = stats::binomial(), data = data)
  out_model <- stats::lm(Y ~ M + T, data = data)

  set.seed(42)
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
    mediator.family = "binomial",
    outcome.family = "gaussian"
  )
  fast <- data.table::fread(output_csv)
  expect_equal(nrow(fast), 1)

  check_effect("d0", fast, med_result$d0, med_result$d0.ci, med_result$d0.p)
  check_effect("d1", fast, med_result$d1, med_result$d1.ci, med_result$d1.p)
  check_effect("z0", fast, med_result$z0, med_result$z0.ci, med_result$z0.p)
  check_effect("z1", fast, med_result$z1, med_result$z1.ci, med_result$z1.p)
  check_effect("tau", fast, med_result$tau.coef, med_result$tau.ci, med_result$tau.p)
})

test_that("fastmed roughly matches mediation::mediate() for Binomial/Binomial (asymptotic)", {
  skip_if_not_installed("mediation")

  data <- generate_mediation_data(
    n = 400,
    treat_family = "binary",
    mediator_family = "binomial",
    outcome_family = "binomial",
    true_acme = 0.3,
    true_ade = 0.4,
    seed = 12
  )

  sims <- 1000

  med_model <- stats::glm(M ~ T, family = stats::binomial(), data = data)
  out_model <- stats::glm(Y ~ M + T, family = stats::binomial(), data = data)

  set.seed(42)
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
    mediator.family = "binomial",
    outcome.family = "binomial"
  )
  fast <- data.table::fread(output_csv)
  expect_equal(nrow(fast), 1)

  check_effect("d0", fast, med_result$d0, med_result$d0.ci, med_result$d0.p)
  check_effect("d1", fast, med_result$d1, med_result$d1.ci, med_result$d1.p)
  check_effect("z0", fast, med_result$z0, med_result$z0.ci, med_result$z0.p)
  check_effect("z1", fast, med_result$z1, med_result$z1.ci, med_result$z1.p)
  check_effect("tau", fast, med_result$tau.coef, med_result$tau.ci, med_result$tau.p)
})

test_that("fastmed roughly matches mediation::mediate() for Poisson/Poisson (asymptotic)", {
  skip_if_not_installed("mediation")

  data <- generate_mediation_data(
    n = 400,
    treat_family = "binary",
    mediator_family = "poisson",
    outcome_family = "poisson",
    true_acme = 0.3,
    true_ade = 0.4,
    seed = 13
  )

  sims <- 1000

  med_model <- stats::glm(M ~ T, family = stats::poisson(), data = data)
  out_model <- stats::glm(Y ~ M + T, family = stats::poisson(), data = data)

  set.seed(42)
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
    mediator.family = "poisson",
    outcome.family = "poisson"
  )
  fast <- data.table::fread(output_csv)
  expect_equal(nrow(fast), 1)

  check_effect("d0", fast, med_result$d0, med_result$d0.ci, med_result$d0.p)
  check_effect("d1", fast, med_result$d1, med_result$d1.ci, med_result$d1.p)
  check_effect("z0", fast, med_result$z0, med_result$z0.ci, med_result$z0.p)
  check_effect("z1", fast, med_result$z1, med_result$z1.ci, med_result$z1.p)
  check_effect("tau", fast, med_result$tau.coef, med_result$tau.ci, med_result$tau.p)
})

test_that("fastmed matches mediation::mediate() for Poisson/Poisson with large mediator lambda (asymptotic)", {
  skip_if_not_installed("mediation")

  set.seed(913)
  n <- 400
  T <- stats::rbinom(n, 1, 0.5)
  lambda_M <- exp(3.0 + 0.4 * T) # >= ~20
  M <- stats::rpois(n, lambda_M)
  lambda_Y <- exp(0.2 + 0.03 * M + 0.4 * T)
  Y <- stats::rpois(n, lambda_Y)
  data <- data.frame(T = T, M = M, Y = Y)

  sims <- 200

  med_model <- stats::glm(M ~ T, family = stats::poisson(), data = data)
  out_model <- stats::glm(Y ~ M + T, family = stats::poisson(), data = data)

  set.seed(42)
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
    mediator.family = "poisson",
    outcome.family = "poisson"
  )
  fast <- data.table::fread(output_csv)
  expect_equal(nrow(fast), 1)

  check_effect("d0", fast, med_result$d0, med_result$d0.ci, med_result$d0.p)
  check_effect("d1", fast, med_result$d1, med_result$d1.ci, med_result$d1.p)
  check_effect("z0", fast, med_result$z0, med_result$z0.ci, med_result$z0.p)
  check_effect("z1", fast, med_result$z1, med_result$z1.ci, med_result$z1.p)
  check_effect("tau", fast, med_result$tau.coef, med_result$tau.ci, med_result$tau.p)
})

test_that("fastmed roughly matches mediation::mediate() for Gaussian/Poisson (asymptotic)", {
  skip_if_not_installed("mediation")

  data <- generate_mediation_data(
    n = 400,
    treat_family = "binary",
    mediator_family = "gaussian",
    outcome_family = "poisson",
    true_acme = 0.3,
    true_ade = 0.4,
    seed = 14
  )

  sims <- 500

  med_model <- stats::lm(M ~ T, data = data)
  out_model <- stats::glm(Y ~ M + T, family = stats::poisson(), data = data)

  set.seed(42)
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
    outcome.family = "poisson"
  )
  fast <- data.table::fread(output_csv)
  expect_equal(nrow(fast), 1)

  check_effect("d0", fast, med_result$d0, med_result$d0.ci, med_result$d0.p)
  check_effect("d1", fast, med_result$d1, med_result$d1.ci, med_result$d1.p)
  check_effect("z0", fast, med_result$z0, med_result$z0.ci, med_result$z0.p)
  check_effect("z1", fast, med_result$z1, med_result$z1.ci, med_result$z1.p)
  check_effect("tau", fast, med_result$tau.coef, med_result$tau.ci, med_result$tau.p)
})

test_that("fastmed roughly matches mediation::mediate() for Poisson/Gaussian (asymptotic)", {
  skip_if_not_installed("mediation")

  data <- generate_mediation_data(
    n = 400,
    treat_family = "binary",
    mediator_family = "poisson",
    outcome_family = "gaussian",
    true_acme = 0.3,
    true_ade = 0.4,
    seed = 15
  )

  sims <- 500

  med_model <- stats::glm(M ~ T, family = stats::poisson(), data = data)
  out_model <- stats::lm(Y ~ M + T, data = data)

  set.seed(42)
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
    mediator.family = "poisson",
    outcome.family = "gaussian"
  )
  fast <- data.table::fread(output_csv)
  expect_equal(nrow(fast), 1)

  check_effect("d0", fast, med_result$d0, med_result$d0.ci, med_result$d0.p)
  check_effect("d1", fast, med_result$d1, med_result$d1.ci, med_result$d1.p)
  check_effect("z0", fast, med_result$z0, med_result$z0.ci, med_result$z0.p)
  check_effect("z1", fast, med_result$z1, med_result$z1.ci, med_result$z1.p)
  check_effect("tau", fast, med_result$tau.coef, med_result$tau.ci, med_result$tau.p)
})

test_that("fastmed roughly matches mediation::mediate() for Binomial/Poisson (asymptotic)", {
  skip_if_not_installed("mediation")

  data <- generate_mediation_data(
    n = 500,
    treat_family = "binary",
    mediator_family = "binomial",
    outcome_family = "poisson",
    true_acme = 0.3,
    true_ade = 0.4,
    seed = 16
  )

  sims <- 500

  med_model <- stats::glm(M ~ T, family = stats::binomial(), data = data)
  out_model <- stats::glm(Y ~ M + T, family = stats::poisson(), data = data)

  set.seed(42)
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
    mediator.family = "binomial",
    outcome.family = "poisson"
  )
  fast <- data.table::fread(output_csv)
  expect_equal(nrow(fast), 1)

  check_effect("d0", fast, med_result$d0, med_result$d0.ci, med_result$d0.p)
  check_effect("d1", fast, med_result$d1, med_result$d1.ci, med_result$d1.p)
  check_effect("z0", fast, med_result$z0, med_result$z0.ci, med_result$z0.p)
  check_effect("z1", fast, med_result$z1, med_result$z1.ci, med_result$z1.p)
  check_effect("tau", fast, med_result$tau.coef, med_result$tau.ci, med_result$tau.p)
})

test_that("fastmed roughly matches mediation::mediate() for Poisson/Binomial (asymptotic)", {
  skip_if_not_installed("mediation")

  data <- generate_mediation_data(
    n = 500,
    treat_family = "binary",
    mediator_family = "poisson",
    outcome_family = "binomial",
    true_acme = 0.3,
    true_ade = 0.4,
    seed = 17
  )

  sims <- 500

  med_model <- stats::glm(M ~ T, family = stats::poisson(), data = data)
  out_model <- stats::glm(Y ~ M + T, family = stats::binomial(), data = data)

  set.seed(42)
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
    mediator.family = "poisson",
    outcome.family = "binomial"
  )
  fast <- data.table::fread(output_csv)
  expect_equal(nrow(fast), 1)

  check_effect("d0", fast, med_result$d0, med_result$d0.ci, med_result$d0.p)
  check_effect("d1", fast, med_result$d1, med_result$d1.ci, med_result$d1.p)
  check_effect("z0", fast, med_result$z0, med_result$z0.ci, med_result$z0.p)
  check_effect("z1", fast, med_result$z1, med_result$z1.ci, med_result$z1.p)
  check_effect("tau", fast, med_result$tau.coef, med_result$tau.ci, med_result$tau.p)
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

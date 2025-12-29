test_that("mediation_analysis works correctly with default prefixes", {
  # Create example data
  set.seed(123)
  test_data_default <- data.table::data.table(
    CECUM1 = rnorm(50),
    CECUM2 = rnorm(50),
    SERUM1 = rnorm(50),
    SERUM2 = rnorm(50),
    CORTEX1 = rnorm(50),
    CORTEX2 = rnorm(50)
  )

  # Use local_tempfile to create a temporary output CSV file
  output_csv <- withr::local_tempfile(fileext = ".csv")

  # Run mediation_analysis function
  mediation_analysis(
    data = test_data_default,
    columns = list(exposure = c("CECUM"), mediator = c("SERUM"), outcome = c("CORTEX")),
    nrep = 100, # Reduced for faster testing
    output_file = output_csv, # Use the temporary file name
    num_threads = 1 # Adjust thread count for testing environment
  )

  # Read results
  results <- data.table::fread(output_csv)

  # Expected result: 2 (CECUM) * 2 (SERUM) * 2 (CORTEX) = 8 combinations
  expect_equal(nrow(results), 8)
})

test_that("mediation_analysis works correctly with custom prefixes", {
  # Create example data
  set.seed(123)
  test_data_custom <- data.table::data.table(
    EXP1 = rnorm(50),
    EXP2 = rnorm(50),
    MED1 = rnorm(50),
    MED2 = rnorm(50),
    OUT1 = rnorm(50),
    OUT2 = rnorm(50)
  )

  # Use local_tempfile to create a temporary output CSV file
  output_csv <- withr::local_tempfile(fileext = ".csv")

  # Run mediation_analysis function
  mediation_analysis(
    data = test_data_custom,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 100, # Reduced for faster testing
    output_file = output_csv, # Use the temporary file name
    num_threads = 1 # Adjust thread count for testing environment
  )

  # Read results
  results <- data.table::fread(output_csv)

  # Expected result: 2 (EXP) * 2 (MED) * 2 (OUT) = 8 combinations
  expect_equal(nrow(results), 8)
})

test_that("mediation_analysis handles large datasets correctly", {
  # Create large example data
  set.seed(123)
  test_data_large <- data.table::data.table(
    CECUM1 = rnorm(1000),
    CECUM2 = rnorm(1000),
    SERUM1 = rnorm(1000),
    SERUM2 = rnorm(1000),
    CORTEX1 = rnorm(1000),
    CORTEX2 = rnorm(1000)
  )

  # Use local_tempfile to create a temporary output CSV file
  output_csv <- withr::local_tempfile(fileext = ".csv")

  # Run mediation_analysis function
  mediation_analysis(
    data = test_data_large,
    columns = list(exposure = c("CECUM"), mediator = c("SERUM"), outcome = c("CORTEX")),
    nrep = 100, # Reduced for faster testing
    output_file = output_csv, # Use the temporary file name
    num_threads = 1 # Adjust thread count for testing environment
  )

  # Read results
  results <- data.table::fread(output_csv)

  # Expected result: 2 (CECUM) * 2 (SERUM) * 2 (CORTEX) = 8 combinations
  expect_equal(nrow(results), 8)
})

test_that("mediation_analysis works correctly with multiple threads", {
  # Create large example data
  set.seed(123)
  test_data_large <- data.table::data.table(
    CECUM1 = rnorm(1000),
    CECUM2 = rnorm(1000),
    SERUM1 = rnorm(1000),
    SERUM2 = rnorm(1000),
    CORTEX1 = rnorm(1000),
    CORTEX2 = rnorm(1000)
  )

  # Use local_tempfile to create a temporary output CSV file
  output_csv <- withr::local_tempfile(fileext = ".csv")

  # Run mediation_analysis function
  mediation_analysis(
    data = test_data_large,
    columns = list(exposure = c("CECUM"), mediator = c("SERUM"), outcome = c("CORTEX")),
    nrep = 100, # Reduced for faster testing
    output_file = output_csv, # Use the temporary file name
    num_threads = 4 # Adjust thread count for testing environment
  )

  # Read results
  results <- data.table::fread(output_csv)

  # Expected result: 2 (CECUM) * 2 (SERUM) * 2 (CORTEX) = 8 combinations
  expect_equal(nrow(results), 8)
})

test_that("asymptotic method errors when n <= p", {
  test_data_small <- data.table::data.table(
    EXP1 = c(0, 1),
    MED1 = c(0.1, 0.2),
    OUT1 = c(0.3, 0.4)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")

  expect_error(
    mediation_analysis(
      data = test_data_small,
      columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
      nrep = 10,
      output_file = output_csv,
      num_threads = 1,
      pert = "asymptotic"
    ),
    regexp = "Insufficient observations for mediator model"
  )
})

test_that("singular design matrix yields NA output row", {
  set.seed(123)
  test_data_singular <- data.table::data.table(
    EXP1 = rep(0, 20),
    MED1 = rnorm(20),
    OUT1 = rnorm(20)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")

  expect_error(
    mediation_analysis(
      data = test_data_singular,
      columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
      nrep = 10,
      output_file = output_csv,
      num_threads = 1,
      pert = "asymptotic"
    ),
    NA
  )

  results <- data.table::fread(output_csv)
  expect_equal(nrow(results), 1)
  expect_equal(results$Combination, "EXP1_MED1_OUT1")
  expect_true(is.na(results$d0_estimate))
})

test_that("p-values use mediate sign test", {
  set.seed(123)
  n <- 200
  x <- rep(c(0, 1), each = n / 2)
  m <- 100 + 50 * x + rnorm(n, sd = 0.01)
  y <- 10 + 20 * m + 10 * x + rnorm(n, sd = 0.01)

  test_data_strong <- data.table::data.table(EXP1 = x, MED1 = m, OUT1 = y)
  output_csv <- withr::local_tempfile(fileext = ".csv")

  nrep <- 100
  mediation_analysis(
    data = test_data_strong,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = nrep,
    output_file = output_csv,
    num_threads = 1,
    pert = "asymptotic",
    seed = 1
  )

  results <- data.table::fread(output_csv)
  expect_equal(nrow(results), 1)

  p_value <- results[["d0_p"]]
  expect_true(is.numeric(p_value))
  expect_equal(length(p_value), 1)
  expect_gte(p_value, 0)
  expect_lte(p_value, 1)
  expect_lte(p_value, 0.05)
})

test_that("NA data is rejected with clear error", {
  test_data <- data.table::data.table(
    EXP1 = c(rnorm(49), NA_real_),
    MED1 = rnorm(50),
    OUT1 = rnorm(50)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")

  expect_error(
    mediation_analysis(
      data = test_data,
      columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
      nrep = 10,
      output_file = output_csv,
      num_threads = 1
    ),
    regexp = "missing values"
  )
})

test_that("seed parameter produces reproducible results", {
  set.seed(123)
  test_data <- data.table::data.table(
    EXP1 = rnorm(80),
    MED1 = rnorm(80),
    OUT1 = rnorm(80)
  )

  output_csv1 <- withr::local_tempfile(fileext = ".csv")
  output_csv2 <- withr::local_tempfile(fileext = ".csv")

  mediation_analysis(
    data = test_data,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 30,
    output_file = output_csv1,
    num_threads = 1,
    seed = 42
  )

  mediation_analysis(
    data = test_data,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 30,
    output_file = output_csv2,
    num_threads = 1,
    seed = 42
  )

  results1 <- data.table::fread(output_csv1)
  results2 <- data.table::fread(output_csv2)

  numeric_cols <- names(results1)[vapply(results1, is.numeric, logical(1))]
  for (col in numeric_cols) {
    expect_equal(results1[[col]], results2[[col]], tolerance = 1e-10)
  }
  expect_equal(results1$Combination, results2$Combination)
})

test_that("results identical across thread counts with same seed", {
  set.seed(123)
  test_data <- data.table::data.table(
    EXP1 = rnorm(80),
    EXP2 = rnorm(80),
    MED1 = rnorm(80),
    OUT1 = rnorm(80)
  )

  output_csv1 <- withr::local_tempfile(fileext = ".csv")
  output_csv2 <- withr::local_tempfile(fileext = ".csv")

  mediation_analysis(
    data = test_data,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 30,
    output_file = output_csv1,
    num_threads = 1,
    seed = 7
  )

  mediation_analysis(
    data = test_data,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 30,
    output_file = output_csv2,
    num_threads = 4,
    seed = 7
  )

  results1 <- data.table::fread(output_csv1)
  results2 <- data.table::fread(output_csv2)

  numeric_cols <- names(results1)[vapply(results1, is.numeric, logical(1))]
  for (col in numeric_cols) {
    expect_equal(results1[[col]], results2[[col]], tolerance = 1e-10)
  }
  expect_equal(results1$Combination, results2$Combination)
})

test_that("CSV escaping handles special characters in Combination", {
  set.seed(123)
  test_data <- data.table::data.table(
    `EXP,1` = rnorm(50),
    MED1 = rnorm(50),
    OUT1 = rnorm(50)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")

  mediation_analysis(
    data = test_data,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 10,
    output_file = output_csv,
    num_threads = 1
  )

  results <- data.table::fread(output_csv)
  expect_true("Combination" %in% names(results))
  expect_equal(nrow(results), 1)
})

test_that("non-numeric columns are rejected", {
  test_data <- data.table::data.table(
    EXP1 = rnorm(50),
    MED1 = rnorm(50),
    OUT1 = rnorm(50),
    Category = sample(c("A", "B"), 50, replace = TRUE)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")

  expect_error(
    mediation_analysis(
      data = test_data,
      columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
      nrep = 10,
      output_file = output_csv,
      num_threads = 1
    ),
    regexp = "Non-numeric columns"
  )
})

test_that("missing column type is rejected", {
  test_data <- data.table::data.table(
    EXP1 = rnorm(50),
    MED1 = rnorm(50),
    OUT1 = rnorm(50)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")

  expect_error(
    mediation_analysis(
      data = test_data,
      columns = list(exposure = c("EXP"), mediator = c("MED")),
      nrep = 10,
      output_file = output_csv,
      num_threads = 1
    ),
    regexp = "3 elements"
  )
})

test_that("no matching columns rejected with clear error", {
  test_data <- data.table::data.table(
    EXP1 = rnorm(50),
    MED1 = rnorm(50),
    OUT1 = rnorm(50)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")

  expect_error(
    mediation_analysis(
      data = test_data,
      columns = list(exposure = c("WRONG"), mediator = c("MED"), outcome = c("OUT")),
      nrep = 10,
      output_file = output_csv,
      num_threads = 1
    ),
    regexp = "No columns found matching exposure prefixes"
  )
})

test_that("startsWith prefix matching works", {
  set.seed(123)
  test_data <- data.table::data.table(
    Exposure_A1 = rnorm(50),
    Exposure_A2 = rnorm(50),
    ExposureBad = rnorm(50),
    Mediator_X1 = rnorm(50),
    Outcome_Y1 = rnorm(50)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")

  mediation_analysis(
    data = test_data,
    columns = list(exposure = c("Exposure_A"), mediator = c("Mediator_X"), outcome = c("Outcome_Y")),
    nrep = 10,
    output_file = output_csv,
    num_threads = 1
  )

  results <- data.table::fread(output_csv)
  expect_equal(nrow(results), 2)
})

test_that("chunked analysis matches non-chunked (same seed)", {
  set.seed(123)
  test_data <- data.table::data.table(
    EXP1 = rnorm(60),
    EXP2 = rnorm(60),
    MED1 = rnorm(60),
    MED2 = rnorm(60),
    OUT1 = rnorm(60),
    OUT2 = rnorm(60)
  )

  output_csv1 <- withr::local_tempfile(fileext = ".csv")
  output_csv2 <- withr::local_tempfile(fileext = ".csv")

  mediation_analysis(
    data = test_data,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 20,
    output_file = output_csv1,
    num_threads = 2,
    seed = 99,
    chunk_size = 100
  )

  mediation_analysis(
    data = test_data,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 20,
    output_file = output_csv2,
    num_threads = 2,
    seed = 99,
    chunk_size = 2
  )

  results1 <- data.table::fread(output_csv1)
  results2 <- data.table::fread(output_csv2)

  numeric_cols <- names(results1)[vapply(results1, is.numeric, logical(1))]
  for (col in numeric_cols) {
    expect_equal(results1[[col]], results2[[col]], tolerance = 1e-10)
  }
  expect_equal(results1$Combination, results2$Combination)
})

test_that("writes a valid CSV with a single header", {
  set.seed(123)
  test_data <- data.table::data.table(
    EXP1 = rnorm(50),
    EXP2 = rnorm(50),
    EXP3 = rnorm(50),
    MED1 = rnorm(50),
    OUT1 = rnorm(50)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")

  mediation_analysis(
    data = test_data,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 10,
    output_file = output_csv,
    num_threads = 2,
    chunk_size = 1
  )

  results <- data.table::fread(output_csv)
  expect_equal(nrow(results), 3)

  lines <- readLines(output_csv)
  header_count <- sum(grepl("^Combination,", lines))
  expect_equal(header_count, 1)
})

test_that("asymptotic sigma matches lm residual SE", {
  set.seed(123)
  n <- 60
  x <- rnorm(n)
  m <- 0.5 + 2 * x + rnorm(n, sd = 0.3)
  y <- -1 + 3 * m + 0.7 * x + rnorm(n, sd = 0.4)

  sigmas <- fastmed:::fastmed_test_ols_sigmas_cpp(x, m, y)

  expect_equal(sigmas$df_med, n - 2)
  expect_equal(sigmas$df_out, n - 3)

  sigma_med_r <- sigma(stats::lm(m ~ x))
  sigma_out_r <- sigma(stats::lm(y ~ m + x))

  expect_equal(sigmas$sigma_med, sigma_med_r, tolerance = 1e-10)
  expect_equal(sigmas$sigma_out, sigma_out_r, tolerance = 1e-10)
})

test_that("p-value deterministic fixtures match definition", {
  p0 <- fastmed:::fastmed_test_p_value_cpp(numeric(0))
  expect_equal(p0, 1.0)

  n <- 100
  p_all_pos <- fastmed:::fastmed_test_p_value_cpp(rep(1, n))
  expect_equal(p_all_pos, 0.0, tolerance = 1e-12)

  p_balanced <- fastmed:::fastmed_test_p_value_cpp(c(rep(-1, n / 2), rep(1, n / 2)))
  expect_equal(p_balanced, 1.0)

  p_zeros_split <- fastmed:::fastmed_test_p_value_cpp(c(0, 0, 0, 1))
  expect_equal(p_zeros_split, 0.0, tolerance = 1e-12)

  p_all_zero <- fastmed:::fastmed_test_p_value_cpp(rep(0, n))
  expect_equal(p_all_zero, 1.0, tolerance = 1e-12)
})

test_that("p-value matches mediation:::pval when available", {
  skip_if_not_installed("mediation")

  pval <- getFromNamespace("pval", "mediation")
  sims <- c(rep(1, 100), -1)
  est <- mean(sims)

  expected <- {
    fmls <- names(formals(pval))
    if ("c0" %in% fmls) {
      pval(sims, c0 = est)
    } else if ("xhat" %in% fmls) {
      pval(sims, xhat = est)
    } else {
      pval(sims)
    }
  }

  got <- fastmed:::fastmed_test_p_value_cpp(sims)
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("statistics summary matches R reference implementation", {
  quantile_ref <- function(samples, prob) {
    as.numeric(stats::quantile(samples, prob, type = 7, names = FALSE))
  }

  p_value_ref <- function(samples) {
    n <- length(samples)
    if (n == 0) {
      return(1.0)
    }

    est <- mean(samples)
    if (est == 0) {
      return(1.0)
    }

    pos <- sum(samples > 0)
    neg <- sum(samples < 0)

    p <- 2 * min(pos, neg) / n
    min(p, 1.0)
  }

  samples <- c(seq(-2, 2, length.out = 101), 0, 0, 0)
  stats_cpp <- fastmed:::fastmed_test_calculate_statistics_cpp(samples)

  expect_equal(stats_cpp$mean, mean(samples), tolerance = 1e-12)
  expect_equal(stats_cpp$percentile_2_5, quantile_ref(samples, 0.025), tolerance = 1e-12)
  expect_equal(stats_cpp$percentile_97_5, quantile_ref(samples, 0.975), tolerance = 1e-12)
  expect_equal(stats_cpp$p_value, p_value_ref(samples), tolerance = 1e-12)
})

test_that("p-value/CI consistency holds on symmetric fixtures (scoped)", {
  stats_pos <- fastmed:::fastmed_test_calculate_statistics_cpp(rep(c(0.9, 1.1), 50))
  ci_excludes_zero_pos <- stats_pos$percentile_2_5 > 0 || stats_pos$percentile_97_5 < 0
  expect_equal(stats_pos$p_value < 0.05, ci_excludes_zero_pos)

  stats_zero <- fastmed:::fastmed_test_calculate_statistics_cpp(rep(c(-1, 1), 50))
  ci_excludes_zero_zero <- stats_zero$percentile_2_5 > 0 || stats_zero$percentile_97_5 < 0
  expect_equal(stats_zero$p_value < 0.05, ci_excludes_zero_zero)
})

test_that("bootstrap retries advance RNG state (no reseed on retry)", {
  n <- 50
  draws <- fastmed:::fastmed_test_two_bootstrap_samples(
    n = n,
    base_seed = 123,
    global_combination_idx = 0,
    rep_idx = 0
  )

  expect_equal(dim(draws), c(n, 2))
  expect_true(all(draws >= 0))
  expect_true(all(draws < n))
  expect_false(all(draws[, 1] == draws[, 2]))

  draws2 <- fastmed:::fastmed_test_two_bootstrap_samples(
    n = n,
    base_seed = 123,
    global_combination_idx = 0,
    rep_idx = 0
  )
  expect_equal(draws, draws2)
})

test_that("output schema and row ordering are deterministic", {
  set.seed(123)
  test_data <- data.table::data.table(
    EXP2 = rnorm(50),
    EXP1 = rnorm(50),
    MEDB = rnorm(50),
    MEDA = rnorm(50),
    OUT2 = rnorm(50),
    OUT1 = rnorm(50)
  )

  output_csv <- withr::local_tempfile(fileext = ".csv")

  mediation_analysis(
    data = test_data,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 10,
    output_file = output_csv,
    num_threads = 4,
    seed = 1
  )

  results <- data.table::fread(output_csv)
  expect_equal(nrow(results), 8)

  expected_names <- c(
    "Combination",
    "d0_estimate", "d0_ci_lower", "d0_ci_upper", "d0_p",
    "d1_estimate", "d1_ci_lower", "d1_ci_upper", "d1_p",
    "z0_estimate", "z0_ci_lower", "z0_ci_upper", "z0_p",
    "z1_estimate", "z1_ci_lower", "z1_ci_upper", "z1_p",
    "tau_estimate", "tau_ci_lower", "tau_ci_upper", "tau_p"
  )
  expect_equal(names(results), expected_names)

  expected_order <- c(
    "EXP2_MEDB_OUT2",
    "EXP1_MEDB_OUT2",
    "EXP2_MEDA_OUT2",
    "EXP1_MEDA_OUT2",
    "EXP2_MEDB_OUT1",
    "EXP1_MEDB_OUT1",
    "EXP2_MEDA_OUT1",
    "EXP1_MEDA_OUT1"
  )
  expect_equal(results$Combination, expected_order)
})

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

test_that("asymptotic method errors on singular design matrix", {
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
    regexp = "Cholesky decomposition failed"
  )
})

test_that("p-values use sign test with continuity correction", {
  set.seed(123)
  n <- 200
  x <- rep(c(0, 1), each = n / 2)
  m <- 10 * x + rnorm(n, sd = 0.1)
  y <- 10 * m + rnorm(n, sd = 0.1)

  test_data_strong <- data.table::data.table(EXP1 = x, MED1 = m, OUT1 = y)
  output_csv <- withr::local_tempfile(fileext = ".csv")

  nrep <- 100
  mediation_analysis(
    data = test_data_strong,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = nrep,
    output_file = output_csv,
    num_threads = 1,
    pert = "asymptotic"
  )

  results <- data.table::fread(output_csv)
  expect_equal(nrow(results), 1)

  p_value <- results[["ACME_p-value"]]
  expect_true(is.numeric(p_value))
  expect_equal(length(p_value), 1)
  expect_gt(p_value, 0)
  expect_lt(p_value, 1)
  expect_gte(p_value, 2 / (nrep + 2))
  expect_lte(p_value, 0.2)
})

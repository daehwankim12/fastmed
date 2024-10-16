# tests/testthat/test-mediation_analysis.R
test_that("mediation_analysis works correctly with default prefixes", {
  # Create example data
  set.seed(123)
  test_data_default <- data.table(
    CECUM1 = rnorm(50),
    CECUM2 = rnorm(50),
    SERUM1 = rnorm(50),
    SERUM2 = rnorm(50),
    CORTEX1 = rnorm(50),
    CORTEX2 = rnorm(50)
  )

  # Temporary output CSV file
  output_csv <- tempfile(fileext = ".csv")

  # Run mediation_analysis function
  mediation_analysis(
    data = test_data_default,
    columns = list(exposure = c("CECUM"), mediator = c("SERUM"), outcome = c("CORTEX")),
    nrep = 100, # Reduced for faster testing
    output_file = output_csv,
    num_threads = 2 # Adjust thread count for testing environment
  )

  # Read results
  results <- fread(output_csv)

  # Expected result: 2 (CECUM) * 2 (SERUM) * 2 (CORTEX) = 8 combinations
  expect_equal(nrow(results), 8)

  # Remove temporary file
  file.remove(output_csv)
})

test_that("mediation_analysis works correctly with custom prefixes", {
  # Create example data
  set.seed(123)
  test_data_custom <- data.table(
    EXP1 = rnorm(50),
    EXP2 = rnorm(50),
    MED1 = rnorm(50),
    MED2 = rnorm(50),
    OUT1 = rnorm(50),
    OUT2 = rnorm(50)
  )

  # Temporary output CSV file
  output_csv <- tempfile(fileext = ".csv")

  # Run mediation_analysis function
  mediation_analysis(
    data = test_data_custom,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = 100, # Reduced for faster testing
    output_file = output_csv,
    num_threads = 2 # Adjust thread count for testing environment
  )

  # Read results
  results <- fread(output_csv)

  # Expected result: 2 (EXP) * 2 (MED) * 2 (OUT) = 8 combinations
  expect_equal(nrow(results), 8)

  # Remove temporary file
  file.remove(output_csv)
})

test_that("mediation_analysis handles large datasets correctly", {
  # Create large example data
  set.seed(123)
  test_data_large <- data.table(
    CECUM1 = rnorm(1000),
    CECUM2 = rnorm(1000),
    SERUM1 = rnorm(1000),
    SERUM2 = rnorm(1000),
    CORTEX1 = rnorm(1000),
    CORTEX2 = rnorm(1000)
  )

  # Temporary output CSV file
  output_csv <- tempfile(fileext = ".csv")

  # Run mediation_analysis function
  mediation_analysis(
    data = test_data_large,
    columns = list(exposure = c("CECUM"), mediator = c("SERUM"), outcome = c("CORTEX")),
    nrep = 100, # Reduced for faster testing
    output_file = output_csv,
    num_threads = 2 # Adjust thread count for testing environment
  )

  # Read results
  results <- fread(output_csv)

  # Expected result: 2 (CECUM) * 2 (SERUM) * 2 (CORTEX) = 8 combinations
  expect_equal(nrow(results), 8)

  # Remove temporary file
  file.remove(output_csv)
})

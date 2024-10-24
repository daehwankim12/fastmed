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

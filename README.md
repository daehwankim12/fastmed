# fastmed

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg) ![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/fastmed) [![R-CMD-check](https://github.com/daehwankim12/fastmed/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/daehwankim12/fastmed/actions/workflows/R-CMD-check.yaml) 

## Overview

**fastmed** is an R package designed to perform efficient and scalable mediation analysis using multiple linear regressions and bootstrap resampling. Leveraging the power of C++ through Rcpp and parallel processing with RcppParallel, `fastmed` is optimized to handle large datasets and provide real-time results output to CSV files.

## Features

-   **Efficient Mediation Analysis:** Conduct mediation analysis with multiple exposure, mediator, and outcome variables.
-   **Parallel Processing:** Utilize multiple CPU cores to accelerate computations.
-   **Bootstrap Resampling:** Estimate direct, indirect, and total effects with confidence intervals and p-values.
-   **Customizable Column Prefixes:** Easily specify prefixes to identify exposure, mediator, and outcome variables.
-   **Scalable for Large Datasets:** Handle large-scale data efficiently by processing in chunks.
-   **Real-time CSV Output:** Save results directly to CSV files during analysis.

## Installation

You can install the development version of **fastmed** from GitHub using the `devtools` package:

``` r
# Install devtools if not already installed
install.packages("devtools")

# Install fastmed from GitHub
devtools::install_github("daehwankim12/fastmed")
```

## Usage

Below is an example of how to use the `mediation_analysis` function in fastmed.

### Example

``` r
library(fastmed)
library(data.table)

# Generate example data
set.seed(123)
my_data <- data.table(
  Exposure_A1 = rnorm(1000),
  Exposure_A2 = rnorm(1000),
  Mediator_X1 = rnorm(1000),
  Mediator_X2 = rnorm(1000),
  Outcome_Y1 = rnorm(1000),
  Outcome_Y2 = rnorm(1000)
)

# Define column prefixes
columns <- list(
  exposure = "Exposure_A",
  mediator = "Mediator_X",
  outcome = "Outcome_Y"
)

# Specify output CSV file path
output_csv <- "mediation_results.csv"

# Perform mediation analysis
mediation_analysis(
  data = my_data,
  columns = columns,
  nrep = 500,            # Number of bootstrap replicates
  output_file = output_csv,
  num_threads = 4        # Number of threads for parallel processing
)

# View results
results <- fread(output_csv)
print(results)
```

### Parameters

-   `data`: A data.table or data.frame containing the dataset.
-   `columns`: A list with three named elements: exposure, mediator, and outcome. Each should be a character string specifying the prefix of the respective columns.
-   `nrep`: (Optional) Number of bootstrap replicates. Default is 1000.
-   `output_file`: Path to the output CSV file where results will be saved.
-   `num_threads`: (Optional) Number of threads for parallel processing. Defaults to the number of available cores.

### Output

The output CSV file (`output_csv` in the example) will contain the following columns for each combination of exposure, mediator, and outcome variables:

-   `combination`: Identifier for the combination (e.g., Exposure_A1_Mediator_X1_Outcome_Y1)
-   `indirect_t1_estimate`: Indirect effect estimate at t=1
-   `indirect_t1_std_err`: Standard error of the indirect effect at t=1
-   `indirect_t1_lcb`: Lower confidence bound for the indirect effect at t=1
-   `indirect_t1_ucb`: Upper confidence bound for the indirect effect at t=1
-   `indirect_t1_p_value`: P-value for the indirect effect at t=1
-   `indirect_t0_estimate`: Indirect effect estimate at t=0
-   `indirect_t0_std_err`: Standard error of the indirect effect at t=0
-   `indirect_t0_lcb`: Lower confidence bound for the indirect effect at t=0
-   `indirect_t0_ucb`: Upper confidence bound for the indirect effect at t=0
-   `indirect_t0_p_value`: P-value for the indirect effect at t=0
-   `direct_t1_estimate`: Direct effect estimate at t=1
-   `direct_t1_std_err`: Standard error of the direct effect at t=1
-   `direct_t1_lcb`: Lower confidence bound for the direct effect at t=1
-   `direct_t1_ucb`: Upper confidence bound for the direct effect at t=1
-   `direct_t1_p_value`: P-value for the direct effect at t=1
-   `direct_t0_estimate`: Direct effect estimate at t=0
-   `direct_t0_std_err`: Standard error of the direct effect at t=0
-   `direct_t0_lcb`: Lower confidence bound for the direct effect at t=0
-   `direct_t0_ucb`: Upper confidence bound for the direct effect at t=0
-   `direct_t0_p_value`: P-value for the direct effect at t=0
-   `total_effect_estimate`: Total effect estimate
-   `total_effect_std_err`: Standard error of the total effect
-   `total_effect_lcb`: Lower confidence bound for the total effect
-   `total_effect_ucb`: Upper confidence bound for the total effect
-   `total_effect_p_value`: P-value for the total effect

## Development

If you wish to contribute to fastmed, please follow these guidelines:

1.  Fork the repository on GitHub.
2.  Create a new branch for your feature or bugfix.
3.  Commit your changes with clear messages.
4.  Push your branch to your forked repository.
5.  Submit a pull request detailing your changes.

## Testing

fastmed includes a suite of tests to ensure functionality. To run the tests, use the following commands:

``` r
library(devtools)
library(testthat)

# Navigate to the package directory
setwd("path/to/fastmed")

# Run tests
devtools::test()
```

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgements

-   Rcpp
-   RcppParallel
-   data.table

## Contact

For any questions or suggestions, please open an issue on the GitHub repository.

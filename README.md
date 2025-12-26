# fastmed

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg) ![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/fastmed) [![R-CMD-check](https://github.com/daehwankim12/fastmed/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/daehwankim12/fastmed/actions/workflows/R-CMD-check.yaml) 

## Overview

**fastmed** is an R package designed to perform efficient and scalable mediation analysis using multiple linear regressions and bootstrap resampling. Leveraging the power of C++ through Rcpp and parallel processing with RcppParallel, `fastmed` is optimized to handle large datasets and provide real-time results output to CSV files.

## Features

- **Efficient Mediation Analysis:** Conduct mediation analysis with multiple exposure, mediator, and outcome variables.
- **Parallel Processing:** Utilize multiple CPU cores to accelerate computations.
- **Bootstrap Resampling:** Estimate direct, indirect, and total effects with confidence intervals and p-values.
- **Customizable Column Prefixes:** Easily specify prefixes to identify exposure, mediator, and outcome variables.
- **Scalable for Large Datasets:** Handle large-scale data efficiently by processing in chunks.
- **Real-time CSV Output:** Save results directly to CSV files during analysis.

## Installation

You can install the development version of **fastmed** from GitHub using the `remotes` package:

```r
# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install fastmed from GitHub
remotes::install_github("daehwankim12/fastmed")
```

## Usage

Here's a detailed example of how to use the `mediation_analysis` function in fastmed:

```r
library(fastmed)
library(data.table)

# Generate example data
set.seed(123)
n <- 1000
my_data <- data.table(
  Exposure_A1 = rnorm(n),
  Exposure_A2 = rnorm(n),
  Mediator_X1 = rnorm(n),
  Mediator_X2 = rnorm(n),
  Outcome_Y1 = rnorm(n),
  Outcome_Y2 = rnorm(n)
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
  num_threads = 4,       # Number of threads for parallel processing
  seed = 42              # Optional: reproducible results
)

# View results
results <- fread(output_csv)
print(results)
```

### Parameters

- `data`: A data.table or data.frame containing the dataset.
- `columns`: A list with three named elements: exposure, mediator, and outcome. Each should be a character string specifying the prefix of the respective columns.
- `nrep`: (Optional) Number of bootstrap replicates. Default is 1000.
- `output_file`: Path to the output CSV file where results will be saved.
- `num_threads`: (Optional) Number of threads for parallel processing. Defaults to the number of available cores.
- `pert`: (Optional) Uncertainty method, one of `"asymptotic"` or `"bootstrap"`. Default is `"asymptotic"`.
- `seed`: (Optional) Non-negative integer seed for reproducible results across thread counts (within the same build/runtime environment).
- `chunk_size`: (Optional) Maximum number of combinations to buffer per chunk inside the C++ backend; smaller values reduce peak memory usage for very large analyses.
- `output.format`: (Optional) Output CSV schema. `"mediate"` (default) writes mediate-style effect columns; `"legacy"` writes the previous schema.

### Output

The output CSV file will contain detailed results for each combination of exposure, mediator, and outcome variables, including mean estimates, 95% percentile confidence intervals, and p-values for indirect, direct, and total effects.

Default schema (`output.format = "mediate"`):

```
Combination,d0_estimate,d0_ci_lower,d0_ci_upper,d0_p,d1_estimate,d1_ci_lower,d1_ci_upper,d1_p,z0_estimate,z0_ci_lower,z0_ci_upper,z0_p,z1_estimate,z1_ci_lower,z1_ci_upper,z1_p,tau_estimate,tau_ci_lower,tau_ci_upper,tau_p
```

### p-value definition

`fastmed` reports a two-sided, sign-based p-value for the null hypothesis H₀: effect = 0 from the replicate effect samples:

- Let `pos = #{v > 0}`, `neg = #{v < 0}`, `zero = n - pos - neg`
- Split zeros evenly: `pos_eff = pos + 0.5 * zero`
- Continuity correction: `prop = (pos_eff + 1) / (n + 2)`
- Two-sided p-value: `p = 2 * min(prop, 1 - prop)`

### Reproducibility

For schedule-independent reproducibility across different `num_threads` values, pass an explicit `seed`. Output row order is deterministic and follows the generated `(exposure, mediator, outcome)` combination order.

## Performance Considerations

- The package is optimized for parallel processing. Increase `num_threads` to utilize more CPU cores.
- For very large analyses, use `chunk_size` to limit peak memory usage while writing results incrementally.
- Monitor memory usage, especially when increasing `nrep` for bootstrap resampling.

## Troubleshooting

If you encounter issues:

1. Ensure you have the latest version of fastmed installed.
2. Check that all dependencies are up to date.
3. For performance issues, try adjusting `num_threads` or reducing `nrep`.
4. If you encounter a bug, please [open an issue](https://github.com/daehwankim12/fastmed/issues) with a reproducible example.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

fastmed builds upon several powerful R packages:

- [Rcpp](https://www.rcpp.org/) for C++ integration
- [RcppParallel](https://rcppcore.github.io/RcppParallel/) for parallel processing
- [data.table](https://rdatatable.gitlab.io/data.table/) for efficient data manipulation

## Citation

If you use fastmed in your research, please cite it as follows:

```
Kim, D. (2024). fastmed: Fast Mediation Analysis in R. R package version 0.1.0.
https://github.com/daehwankim12/fastmed
```

## Contact

For questions, suggestions, or collaborations, please [open an issue](https://github.com/daehwankim12/fastmed/issues) on the GitHub repository or contact the package maintainer at [kdh5358@snu.ac.kr](mailto:kdh5358@snu.ac.kr).

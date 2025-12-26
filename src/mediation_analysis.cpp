// mediation_analysis.cpp
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppParallel.h>

#include "glm_fit.h"
#include "mediation_worker.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppEigen, RcppParallel)]]

using Eigen::MatrixXd;
using Rcpp::as;
using Rcpp::CharacterVector;
using Rcpp::IntegerVector;
using Rcpp::NumericMatrix;

// [[Rcpp::export]]
void mediation_analysis_cpp(NumericMatrix data,
                            CharacterVector column_names,
                            IntegerVector exposure_col_idx,
                            IntegerVector mediator_col_idx,
                            IntegerVector outcome_col_idx,
                            int nrep,
                            std::string output_file,
                            std::string pert = "asymptotic",
                            uint64_t base_seed = 0,
                            std::string mediator_family = "auto",
                            std::string outcome_family = "auto",
                            bool replace_outcome = false,
                            std::string output_format = "mediate",
                            int chunk_size = 1024,
                            int grain_size = 1) {
    try {
        if (data.nrow() == 0 || data.ncol() == 0) {
            throw std::invalid_argument("Data matrix is empty");
        }
        if (column_names.size() != data.ncol()) {
            throw std::invalid_argument(
                "Column names size does not match data columns");
        }

        const int n = data.nrow();
        const int p_med = 2;
        const int p_out = 3;
        if (n <= p_med) {
            throw std::runtime_error(
                "Insufficient observations for mediator model: n must be > p_med");
        }
        if (n <= p_out) {
            throw std::runtime_error(
                "Insufficient observations for outcome model: n must be > p_out");
        }
        if (nrep <= 0) {
            throw std::invalid_argument(
                "Number of bootstrap replicates must be positive");
        }
        if (pert != "asymptotic" && pert != "bootstrap") {
            throw std::invalid_argument(
                "Invalid perturbation method. Use 'asymptotic' or 'bootstrap'");
        }
        if (chunk_size <= 0) {
            throw std::invalid_argument("chunk_size must be positive");
        }
        if (grain_size <= 0) {
            throw std::invalid_argument("grain_size must be positive");
        }

        std::vector<std::string> column_names_cpp =
            as<std::vector<std::string>>(column_names);

        std::vector<int> exposure_col_idx_cpp = as<std::vector<int>>(exposure_col_idx);
        std::vector<int> mediator_col_idx_cpp = as<std::vector<int>>(mediator_col_idx);
        std::vector<int> outcome_col_idx_cpp = as<std::vector<int>>(outcome_col_idx);

        if (exposure_col_idx_cpp.empty() || mediator_col_idx_cpp.empty() ||
            outcome_col_idx_cpp.empty()) {
            throw std::invalid_argument(
                "Exposure/mediator/outcome index lists must be non-empty");
        }
        for (int v : exposure_col_idx_cpp) {
            if (v < 0 || v >= data.ncol()) {
                throw std::invalid_argument("exposure_col_idx out of bounds");
            }
        }
        for (int v : mediator_col_idx_cpp) {
            if (v < 0 || v >= data.ncol()) {
                throw std::invalid_argument("mediator_col_idx out of bounds");
            }
        }
        for (int v : outcome_col_idx_cpp) {
            if (v < 0 || v >= data.ncol()) {
                throw std::invalid_argument("outcome_col_idx out of bounds");
            }
        }

        Eigen::Map<const MatrixXd> data_map(data.begin(), data.nrow(), data.ncol());

        const size_t e = exposure_col_idx_cpp.size();
        const size_t m = mediator_col_idx_cpp.size();
        const size_t o = outcome_col_idx_cpp.size();

        std::vector<GlmFamily> mediator_fams(m);
        std::vector<GlmFamily> outcome_fams(o);
        std::vector<char> mediator_fams_ok(m, 1);
        std::vector<char> outcome_fams_ok(o, 1);

        if (mediator_family == "auto") {
            for (size_t j = 0; j < m; ++j) {
                try {
                    mediator_fams[j] =
                        parse_family_or_auto(mediator_family,
                                             data_map.col(mediator_col_idx_cpp[j]));
                } catch (...) {
                    mediator_fams_ok[j] = 0;
                    mediator_fams[j] = GlmFamily::Gaussian;
                }
            }
        } else {
            const GlmFamily fam_m =
                parse_family_or_auto(mediator_family,
                                     data_map.col(mediator_col_idx_cpp[0]));
            for (size_t j = 0; j < m; ++j) {
                mediator_fams[j] = fam_m;
            }
        }

        if (outcome_family == "auto") {
            for (size_t j = 0; j < o; ++j) {
                try {
                    outcome_fams[j] =
                        parse_family_or_auto(outcome_family,
                                             data_map.col(outcome_col_idx_cpp[j]));
                } catch (...) {
                    outcome_fams_ok[j] = 0;
                    outcome_fams[j] = GlmFamily::Gaussian;
                }
            }
        } else {
            const GlmFamily fam_y =
                parse_family_or_auto(outcome_family,
                                     data_map.col(outcome_col_idx_cpp[0]));
            for (size_t j = 0; j < o; ++j) {
                outcome_fams[j] = fam_y;
            }
        }

        auto checked_mul = [](size_t a, size_t b) -> size_t {
            if (a == 0 || b == 0) {
                return 0;
            }
            if (a > std::numeric_limits<size_t>::max() / b) {
                throw std::overflow_error("Too many combinations");
            }
            return a * b;
        };

        const size_t em = checked_mul(e, m);
        const size_t num_combinations = checked_mul(em, o);

        std::ofstream output_stream;
        output_stream.open(output_file, std::ios::out | std::ios::trunc);
        if (!output_stream.is_open()) {
            throw std::runtime_error("Failed to open output file: " + output_file);
        }

        const bool legacy_output_schema = [&]() {
            if (output_format == "legacy") {
                return true;
            }
            if (output_format == "mediate") {
                return false;
            }
            throw std::invalid_argument("output_format must be one of: mediate, legacy");
        }();

        std::string header;
        if (legacy_output_schema) {
            header =
                "Combination,ACME_Mean,ACME_2.5%,ACME_97.5%,ACME_p-value,"
                "ADE_Mean,ADE_2.5%,ADE_97.5%,ADE_p-value,"
                "Total_Effect_Mean,Total_Effect_2.5%,Total_Effect_97.5%,Total_"
                "Effect_p-value\n";
        } else {
            header =
                "Combination,"
                "d0_estimate,d0_ci_lower,d0_ci_upper,d0_p,"
                "d1_estimate,d1_ci_lower,d1_ci_upper,d1_p,"
                "z0_estimate,z0_ci_lower,z0_ci_upper,z0_p,"
                "z1_estimate,z1_ci_lower,z1_ci_upper,z1_p,"
                "tau_estimate,tau_ci_lower,tau_ci_upper,tau_p\n";
        }
        output_stream << header;

        const size_t chunk_size_cpp = static_cast<size_t>(chunk_size);
        const size_t grain_size_cpp = static_cast<size_t>(grain_size);

        for (size_t chunk_begin = 0; chunk_begin < num_combinations;
             chunk_begin += chunk_size_cpp) {
            size_t chunk_end = std::min(chunk_begin + chunk_size_cpp, num_combinations);
            std::vector<std::string> chunk_results(chunk_end - chunk_begin);

            MediationWorker worker(data_map,
                                   column_names_cpp,
                                   nrep,
                                   exposure_col_idx_cpp,
                                   mediator_col_idx_cpp,
                                   outcome_col_idx_cpp,
                                   mediator_fams,
                                   outcome_fams,
                                   mediator_fams_ok,
                                   outcome_fams_ok,
                                   pert,
                                   replace_outcome,
                                   legacy_output_schema,
                                   base_seed,
                                   chunk_begin,
                                   chunk_results);

            RcppParallel::parallelFor(chunk_begin, chunk_end, worker, grain_size_cpp);

            for (const auto& line : chunk_results) {
                output_stream << line;
            }
        }
        output_stream.flush();

    } catch (const std::exception& e) {
        Rcpp::stop("Error in mediation_analysis_cpp: %s", e.what());
    }
}


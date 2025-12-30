// mediation_analysis.cpp
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppParallel.h>

#include <R_ext/Random.h>

#include "glm_fit.h"
#include "match_mediation_rng.h"
#include "mediation_worker.h"

#include <algorithm>
#include <cmath>
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

namespace {
std::size_t estimate_match_mediation_rng_bytes(bool bootstrap,
                                               bool ok,
                                               GlmFamily fam_m,
                                               int n,
                                               int nrep) {
    if (!ok || n <= 0 || nrep <= 0) {
        return 0;
    }

    const std::size_t nn = static_cast<std::size_t>(n);
    const std::size_t rr = static_cast<std::size_t>(nrep);
    std::size_t bytes = 0;

    if (!bootstrap) {
        bytes += sizeof(double) * rr * 5;  // coef normals (2 + 3)
        if (fam_m == GlmFamily::Gaussian) {
            bytes += sizeof(double) * rr * nn;  // mediator noise normals
        } else {
            bytes += sizeof(double) * rr * nn * 2;  // mediator uniforms (M1 + M0)
        }
        return bytes;
    }

    bytes += sizeof(int) * rr * nn;  // bootstrap indices
    if (fam_m == GlmFamily::Gaussian) {
        bytes += sizeof(double) * (rr + 1) * nn;  // t0 + reps
    } else {
        bytes += sizeof(double) * (rr + 1) * nn * 2;  // t0 + reps (M1 + M0)
    }
    return bytes;
}

inline void decode_combination_indices(std::size_t idx,
                                      std::size_t e,
                                      std::size_t m,
                                      std::size_t& exp_list_idx,
                                      std::size_t& med_list_idx,
                                      std::size_t& out_list_idx) {
    exp_list_idx = idx % e;
    med_list_idx = (idx / e) % m;
    out_list_idx = idx / (e * m);
}

MatchMediationRngBlock build_match_mediation_rng_block_asymptotic(
    std::size_t block_begin,
    std::size_t block_end,
    std::size_t e,
    std::size_t m,
    const std::vector<GlmFamily>& mediator_fams,
    const std::vector<char>& mediator_fams_ok,
    const std::vector<char>& outcome_fams_ok,
    int n,
    int nrep) {
    MatchMediationRngBlock block;
    block.begin = block_begin;
    block.end = block_end;
    block.n = n;
    block.nrep = nrep;
    block.bootstrap = false;
    block.slices.resize(block_end - block_begin);

    std::size_t ok_count = 0;
    std::size_t gaussian_count = 0;
    std::size_t other_count = 0;
    for (std::size_t idx = block_begin; idx < block_end; ++idx) {
        std::size_t exp_list_idx = 0;
        std::size_t med_list_idx = 0;
        std::size_t out_list_idx = 0;
        decode_combination_indices(idx, e, m, exp_list_idx, med_list_idx, out_list_idx);
        const bool ok =
            (mediator_fams_ok[med_list_idx] != 0) && (outcome_fams_ok[out_list_idx] != 0);
        if (!ok) {
            continue;
        }
        ++ok_count;
        if (mediator_fams[med_list_idx] == GlmFamily::Gaussian) {
            ++gaussian_count;
        } else {
            ++other_count;
        }
    }

    const std::size_t nn = static_cast<std::size_t>(n);
    const std::size_t rr = static_cast<std::size_t>(nrep);
    block.coef_normals.reserve(ok_count * rr * 5);
    block.noise_normals.reserve(gaussian_count * rr * nn);
    block.noise_uniforms.reserve(other_count * rr * nn * 2);

    for (std::size_t idx = block_begin; idx < block_end; ++idx) {
        const std::size_t local = idx - block_begin;
        MatchMediationRngSlice& slice = block.slices[local];

        std::size_t exp_list_idx = 0;
        std::size_t med_list_idx = 0;
        std::size_t out_list_idx = 0;
        decode_combination_indices(idx, e, m, exp_list_idx, med_list_idx, out_list_idx);
        const bool ok =
            (mediator_fams_ok[med_list_idx] != 0) && (outcome_fams_ok[out_list_idx] != 0);
        if (!ok) {
            continue;
        }

        slice.coef_offset = block.coef_normals.size();
        // mediator model coef draws (nrep x 2), row-major RNG consumption
        for (int rep_idx = 0; rep_idx < nrep; ++rep_idx) {
            for (int j = 0; j < 2; ++j) {
                block.coef_normals.push_back(norm_rand());
            }
        }
        // outcome model coef draws (nrep x 3), row-major RNG consumption
        for (int rep_idx = 0; rep_idx < nrep; ++rep_idx) {
            for (int j = 0; j < 3; ++j) {
                block.coef_normals.push_back(norm_rand());
            }
        }

        const GlmFamily fam_m = mediator_fams[med_list_idx];
        if (fam_m == GlmFamily::Gaussian) {
            slice.noise_kind = MatchMediationNoiseKind::Normals;
            slice.noise_offset = block.noise_normals.size();
            // E matrix (nrep x n) consumed column-major: j outer, rep inner
            for (int j = 0; j < n; ++j) {
                for (int rep_idx = 0; rep_idx < nrep; ++rep_idx) {
                    block.noise_normals.push_back(norm_rand());
                }
            }
        } else {
            slice.noise_kind = MatchMediationNoiseKind::Uniforms;
            slice.noise_offset = block.noise_uniforms.size();
            // M1 uniforms (nrep x n), column-major: j outer, rep inner
            for (int j = 0; j < n; ++j) {
                for (int rep_idx = 0; rep_idx < nrep; ++rep_idx) {
                    block.noise_uniforms.push_back(unif_rand());
                }
            }
            // M0 uniforms (nrep x n), column-major: j outer, rep inner
            for (int j = 0; j < n; ++j) {
                for (int rep_idx = 0; rep_idx < nrep; ++rep_idx) {
                    block.noise_uniforms.push_back(unif_rand());
                }
            }
        }
    }

    return block;
}

MatchMediationRngBlock build_match_mediation_rng_block_bootstrap(
    std::size_t block_begin,
    std::size_t block_end,
    std::size_t e,
    std::size_t m,
    const std::vector<GlmFamily>& mediator_fams,
    const std::vector<char>& mediator_fams_ok,
    const std::vector<char>& outcome_fams_ok,
    int n,
    int nrep) {
    MatchMediationRngBlock block;
    block.begin = block_begin;
    block.end = block_end;
    block.n = n;
    block.nrep = nrep;
    block.bootstrap = true;
    block.slices.resize(block_end - block_begin);

    std::size_t ok_count = 0;
    std::size_t gaussian_count = 0;
    std::size_t other_count = 0;
    for (std::size_t idx = block_begin; idx < block_end; ++idx) {
        std::size_t exp_list_idx = 0;
        std::size_t med_list_idx = 0;
        std::size_t out_list_idx = 0;
        decode_combination_indices(idx, e, m, exp_list_idx, med_list_idx, out_list_idx);
        const bool ok =
            (mediator_fams_ok[med_list_idx] != 0) && (outcome_fams_ok[out_list_idx] != 0);
        if (!ok) {
            continue;
        }
        ++ok_count;
        if (mediator_fams[med_list_idx] == GlmFamily::Gaussian) {
            ++gaussian_count;
        } else {
            ++other_count;
        }
    }

    const std::size_t nn = static_cast<std::size_t>(n);
    const std::size_t rr = static_cast<std::size_t>(nrep);
    block.bootstrap_indices.reserve(ok_count * rr * nn);
    block.noise_normals.reserve(gaussian_count * (rr + 1) * nn);
    block.noise_uniforms.reserve(other_count * (rr + 1) * nn * 2);

    for (std::size_t idx = block_begin; idx < block_end; ++idx) {
        const std::size_t local = idx - block_begin;
        MatchMediationRngSlice& slice = block.slices[local];

        std::size_t exp_list_idx = 0;
        std::size_t med_list_idx = 0;
        std::size_t out_list_idx = 0;
        decode_combination_indices(idx, e, m, exp_list_idx, med_list_idx, out_list_idx);
        const bool ok =
            (mediator_fams_ok[med_list_idx] != 0) && (outcome_fams_ok[out_list_idx] != 0);
        if (!ok) {
            continue;
        }

        slice.idx_offset = block.bootstrap_indices.size();
        // indices consumed column-major: j outer, rep inner
        for (int j = 0; j < n; ++j) {
            for (int rep_idx = 0; rep_idx < nrep; ++rep_idx) {
                block.bootstrap_indices.push_back(
                    static_cast<int>(R_unif_index(static_cast<double>(n))));
            }
        }

        const GlmFamily fam_m = mediator_fams[med_list_idx];
        if (fam_m == GlmFamily::Gaussian) {
            slice.noise_kind = MatchMediationNoiseKind::Normals;
            slice.noise_offset = block.noise_normals.size();
            // t0 noise: i outer
            for (int i = 0; i < n; ++i) {
                block.noise_normals.push_back(norm_rand());
            }
            // rep noise: rep outer, i inner
            for (int rep_idx = 0; rep_idx < nrep; ++rep_idx) {
                for (int i = 0; i < n; ++i) {
                    block.noise_normals.push_back(norm_rand());
                }
            }
        } else {
            slice.noise_kind = MatchMediationNoiseKind::Uniforms;
            slice.noise_offset = block.noise_uniforms.size();
            // t0 noise: i outer, (M1, M0) per i
            for (int i = 0; i < n; ++i) {
                block.noise_uniforms.push_back(unif_rand());  // M1
                block.noise_uniforms.push_back(unif_rand());  // M0
            }
            // rep noise: rep outer, i inner, (M1, M0) per i
            for (int rep_idx = 0; rep_idx < nrep; ++rep_idx) {
                for (int i = 0; i < n; ++i) {
                    block.noise_uniforms.push_back(unif_rand());  // M1
                    block.noise_uniforms.push_back(unif_rand());  // M0
                }
            }
        }
    }

    return block;
}
}  // namespace

// [[Rcpp::export]]
void mediation_analysis_cpp(NumericMatrix data,
                            CharacterVector column_names,
                            IntegerVector exposure_col_idx,
                            IntegerVector mediator_col_idx,
                            IntegerVector outcome_col_idx,
                            int nrep,
                            std::string output_file,
                            Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue,
                            std::string pert = "asymptotic",
                            uint64_t base_seed = 0,
                            std::string mediator_family = "auto",
                            std::string outcome_family = "auto",
                            bool replace_outcome = false,
                            std::string output_format = "mediate",
                            double treat_value = 1.0,
                            double control_value = 0.0,
                            int chunk_size = 1024,
                            int grain_size = 1,
                            bool overwrite = true,
                            bool excel_safe_csv = false,
                            bool match_mediation = false) {
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
        const int max_nrep = std::numeric_limits<int>::max() / 3;
        if (nrep > max_nrep) {
            throw std::invalid_argument("nrep is too large (max " +
                                        std::to_string(max_nrep) + ")");
        }
        if (pert != "asymptotic" && pert != "bootstrap") {
            throw std::invalid_argument(
                "Invalid perturbation method. Use 'asymptotic' or 'bootstrap'");
        }
        if (!std::isfinite(treat_value) || !std::isfinite(control_value)) {
            throw std::invalid_argument("treat_value and control_value must be finite");
        }
        if (treat_value == control_value) {
            throw std::invalid_argument(
                "treat_value and control_value must be different");
        }
        if (chunk_size <= 0) {
            throw std::invalid_argument("chunk_size must be positive");
        }
        if (grain_size <= 0) {
            throw std::invalid_argument("grain_size must be positive");
        }
        if (match_mediation && base_seed > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
            throw std::invalid_argument("match_mediation requires base_seed <= .Machine$integer.max");
        }

        std::vector<std::string> column_names_cpp =
            as<std::vector<std::string>>(column_names);

        Eigen::VectorXd weights_cpp = Eigen::VectorXd::Ones(n);
        if (weights.isNotNull()) {
            Rcpp::NumericVector w(weights.get());
            if (w.size() != n) {
                throw std::invalid_argument("weights must have length nrow(data)");
            }
            double wsum = 0.0;
            for (int i = 0; i < n; ++i) {
                const double wi = w[i];
                if (!std::isfinite(wi)) {
                    throw std::invalid_argument("weights must be finite");
                }
                if (wi < 0.0) {
                    throw std::invalid_argument("weights must be non-negative");
                }
                weights_cpp[i] = wi;
                wsum += wi;
            }
            if (!(wsum > 0.0) || !std::isfinite(wsum)) {
                throw std::invalid_argument("weights must sum to a positive finite value");
            }
        }

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

        if (!overwrite) {
            std::ifstream existing(output_file);
            if (existing.good()) {
                throw std::runtime_error("Output file exists and overwrite is FALSE: " +
                                         output_file);
            }
        }

        std::ofstream output_stream;
        output_stream.open(output_file, std::ios::out | std::ios::trunc);
        if (!output_stream.is_open()) {
            throw std::runtime_error("Failed to open output file: " + output_file);
        }

        // When match_mediation is enabled we rely on R's RNG stream. Any seeding
        // should be handled on the R side (e.g., via set.seed()) prior to
        // calling into this routine.

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

            if (match_mediation) {
                const bool bootstrap = (pert == "bootstrap");
                const std::size_t max_rng_bytes = 256ULL * 1024ULL * 1024ULL;

                const std::size_t e_size = exposure_col_idx_cpp.size();
                const std::size_t m_size = mediator_col_idx_cpp.size();

                MediationWorker serial_worker(data_map,
                                              weights_cpp,
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
                                              excel_safe_csv,
                                              treat_value,
                                              control_value,
                                              base_seed,
                                              match_mediation,
                                              chunk_begin,
                                              chunk_results,
                                              nullptr);

                Eigen::VectorXd off_m = Eigen::VectorXd::Zero(n);
                Eigen::VectorXd off_y = Eigen::VectorXd::Zero(n);
                Eigen::MatrixXd X_med(n, 2);
                Eigen::MatrixXd X_out(n, 3);
                X_med.col(0).setOnes();
                X_out.col(0).setOnes();

                for (std::size_t cursor = chunk_begin; cursor < chunk_end;) {
                    std::size_t exp_list_idx = 0;
                    std::size_t med_list_idx = 0;
                    std::size_t out_list_idx = 0;
                    decode_combination_indices(cursor,
                                              e_size,
                                              m_size,
                                              exp_list_idx,
                                              med_list_idx,
                                              out_list_idx);

                    const bool ok = (mediator_fams_ok[med_list_idx] != 0) &&
                                    (outcome_fams_ok[out_list_idx] != 0);
                    const GlmFamily fam_m = mediator_fams[med_list_idx];

                    if (ok && fam_m == GlmFamily::Poisson) {
                        chunk_results[cursor - chunk_begin] =
                            serial_worker.process_combination_serial(
                                cursor, X_med, X_out, off_m, off_y);
                        ++cursor;
                        Rcpp::checkUserInterrupt();
                        continue;
                    }

                    const std::size_t block_begin = cursor;
                    std::size_t block_end = cursor;
                    std::size_t rng_bytes = 0;

                    while (block_end < chunk_end) {
                        decode_combination_indices(block_end,
                                                  e_size,
                                                  m_size,
                                                  exp_list_idx,
                                                  med_list_idx,
                                                  out_list_idx);

                        const bool ok_block = (mediator_fams_ok[med_list_idx] != 0) &&
                                              (outcome_fams_ok[out_list_idx] != 0);
                        const GlmFamily fam_m_block = mediator_fams[med_list_idx];
                        if (ok_block && fam_m_block == GlmFamily::Poisson) {
                            break;
                        }

                        const std::size_t add = estimate_match_mediation_rng_bytes(
                            bootstrap, ok_block, fam_m_block, n, nrep);

                        if (block_end > block_begin && rng_bytes + add > max_rng_bytes) {
                            break;
                        }
                        rng_bytes += add;
                        ++block_end;
                    }

                    MatchMediationRngBlock rng_block;
                    {
                        Rcpp::RNGScope rng_scope;
                        if (bootstrap) {
                            rng_block = build_match_mediation_rng_block_bootstrap(
                                block_begin,
                                block_end,
                                e_size,
                                m_size,
                                mediator_fams,
                                mediator_fams_ok,
                                outcome_fams_ok,
                                n,
                                nrep);
                        } else {
                            rng_block = build_match_mediation_rng_block_asymptotic(
                                block_begin,
                                block_end,
                                e_size,
                                m_size,
                                mediator_fams,
                                mediator_fams_ok,
                                outcome_fams_ok,
                                n,
                                nrep);
                        }
                    }

                    MediationWorker worker(data_map,
                                           weights_cpp,
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
                                           excel_safe_csv,
                                           treat_value,
                                           control_value,
                                           base_seed,
                                           match_mediation,
                                           chunk_begin,
                                           chunk_results,
                                           &rng_block);

                    RcppParallel::parallelFor(block_begin, block_end, worker, grain_size_cpp);
                    cursor = block_end;
                    Rcpp::checkUserInterrupt();
                }
            } else {
                MediationWorker worker(data_map,
                                       weights_cpp,
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
                                       excel_safe_csv,
                                       treat_value,
                                       control_value,
                                       base_seed,
                                       match_mediation,
                                       chunk_begin,
                                       chunk_results,
                                       nullptr);
                RcppParallel::parallelFor(chunk_begin, chunk_end, worker, grain_size_cpp);
            }

            for (const auto& line : chunk_results) {
                output_stream << line;
            }

            Rcpp::checkUserInterrupt();
        }
        output_stream.flush();

    } catch (const std::exception& e) {
        Rcpp::stop("Error in mediation_analysis_cpp: %s", e.what());
    }
}

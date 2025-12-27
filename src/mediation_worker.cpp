#include "mediation_worker.h"

#include "fastmed_csv.h"
#include "fastmed_glm_link.h"
#include "fastmed_rng.h"
#include "fastmed_stats.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {
MatrixXd cholesky_lower_or_throw(MatrixXd cov, const std::string& context) {
    cov = 0.5 * (cov + cov.transpose());
    Eigen::LLT<MatrixXd> chol(cov);
    if (chol.info() != Eigen::Success) {
        throw std::runtime_error("Cholesky decomposition failed for " + context);
    }
    MatrixXd L = chol.matrixL();
    if (!L.allFinite() || (L.diagonal().array() <= 0.0).any()) {
        throw std::runtime_error("Cholesky decomposition failed for " + context);
    }
    return L;
}
}  // namespace

MediationWorker::MediationWorker(const Eigen::Map<const MatrixXd>& data_,
                                 const VectorXd& weights_,
                                 const std::vector<std::string>& column_names_,
                                 int nrep_,
                                 const std::vector<int>& exposure_col_idx_,
                                 const std::vector<int>& mediator_col_idx_,
                                 const std::vector<int>& outcome_col_idx_,
                                 const std::vector<GlmFamily>& mediator_fams_,
                                 const std::vector<GlmFamily>& outcome_fams_,
                                 const std::vector<char>& mediator_fams_ok_,
                                 const std::vector<char>& outcome_fams_ok_,
                                 const std::string& pert_method_,
                                 bool replace_outcome_,
                                 bool legacy_output_schema_,
                                 double treat_value_,
                                 double control_value_,
                                 uint64_t base_seed_,
                                 size_t chunk_begin_,
                                 std::vector<std::string>& output_lines_)
    : data(data_),
      weights(weights_),
      column_names(column_names_),
      nrep(nrep_),
      exposure_col_idx(exposure_col_idx_),
      mediator_col_idx(mediator_col_idx_),
      outcome_col_idx(outcome_col_idx_),
      mediator_fams(mediator_fams_),
      outcome_fams(outcome_fams_),
      mediator_fams_ok(mediator_fams_ok_),
      outcome_fams_ok(outcome_fams_ok_),
      pert_method(pert_method_),
      replace_outcome(replace_outcome_),
      legacy_output_schema(legacy_output_schema_),
      treat_value(treat_value_),
      control_value(control_value_),
      base_seed(base_seed_),
      chunk_begin(chunk_begin_),
      output_lines(output_lines_) {}

void MediationWorker::operator()(std::size_t begin, std::size_t end) {
    const int n = static_cast<int>(data.rows());
    VectorXd off_m = VectorXd::Zero(n);
    VectorXd off_y = VectorXd::Zero(n);

    MatrixXd X_med(n, 2);
    X_med.col(0).setOnes();

    MatrixXd X_out(n, 3);
    X_out.col(0).setOnes();

    for (std::size_t idx = begin; idx < end; ++idx) {
        output_lines[idx - chunk_begin] =
            process_combination(idx, X_med, X_out, off_m, off_y);
    }
}

std::string MediationWorker::format_na_row(const std::string& exposure_col,
                                          const std::string& mediator_col,
                                          const std::string& outcome_col) {
    std::string combination = exposure_col + "_" + mediator_col + "_" + outcome_col;

    std::string result;
    result.reserve(combination.size() + 128);
    result += csv_escape(combination);
    const int na_cols = legacy_output_schema ? 12 : 20;
    for (int i = 0; i < na_cols; ++i) {
        result += ",NA";
    }
    result.push_back('\n');
    return result;
}

std::string MediationWorker::process_combination(std::size_t idx,
                                                 MatrixXd& X_med,
                                                 MatrixXd& X_out,
                                                 VectorXd& off_m,
                                                 VectorXd& off_y) {
    const size_t e = exposure_col_idx.size();
    const size_t m = mediator_col_idx.size();

    const size_t exp_list_idx = idx % e;
    const size_t med_list_idx = (idx / e) % m;
    const size_t out_list_idx = idx / (e * m);

    const int exp_idx = exposure_col_idx[exp_list_idx];
    const int med_idx = mediator_col_idx[med_list_idx];
    const int out_idx = outcome_col_idx[out_list_idx];

    const std::string& exposure_col = column_names[exp_idx];
    const std::string& mediator_col = column_names[med_idx];
    const std::string& outcome_col = column_names[out_idx];

    const auto exposure = data.col(exp_idx);
    const auto mediator = data.col(med_idx);
    const auto outcome = data.col(out_idx);

    const uint64_t global_combination_idx = static_cast<uint64_t>(idx);

    if (!mediator_fams_ok[med_list_idx] || !outcome_fams_ok[out_list_idx]) {
        return format_na_row(exposure_col, mediator_col, outcome_col);
    }
    const GlmFamily fam_m = mediator_fams[med_list_idx];
    const GlmFamily fam_y = outcome_fams[out_list_idx];

    try {
        const int n = static_cast<int>(data.rows());

        X_med.col(1) = exposure;

        X_out.col(1) = mediator;
        X_out.col(2) = exposure;

        const int glm_maxit = 25;
        const double glm_epsilon = 1e-8;
        const double glm_qr_tol = 1e-12;

        GlmFit fit_m =
            glm_fit_irls_qr(X_med, mediator, fam_m, weights, off_m, glm_maxit,
                            glm_epsilon, glm_qr_tol);

        GlmFit fit_y =
            glm_fit_irls_qr(X_out, outcome, fam_y, weights, off_y, glm_maxit,
                            glm_epsilon, glm_qr_tol);

        std::vector<BootstrapResult> bootstrap_results;
        BootstrapResult t0_result{};
        const BootstrapResult* t0_ptr = nullptr;
        if (pert_method == "asymptotic") {
            bootstrap_results = perform_bootstrap_asymptotic(
                fit_m, fit_y, fam_m, fam_y, n, global_combination_idx, exposure,
                outcome, replace_outcome);
        } else if (pert_method == "bootstrap") {
            std::mt19937_64 rng_t0(
                derive_seed(base_seed,
                            global_combination_idx,
                            std::numeric_limits<uint64_t>::max()));
            const double sigma2_m =
                (fam_m == GlmFamily::Gaussian) ? fit_m.dispersion : 1.0;
            t0_result = simulate_effect_draw(fit_m.coef,
                                             fit_y.coef,
                                             fam_m,
                                             fam_y,
                                             sigma2_m,
                                             rng_t0,
                                             exposure,
                                             outcome,
                                             weights,
                                             replace_outcome);
            t0_ptr = &t0_result;
            bootstrap_results = perform_bootstrap_resample(
                exposure, mediator, outcome, fam_m, fam_y, n, global_combination_idx,
                replace_outcome);
        } else {
            throw std::invalid_argument("Unknown perturbation method: " + pert_method);
        }

        return format_results(exposure_col,
                              mediator_col,
                              outcome_col,
                              bootstrap_results,
                              t0_ptr);
    } catch (...) {
        return format_na_row(exposure_col, mediator_col, outcome_col);
    }
}

std::vector<BootstrapResult> MediationWorker::perform_bootstrap_asymptotic(
    const GlmFit& fit_m,
    const GlmFit& fit_y,
    GlmFamily fam_m,
    GlmFamily fam_y,
    int n,
    uint64_t global_combination_idx,
    const Eigen::Ref<const VectorXd>& exposure_obs,
    const Eigen::Ref<const VectorXd>& outcome_obs,
    bool replace_outcome_) {
    std::vector<BootstrapResult> results;
    results.reserve(nrep);

    MatrixXd Lm = cholesky_lower_or_throw(fit_m.vcov, "mediator model vcov");
    MatrixXd Ly = cholesky_lower_or_throw(fit_y.vcov, "outcome model vcov");

    std::normal_distribution<double> stdnorm(0.0, 1.0);
    const double sigma2_m =
        (fam_m == GlmFamily::Gaussian) ? fit_m.dispersion : 1.0;

    for (int rep_idx = 0; rep_idx < nrep; ++rep_idx) {
        std::mt19937_64 rng(derive_seed(base_seed, global_combination_idx, rep_idx));

        VectorXd zm(fit_m.coef.size());
        VectorXd zy(fit_y.coef.size());
        for (int i = 0; i < zm.size(); ++i) {
            zm[i] = stdnorm(rng);
        }
        for (int i = 0; i < zy.size(); ++i) {
            zy[i] = stdnorm(rng);
        }

        VectorXd bm = fit_m.coef + Lm * zm;
        VectorXd by = fit_y.coef + Ly * zy;
        results.push_back(simulate_effect_draw(bm,
                                               by,
                                               fam_m,
                                               fam_y,
                                               sigma2_m,
                                               rng,
                                               exposure_obs,
                                               outcome_obs,
                                               weights,
                                               replace_outcome_));
    }

    return results;
}

std::vector<BootstrapResult> MediationWorker::perform_bootstrap_resample(
    const Eigen::Ref<const VectorXd>& exposure,
    const Eigen::Ref<const VectorXd>& mediator,
    const Eigen::Ref<const VectorXd>& outcome,
    GlmFamily fam_m,
    GlmFamily fam_y,
    int n,
    uint64_t global_combination_idx,
    bool replace_outcome_) {
    std::vector<BootstrapResult> results;
    results.reserve(nrep);

    const int MAX_ATTEMPTS = 3 * nrep;
    int rep_idx = 0;
    int total_attempts = 0;
    std::string last_exception_msg;
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

    std::uniform_int_distribution<int> dis(0, n - 1);

    VectorXd off_m = VectorXd::Zero(n);
    VectorXd off_y = VectorXd::Zero(n);

    VectorXd boot_exposure(n);
    VectorXd boot_mediator(n);
    VectorXd boot_outcome(n);
    VectorXd boot_weights(n);

    MatrixXd X_med_boot(n, 2);
    MatrixXd X_out_boot(n, 3);

    while (rep_idx < nrep) {
        if (total_attempts >= MAX_ATTEMPTS) {
            std::ostringstream oss;
            oss << "Bootstrap failed: " << (total_attempts - rep_idx) << " failures, "
                << rep_idx << " successes after " << MAX_ATTEMPTS
                << " attempts. Last error: " << last_exception_msg;
            throw std::runtime_error(oss.str());
        }

        std::mt19937_64 rng(derive_seed(base_seed, global_combination_idx, rep_idx));

        bool success = false;
        while (!success) {
            if (total_attempts >= MAX_ATTEMPTS) {
                std::ostringstream oss;
                oss << "Bootstrap failed: " << (total_attempts - rep_idx)
                    << " failures, " << rep_idx << " successes after " << MAX_ATTEMPTS
                    << " attempts. Last error: " << last_exception_msg;
                throw std::runtime_error(oss.str());
            }
            ++total_attempts;

            for (int j = 0; j < n; ++j) {
                int sample_idx = dis(rng);
                boot_exposure[j] = exposure[sample_idx];
                boot_mediator[j] = mediator[sample_idx];
                boot_outcome[j] = outcome[sample_idx];
                boot_weights[j] = weights[sample_idx];
            }

            X_med_boot.col(0).setOnes();
            X_med_boot.col(1) = boot_exposure;

            X_out_boot.col(0).setOnes();
            X_out_boot.col(1) = boot_mediator;
            X_out_boot.col(2) = boot_exposure;

            const int glm_maxit = 25;
            const double glm_epsilon = 1e-8;
            const double glm_qr_tol = 1e-12;

            try {
                GlmFit fit_m_boot =
                    glm_fit_irls_qr(X_med_boot,
                                    boot_mediator,
                                    fam_m,
                                    boot_weights,
                                    off_m,
                                    glm_maxit, glm_epsilon, glm_qr_tol);
                GlmFit fit_y_boot =
                    glm_fit_irls_qr(X_out_boot,
                                    boot_outcome,
                                    fam_y,
                                    boot_weights,
                                    off_y,
                                    glm_maxit, glm_epsilon, glm_qr_tol);

                const double sigma2_m =
                    (fam_m == GlmFamily::Gaussian) ? fit_m_boot.dispersion : 1.0;

                results.push_back(simulate_effect_draw(fit_m_boot.coef,
                                                       fit_y_boot.coef,
                                                       fam_m,
                                                       fam_y,
                                                       sigma2_m,
                                                       rng,
                                                       exposure,
                                                       outcome,
                                                       weights,
                                                       replace_outcome_));
                success = true;
                ++rep_idx;
            } catch (const std::exception& e) {
                last_exception_msg = e.what();
            }
        }
    }

    return results;
}

BootstrapResult MediationWorker::simulate_effect_draw(
    const VectorXd& bm,
    const VectorXd& by,
    GlmFamily fam_m,
    GlmFamily fam_y,
    double sigma2_m,
    std::mt19937_64& rng,
    const Eigen::Ref<const VectorXd>& exposure_obs,
    const Eigen::Ref<const VectorXd>& outcome_obs,
    const Eigen::Ref<const VectorXd>& weights_obs,
    bool replace_outcome_) {
    const int n = static_cast<int>(exposure_obs.size());
    if (outcome_obs.size() != n) {
        throw std::runtime_error("simulate_effect_draw: observed size mismatch");
    }
    if (weights_obs.size() != n) {
        throw std::runtime_error("simulate_effect_draw: weights size mismatch");
    }

    const double weight_sum = weights_obs.sum();
    if (!(weight_sum > 0.0) || !std::isfinite(weight_sum)) {
        throw std::runtime_error("simulate_effect_draw: non-positive weight sum");
    }

    const double X0 = control_value;
    const double X1 = treat_value;

    double eta_m0 = bm[0] + bm[1] * X0;
    double eta_m1 = bm[0] + bm[1] * X1;

    if (fam_y == GlmFamily::Gaussian && !replace_outcome_) {
        double sum_wM0 = 0.0;
        double sum_wM1 = 0.0;

        if (fam_m == GlmFamily::Gaussian) {
            double sd = std::sqrt(std::max(0.0, sigma2_m));
            std::normal_distribution<double> nd(0.0, sd);
            for (int i = 0; i < n; ++i) {
                double e = nd(rng);
                sum_wM0 += weights_obs[i] * (eta_m0 + e);
                sum_wM1 += weights_obs[i] * (eta_m1 + e);
            }
        } else if (fam_m == GlmFamily::Binomial) {
            double p0 = inv_logit(eta_m0);
            double p1 = inv_logit(eta_m1);
            std::uniform_real_distribution<double> unif(0.0, 1.0);
            for (int i = 0; i < n; ++i) {
                sum_wM0 += weights_obs[i] * ((unif(rng) < p0) ? 1.0 : 0.0);
                sum_wM1 += weights_obs[i] * ((unif(rng) < p1) ? 1.0 : 0.0);
            }
        } else {  // Poisson
            double lam0 = safe_exp(eta_m0);
            double lam1 = safe_exp(eta_m1);
            std::poisson_distribution<int> p0(lam0);
            std::poisson_distribution<int> p1(lam1);
            for (int i = 0; i < n; ++i) {
                sum_wM0 += weights_obs[i] * static_cast<double>(p0(rng));
                sum_wM1 += weights_obs[i] * static_cast<double>(p1(rng));
            }
        }

        const double mean_M0 = sum_wM0 / weight_sum;
        const double mean_M1 = sum_wM1 / weight_sum;

        const double mean_y00 = by[0] + by[1] * mean_M0 + by[2] * X0;
        const double mean_y01 = by[0] + by[1] * mean_M1 + by[2] * X0;
        const double mean_y10 = by[0] + by[1] * mean_M0 + by[2] * X1;
        const double mean_y11 = by[0] + by[1] * mean_M1 + by[2] * X1;

        const double d0 = mean_y01 - mean_y00;   // ACME(control)
        const double d1 = mean_y11 - mean_y10;   // ACME(treated)
        const double z0 = mean_y10 - mean_y00;   // ADE(control)
        const double z1 = mean_y11 - mean_y01;   // ADE(treated)
        const double tau = mean_y11 - mean_y00;  // total

        return {d0, d1, z0, z1, tau};
    }

    double sum_y00 = 0.0;
    double sum_y01 = 0.0;
    double sum_y10 = 0.0;
    double sum_y11 = 0.0;

    if (fam_m == GlmFamily::Gaussian) {
        double sd = std::sqrt(std::max(0.0, sigma2_m));
        std::normal_distribution<double> nd(0.0, sd);
        for (int i = 0; i < n; ++i) {
            double e = nd(rng);
            const double M0_i = eta_m0 + e;
            const double M1_i = eta_m1 + e;

            const double eta00 = by[0] + by[1] * M0_i + by[2] * X0;
            const double eta01 = by[0] + by[1] * M1_i + by[2] * X0;
            const double eta10 = by[0] + by[1] * M0_i + by[2] * X1;
            const double eta11 = by[0] + by[1] * M1_i + by[2] * X1;

            double y00 = linkinv(eta00, fam_y);
            double y01 = linkinv(eta01, fam_y);
            double y10 = linkinv(eta10, fam_y);
            double y11 = linkinv(eta11, fam_y);

            if (replace_outcome_) {
                const double xi = exposure_obs[i];
                if (std::fabs(xi - X0) < 1e-12) {
                    y00 = outcome_obs[i];
                }
                if (std::fabs(xi - X1) < 1e-12) {
                    y11 = outcome_obs[i];
                }
            }

            const double wi = weights_obs[i];
            sum_y00 += wi * y00;
            sum_y01 += wi * y01;
            sum_y10 += wi * y10;
            sum_y11 += wi * y11;
        }
    } else if (fam_m == GlmFamily::Binomial) {
        double p0 = inv_logit(eta_m0);
        double p1 = inv_logit(eta_m1);
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        for (int i = 0; i < n; ++i) {
            const double M0_i = (unif(rng) < p0) ? 1.0 : 0.0;
            const double M1_i = (unif(rng) < p1) ? 1.0 : 0.0;

            const double eta00 = by[0] + by[1] * M0_i + by[2] * X0;
            const double eta01 = by[0] + by[1] * M1_i + by[2] * X0;
            const double eta10 = by[0] + by[1] * M0_i + by[2] * X1;
            const double eta11 = by[0] + by[1] * M1_i + by[2] * X1;

            double y00 = linkinv(eta00, fam_y);
            double y01 = linkinv(eta01, fam_y);
            double y10 = linkinv(eta10, fam_y);
            double y11 = linkinv(eta11, fam_y);

            if (replace_outcome_) {
                const double xi = exposure_obs[i];
                if (std::fabs(xi - X0) < 1e-12) {
                    y00 = outcome_obs[i];
                }
                if (std::fabs(xi - X1) < 1e-12) {
                    y11 = outcome_obs[i];
                }
            }

            const double wi = weights_obs[i];
            sum_y00 += wi * y00;
            sum_y01 += wi * y01;
            sum_y10 += wi * y10;
            sum_y11 += wi * y11;
        }
    } else {  // Poisson
        double lam0 = safe_exp(eta_m0);
        double lam1 = safe_exp(eta_m1);
        std::poisson_distribution<int> p0(lam0);
        std::poisson_distribution<int> p1(lam1);
        for (int i = 0; i < n; ++i) {
            const double M0_i = static_cast<double>(p0(rng));
            const double M1_i = static_cast<double>(p1(rng));

            const double eta00 = by[0] + by[1] * M0_i + by[2] * X0;
            const double eta01 = by[0] + by[1] * M1_i + by[2] * X0;
            const double eta10 = by[0] + by[1] * M0_i + by[2] * X1;
            const double eta11 = by[0] + by[1] * M1_i + by[2] * X1;

            double y00 = linkinv(eta00, fam_y);
            double y01 = linkinv(eta01, fam_y);
            double y10 = linkinv(eta10, fam_y);
            double y11 = linkinv(eta11, fam_y);

            if (replace_outcome_) {
                const double xi = exposure_obs[i];
                if (std::fabs(xi - X0) < 1e-12) {
                    y00 = outcome_obs[i];
                }
                if (std::fabs(xi - X1) < 1e-12) {
                    y11 = outcome_obs[i];
                }
            }

            const double wi = weights_obs[i];
            sum_y00 += wi * y00;
            sum_y01 += wi * y01;
            sum_y10 += wi * y10;
            sum_y11 += wi * y11;
        }
    }

    const double mean_y00 = sum_y00 / weight_sum;
    const double mean_y01 = sum_y01 / weight_sum;
    const double mean_y10 = sum_y10 / weight_sum;
    const double mean_y11 = sum_y11 / weight_sum;

    const double d0 = mean_y01 - mean_y00;   // ACME(control)
    const double d1 = mean_y11 - mean_y10;   // ACME(treated)
    const double z0 = mean_y10 - mean_y00;   // ADE(control)
    const double z1 = mean_y11 - mean_y01;   // ADE(treated)
    const double tau = mean_y11 - mean_y00;  // total

    return {d0, d1, z0, z1, tau};
}

std::string MediationWorker::format_results(
    const std::string& exposure_col,
    const std::string& mediator_col,
    const std::string& outcome_col,
    const std::vector<BootstrapResult>& results,
    const BootstrapResult* t0) {
    std::string combination = exposure_col + "_" + mediator_col + "_" + outcome_col;

    std::vector<double> d0_samples, d1_samples, z0_samples, z1_samples, tau_samples;
    d0_samples.reserve(nrep);
    d1_samples.reserve(nrep);
    z0_samples.reserve(nrep);
    z1_samples.reserve(nrep);
    tau_samples.reserve(nrep);

    for (const auto& result : results) {
        d0_samples.push_back(result.indirect_effect_0);  // ACME(control)
        d1_samples.push_back(result.indirect_effect_1);  // ACME(treated)
        z0_samples.push_back(result.direct_effect_0);    // ADE(control)
        z1_samples.push_back(result.direct_effect_1);    // ADE(treated)
        tau_samples.push_back(result.total_effect);      // total effect
    }

    const bool has_t0 = (t0 != nullptr);

    auto d0_stats = has_t0
                        ? calculate_statistics_inplace_with_estimate(
                              d0_samples, t0->indirect_effect_0)
                        : calculate_statistics_inplace(d0_samples);
    auto d1_stats = has_t0
                        ? calculate_statistics_inplace_with_estimate(
                              d1_samples, t0->indirect_effect_1)
                        : calculate_statistics_inplace(d1_samples);
    auto z0_stats = has_t0
                        ? calculate_statistics_inplace_with_estimate(
                              z0_samples, t0->direct_effect_0)
                        : calculate_statistics_inplace(z0_samples);
    auto z1_stats = has_t0
                        ? calculate_statistics_inplace_with_estimate(
                              z1_samples, t0->direct_effect_1)
                        : calculate_statistics_inplace(z1_samples);
    auto tau_stats =
        has_t0 ? calculate_statistics_inplace_with_estimate(tau_samples,
                                                           t0->total_effect)
               : calculate_statistics_inplace(tau_samples);

    std::string result;
    result.reserve(combination.size() + 512);
    result += csv_escape(combination);
    result.push_back(',');
    if (legacy_output_schema) {
        append_statistics(result, d0_stats);
        result.push_back(',');
        append_statistics(result, z0_stats);
        result.push_back(',');
        append_statistics(result, tau_stats);
    } else {
        append_statistics(result, d0_stats);
        result.push_back(',');
        append_statistics(result, d1_stats);
        result.push_back(',');
        append_statistics(result, z0_stats);
        result.push_back(',');
        append_statistics(result, z1_stats);
        result.push_back(',');
        append_statistics(result, tau_stats);
    }
    result.push_back('\n');

    return result;
}

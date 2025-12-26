// mediation_analysis.cpp
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppParallel.h>
#include "glm_fit.h"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <locale>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// [[Rcpp::depends(RcppEigen, RcppParallel)]]

using namespace Rcpp;
using namespace RcppParallel;
using namespace Eigen;

const double X0 = 0.0;
const double X1 = 1.0;

struct StatisticsSummary {
    double mean;
    double percentile_2_5;
    double percentile_97_5;
    double p_value;
};

double mean_cpp(const std::vector<double>& data) {
    if (data.empty()) {
        throw std::runtime_error("Cannot calculate mean of empty vector");
    }
    return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

// R default quantile: type = 7
double quantile_type7(std::vector<double> x, double prob) {
    if (x.empty()) {
        throw std::runtime_error("quantile_type7: empty");
    }
    if (!(prob >= 0.0 && prob <= 1.0)) {
        throw std::invalid_argument("quantile_type7: prob must be in [0,1]");
    }

    std::sort(x.begin(), x.end());
    const size_t n = x.size();
    if (n == 1) {
        return x[0];
    }

    const double h = (static_cast<double>(n) - 1.0) * prob + 1.0; // 1-indexed
    const double hf = std::floor(h);
    size_t j = static_cast<size_t>(hf); // 1..n
    const double g = h - hf;

    if (j <= 1) {
        return x[0];
    }
    if (j >= n) {
        return x[n - 1];
    }

    const size_t idx = j - 1; // 0-index
    return (1.0 - g) * x[idx] + g * x[idx + 1];
}

// mediation::pval-style sign test:
// if estimate == 0 => 1; else p = 2*min(#pos,#neg)/N
double pval_mediate(const std::vector<double>& sims, double estimate) {
    if (sims.empty()) {
        return 1.0;
    }
    if (estimate == 0.0) {
        return 1.0;
    }

    size_t pos = 0;
    size_t neg = 0;
    for (double v : sims) {
        if (v > 0) {
            ++pos;
        } else if (v < 0) {
            ++neg;
        }
    }

    double p = 2.0 * static_cast<double>(std::min(pos, neg)) /
               static_cast<double>(sims.size());
    return (p > 1.0 ? 1.0 : p);
}

static inline double inv_logit(double x) {
    if (x >= 0.0) {
        double z = std::exp(-x);
        return 1.0 / (1.0 + z);
    }
    double z = std::exp(x);
    return z / (1.0 + z);
}

static inline double safe_exp(double x) {
    if (x > 700.0) {
        x = 700.0;
    }
    if (x < -700.0) {
        x = -700.0;
    }
    return std::exp(x);
}

static inline double linkinv(double eta, GlmFamily fam) {
    if (fam == GlmFamily::Gaussian) {
        return eta;
    }
    if (fam == GlmFamily::Binomial) {
        return inv_logit(eta);
    }
    return safe_exp(eta);  // Poisson
}

static inline void simulate_mediator_draw(GlmFamily fam_m,
                                         double eta0,
                                         double eta1,
                                         double sigma2_m,
                                         std::mt19937_64& rng,
                                         VectorXd& M0,
                                         VectorXd& M1) {
    const int n = static_cast<int>(M0.size());
    if (M1.size() != n) {
        throw std::runtime_error("simulate_mediator_draw: size mismatch");
    }

    if (fam_m == GlmFamily::Gaussian) {
        double sd = std::sqrt(std::max(0.0, sigma2_m));
        std::normal_distribution<double> nd(0.0, sd);
        for (int i = 0; i < n; ++i) {
            double e = nd(rng);  // SAME noise for both conditions
            M0[i] = eta0 + e;
            M1[i] = eta1 + e;
        }
    } else if (fam_m == GlmFamily::Binomial) {
        double p0 = inv_logit(eta0);
        double p1 = inv_logit(eta1);
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        for (int i = 0; i < n; ++i) {
            M0[i] = (unif(rng) < p0) ? 1.0 : 0.0;
            M1[i] = (unif(rng) < p1) ? 1.0 : 0.0;
        }
    } else {  // Poisson
        double lam0 = safe_exp(eta0);
        double lam1 = safe_exp(eta1);
        std::poisson_distribution<int> p0(lam0);
        std::poisson_distribution<int> p1(lam1);
        for (int i = 0; i < n; ++i) {
            M0[i] = static_cast<double>(p0(rng));
            M1[i] = static_cast<double>(p1(rng));
        }
    }
}

uint64_t derive_seed(uint64_t base_seed,
                     uint64_t global_combination_idx,
                     uint64_t rep_idx) {
    uint64_t x =
        base_seed ^ (global_combination_idx << 32) ^ rep_idx;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

std::string csv_escape(const std::string& field) {
    bool needs_escape = false;
    for (char c : field) {
        if (c == ',' || c == '"' || c == '\n' || c == '\r') {
            needs_escape = true;
            break;
        }
    }
    if (!needs_escape) {
        return field;
    }

    std::string escaped = "\"";
    for (char c : field) {
        if (c == '"') {
            escaped += "\"\"";
        } else {
            escaped += c;
        }
    }
    escaped += "\"";
    return escaped;
}

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

StatisticsSummary calculate_statistics(const std::vector<double>& samples) {
    const double est = mean_cpp(samples);
    const double alpha = 0.05;
    return {est,
            quantile_type7(samples, alpha / 2.0),
            quantile_type7(samples, 1.0 - alpha / 2.0),
            pval_mediate(samples, est)};
}

Eigen::VectorXd linear_regression(const Eigen::MatrixXd& X,
                                  const Eigen::VectorXd& y) {
    if (X.rows() != y.rows()) {
        throw std::runtime_error("Mismatch in number of rows between X and y");
    }
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
            X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    return svd.solve(y);
}

class BufferedFileWriter {
private:
    std::ofstream file;
    std::vector<std::string> buffer;
    std::mutex mutex;
    const size_t max_buffer_size;
    
public:
    BufferedFileWriter(const std::string& filename, size_t buffer_size = 100)
        : file(filename, std::ios::out | std::ios::trunc),
          max_buffer_size(buffer_size) {
        if (!file.is_open()) {
            throw std::runtime_error("Unable to open output file");
        }
        buffer.reserve(max_buffer_size);
    }
    
    ~BufferedFileWriter() {
        try {
            flush();
        } catch (...) {
            // Ignore exceptions in destructor
        }
        file.close();
    }
    
    void write(const std::string& data) {
        std::lock_guard<std::mutex> lock(mutex);
        buffer.push_back(data);
        if (buffer.size() >= max_buffer_size) {
            flush_internal();
        }
    }
    
    void flush() {
        std::lock_guard<std::mutex> lock(mutex);
        flush_internal();
    }
    
private:
    void flush_internal() {
        for (const auto& line : buffer) {
            file << line;
        }
        file.flush();
        buffer.clear();
    }
};

struct BootstrapResult {
    double indirect_effect_0, indirect_effect_1;
    double direct_effect_0, direct_effect_1;
    double total_effect;
};

class MediationWorker : public Worker {
private:
    const MatrixXd& data;
    const std::unordered_map<std::string, int>& column_index_map;
    const int nrep;
    const std::vector<std::string>& exposure_vec;
    const std::vector<std::string>& mediator_vec;
    const std::vector<std::string>& outcome_vec;
    const std::vector<uint64_t>& global_indices;
    std::string pert_method;
    const std::string mediator_family;
    const std::string outcome_family;
    const bool replace_outcome;
    const uint64_t base_seed;
    const size_t chunk_begin;
    std::vector<std::string>& output_lines;
    std::exception_ptr* error;
    std::mutex* error_mutex;
    
public:
    MediationWorker(const MatrixXd& data_,
                    const std::unordered_map<std::string, int>& column_index_map_,
                    int nrep_,
                    const std::vector<std::string>& exposure_vec_,
                    const std::vector<std::string>& mediator_vec_,
                    const std::vector<std::string>& outcome_vec_,
                    const std::vector<uint64_t>& global_indices_,
                    const std::string& pert_method_,
                    const std::string& mediator_family_,
                    const std::string& outcome_family_,
                    bool replace_outcome_,
                    uint64_t base_seed_,
                    size_t chunk_begin_,
                    std::vector<std::string>& output_lines_,
                    std::exception_ptr* error_,
                    std::mutex* error_mutex_)
        : data(data_),
          column_index_map(column_index_map_),
          nrep(nrep_),
          exposure_vec(exposure_vec_),
          mediator_vec(mediator_vec_),
          outcome_vec(outcome_vec_),
          global_indices(global_indices_),
          pert_method(pert_method_),
          mediator_family(mediator_family_),
          outcome_family(outcome_family_),
          replace_outcome(replace_outcome_),
          base_seed(base_seed_),
          chunk_begin(chunk_begin_),
          output_lines(output_lines_),
          error(error_),
          error_mutex(error_mutex_) {}
    
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t idx = begin; idx < end; ++idx) {
            if (idx >= exposure_vec.size()) {
                break;
            }

            try {
                output_lines[idx - chunk_begin] = process_combination(idx);
            } catch (...) {
                std::lock_guard<std::mutex> lock(*error_mutex);
                if (!*error) {
                    *error = std::current_exception();
                }
                break;
            }
        }
    }
    
private:
    std::string process_combination(std::size_t idx) {
        uint64_t global_combination_idx = global_indices.at(idx);

        std::string exposure_col = exposure_vec[idx];
        std::string mediator_col = mediator_vec[idx];
        std::string outcome_col = outcome_vec[idx];
        
        int exp_idx = get_column_index(exposure_col);
        int med_idx = get_column_index(mediator_col);
        int out_idx = get_column_index(outcome_col);
        
        VectorXd exposure = data.col(exp_idx);
        VectorXd mediator = data.col(med_idx);
        VectorXd outcome = data.col(out_idx);
        
        const int n = static_cast<int>(exposure.size());
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

        GlmFamily fam_m = parse_family_or_auto(mediator_family, mediator);
        GlmFamily fam_y = parse_family_or_auto(outcome_family, outcome);

        VectorXd prior_w = VectorXd::Ones(n);
        VectorXd off_m = VectorXd::Zero(n);
        VectorXd off_y = VectorXd::Zero(n);

        // Fit mediator and outcome models
        MatrixXd X_med(n, 2);
        X_med.col(0).setOnes();
        X_med.col(1) = exposure;
        
        MatrixXd X_out(n, 3);
        X_out.col(0).setOnes();
        X_out.col(1) = mediator;
        X_out.col(2) = exposure;

        const int glm_maxit = 25;
        const double glm_epsilon = 1e-8;
        const double glm_qr_tol = 1e-12;

        GlmFit fit_m;
        try {
            fit_m = glm_fit_irls_qr(X_med,
                                    mediator,
                                    fam_m,
                                    prior_w,
                                    off_m,
                                    glm_maxit,
                                    glm_epsilon,
                                    glm_qr_tol);
        } catch (const std::exception& e) {
            std::string msg = e.what();
            if (msg.find("rank deficient") != std::string::npos) {
                throw std::runtime_error(
                    "Cholesky decomposition failed for mediator model design matrix");
            }
            throw;
        }

        GlmFit fit_y;
        try {
            fit_y = glm_fit_irls_qr(X_out,
                                    outcome,
                                    fam_y,
                                    prior_w,
                                    off_y,
                                    glm_maxit,
                                    glm_epsilon,
                                    glm_qr_tol);
        } catch (const std::exception& e) {
            std::string msg = e.what();
            if (msg.find("rank deficient") != std::string::npos) {
                throw std::runtime_error(
                    "Cholesky decomposition failed for outcome model design matrix");
            }
            throw;
        }
        
        std::vector<BootstrapResult> bootstrap_results;
        
        if (pert_method == "asymptotic") {
            bootstrap_results = perform_bootstrap_asymptotic(
                fit_m,
                fit_y,
                fam_m,
                fam_y,
                n,
                global_combination_idx,
                exposure,
                outcome,
                replace_outcome);
        } else if (pert_method == "bootstrap") {
            bootstrap_results =
                perform_bootstrap_resample(exposure,
                                           mediator,
                                           outcome,
                                           fam_m,
                                           fam_y,
                                           n,
                                           global_combination_idx,
                                           replace_outcome);
        } else {
            throw std::invalid_argument("Unknown perturbation method: " +
                                        pert_method);
        }
        
        return format_results(exposure_col, mediator_col, outcome_col,
                              bootstrap_results);
    }
    
    int get_column_index(const std::string& col_name) {
        auto it = column_index_map.find(col_name);
        if (it == column_index_map.end()) {
            throw std::runtime_error("Column not found: " + col_name);
        }
        return it->second;
    }
    
    std::vector<BootstrapResult> perform_bootstrap_asymptotic(
            const GlmFit& fit_m,
            const GlmFit& fit_y,
            GlmFamily fam_m,
            GlmFamily fam_y,
            int n,
            uint64_t global_combination_idx,
            const VectorXd& exposure_obs,
            const VectorXd& outcome_obs,
            bool replace_outcome_) {
        std::vector<BootstrapResult> results;
        results.reserve(nrep);

        MatrixXd Lm = cholesky_lower_or_throw(fit_m.vcov, "mediator model vcov");
        MatrixXd Ly = cholesky_lower_or_throw(fit_y.vcov, "outcome model vcov");

        VectorXd M0(n), M1(n);
        std::normal_distribution<double> stdnorm(0.0, 1.0);
        const double sigma2_m = (fam_m == GlmFamily::Gaussian) ? fit_m.dispersion : 1.0;

        for (int rep_idx = 0; rep_idx < nrep; ++rep_idx) {
            std::mt19937_64 rng(
                derive_seed(base_seed, global_combination_idx, rep_idx));

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
            results.push_back(simulate_effect_draw(
                bm,
                by,
                fam_m,
                fam_y,
                sigma2_m,
                rng,
                exposure_obs,
                outcome_obs,
                replace_outcome_,
                M0,
                M1));
        }
        
        return results;
    }
    
    std::vector<BootstrapResult> perform_bootstrap_resample(
            const VectorXd& exposure, const VectorXd& mediator,
            const VectorXd& outcome,
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

        VectorXd prior_w = VectorXd::Ones(n);
        VectorXd off_m = VectorXd::Zero(n);
        VectorXd off_y = VectorXd::Zero(n);

        VectorXd boot_exposure(n);
        VectorXd boot_mediator(n);
        VectorXd boot_outcome(n);

        MatrixXd X_med_boot(n, 2);
        MatrixXd X_out_boot(n, 3);

        VectorXd M0(n), M1(n);

        while (rep_idx < nrep) {
            if (total_attempts >= MAX_ATTEMPTS) {
                std::ostringstream oss;
                oss << "Bootstrap failed: " << (total_attempts - rep_idx)
                    << " failures, " << rep_idx << " successes after "
                    << MAX_ATTEMPTS << " attempts. Last error: "
                    << last_exception_msg;
                throw std::runtime_error(oss.str());
            }

            std::mt19937_64 rng(
                derive_seed(base_seed, global_combination_idx, rep_idx));

            bool success = false;
            while (!success) {
                if (total_attempts >= MAX_ATTEMPTS) {
                    std::ostringstream oss;
                    oss << "Bootstrap failed: " << (total_attempts - rep_idx)
                        << " failures, " << rep_idx << " successes after "
                        << MAX_ATTEMPTS << " attempts. Last error: "
                        << last_exception_msg;
                    throw std::runtime_error(oss.str());
                }
                ++total_attempts;

                for (int m = 0; m < n; ++m) {
                    int sample_idx = dis(rng);
                    boot_exposure[m] = exposure[sample_idx];
                    boot_mediator[m] = mediator[sample_idx];
                    boot_outcome[m] = outcome[sample_idx];
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
                    GlmFit fit_m_boot = glm_fit_irls_qr(X_med_boot,
                                                       boot_mediator,
                                                       fam_m,
                                                       prior_w,
                                                       off_m,
                                                       glm_maxit,
                                                       glm_epsilon,
                                                       glm_qr_tol);
                    GlmFit fit_y_boot = glm_fit_irls_qr(X_out_boot,
                                                       boot_outcome,
                                                       fam_y,
                                                       prior_w,
                                                       off_y,
                                                       glm_maxit,
                                                       glm_epsilon,
                                                       glm_qr_tol);

                    const double sigma2_m =
                        (fam_m == GlmFamily::Gaussian) ? fit_m_boot.dispersion : 1.0;

                    results.push_back(simulate_effect_draw(fit_m_boot.coef,
                                                          fit_y_boot.coef,
                                                          fam_m,
                                                          fam_y,
                                                          sigma2_m,
                                                          rng,
                                                          boot_exposure,
                                                          boot_outcome,
                                                          replace_outcome_,
                                                          M0,
                                                          M1));
                    success = true;
                    ++rep_idx;
                } catch (const std::exception& e) {
                    last_exception_msg = e.what();
                }
            }
        }
        
        return results;
    }

    BootstrapResult simulate_effect_draw(const VectorXd& bm,
                                        const VectorXd& by,
                                        GlmFamily fam_m,
                                        GlmFamily fam_y,
                                        double sigma2_m,
                                        std::mt19937_64& rng,
                                        const VectorXd& exposure_obs,
                                        const VectorXd& outcome_obs,
                                        bool replace_outcome_,
                                        VectorXd& M0,
                                        VectorXd& M1) {
        const int n = static_cast<int>(M0.size());
        if (M1.size() != n) {
            throw std::runtime_error("simulate_effect_draw: size mismatch");
        }
        if (exposure_obs.size() != n || outcome_obs.size() != n) {
            throw std::runtime_error("simulate_effect_draw: observed size mismatch");
        }

        double eta_m0 = bm[0] + bm[1] * X0;
        double eta_m1 = bm[0] + bm[1] * X1;
        simulate_mediator_draw(fam_m, eta_m0, eta_m1, sigma2_m, rng, M0, M1);

        double sum_y00 = 0.0;
        double sum_y01 = 0.0;
        double sum_y10 = 0.0;
        double sum_y11 = 0.0;

        for (int i = 0; i < n; ++i) {
            double eta00 = by[0] + by[1] * M0[i] + by[2] * X0;
            double eta01 = by[0] + by[1] * M1[i] + by[2] * X0;
            double eta10 = by[0] + by[1] * M0[i] + by[2] * X1;
            double eta11 = by[0] + by[1] * M1[i] + by[2] * X1;

            double y00 = linkinv(eta00, fam_y);
            double y01 = linkinv(eta01, fam_y);
            double y10 = linkinv(eta10, fam_y);
            double y11 = linkinv(eta11, fam_y);

            if (replace_outcome_) {
                if (exposure_obs[i] == X0) {
                    y00 = outcome_obs[i];
                }
                if (exposure_obs[i] == X1) {
                    y11 = outcome_obs[i];
                }
            }

            sum_y00 += y00;
            sum_y01 += y01;
            sum_y10 += y10;
            sum_y11 += y11;
        }

        double mean_y00 = sum_y00 / static_cast<double>(n);
        double mean_y01 = sum_y01 / static_cast<double>(n);
        double mean_y10 = sum_y10 / static_cast<double>(n);
        double mean_y11 = sum_y11 / static_cast<double>(n);

        double d0 = mean_y01 - mean_y00;  // ACME(control)
        double d1 = mean_y11 - mean_y10;  // ACME(treated)
        double z0 = mean_y10 - mean_y00;  // ADE(control)
        double z1 = mean_y11 - mean_y01;  // ADE(treated)
        double tau = mean_y11 - mean_y00; // total

        return {d0, d1, z0, z1, tau};
    }
    
    std::string format_results(const std::string& exposure_col,
                               const std::string& mediator_col,
                               const std::string& outcome_col,
                               const std::vector<BootstrapResult>& results) {
        std::string combination =
            exposure_col + "_" + mediator_col + "_" + outcome_col;
        
        std::vector<double> acme_samples, ade_samples, total_samples;
        acme_samples.reserve(nrep);
        ade_samples.reserve(nrep);
        total_samples.reserve(nrep);
        
        for (const auto& result : results) {
            acme_samples.push_back(result.indirect_effect_0);  // ACME
            ade_samples.push_back(result.direct_effect_0);     // ADE
            total_samples.push_back(result.total_effect);      // Total effect
        }
        
        auto acme_stats = calculate_statistics(acme_samples);
        auto ade_stats = calculate_statistics(ade_samples);
        auto total_stats = calculate_statistics(total_samples);
        
        std::stringstream result_stream;
        result_stream.imbue(std::locale::classic());
        result_stream << std::fixed << std::setprecision(6);
        result_stream << csv_escape(combination) << ",";
        write_statistics(result_stream, acme_stats);
        result_stream << ",";
        write_statistics(result_stream, ade_stats);
        result_stream << ",";
        write_statistics(result_stream, total_stats);
        result_stream << "\n";
        
        return result_stream.str();
    }
    
    void write_statistics(std::stringstream& stream,
                          const StatisticsSummary& stats) {
        stream << stats.mean << "," << stats.percentile_2_5 << ","
               << stats.percentile_97_5 << "," << stats.p_value;
    }
};

// [[Rcpp::export]]
void mediation_analysis_cpp(NumericMatrix data, CharacterVector column_names,
                            DataFrame combinations, int nrep,
                            std::string output_file,
                            std::string pert = "asymptotic",
                            uint64_t base_seed = 0,
                            bool append = false,
                            std::string mediator_family = "auto",
                            std::string outcome_family = "auto",
                            bool replace_outcome = false) {
    try {
        if (data.nrow() == 0 || data.ncol() == 0) {
            throw std::invalid_argument("Data matrix is empty");
        }
        if (column_names.size() != data.ncol()) {
            throw std::invalid_argument(
                    "Column names size does not match data columns");
        }
        if (combinations.nrows() == 0) {
            throw std::invalid_argument("Combinations dataframe is empty");
        }
        if (nrep <= 0) {
            throw std::invalid_argument(
                    "Number of bootstrap replicates must be positive");
        }
        if (pert != "asymptotic" && pert != "bootstrap") {
            throw std::invalid_argument(
                    "Invalid perturbation method. Use 'asymptotic' or 'bootstrap'");
        }
        
        MatrixXd data_eigen = as<MatrixXd>(data);
        
        std::vector<std::string> column_names_cpp =
            as<std::vector<std::string>>(column_names);
        
        std::vector<std::string> exposure_vec_cpp =
            as<std::vector<std::string>>(combinations["exposure"]);
        std::vector<std::string> mediator_vec_cpp =
            as<std::vector<std::string>>(combinations["mediator"]);
        std::vector<std::string> outcome_vec_cpp =
            as<std::vector<std::string>>(combinations["outcome"]);

        std::vector<uint64_t> global_indices;
        global_indices.reserve(combinations.nrows());
        if (combinations.containsElementNamed("global_idx")) {
            SEXP global_idx_col = combinations["global_idx"];
            if (TYPEOF(global_idx_col) == INTSXP) {
                IntegerVector idx = combinations["global_idx"];
                for (int v : idx) {
                    if (v == NA_INTEGER) {
                        throw std::invalid_argument("global_idx contains NA");
                    }
                    if (v < 0) {
                        throw std::invalid_argument(
                            "global_idx must be non-negative");
                    }
                    global_indices.push_back(static_cast<uint64_t>(v));
                }
            } else {
                NumericVector idx = combinations["global_idx"];
                for (double v : idx) {
                    if (NumericVector::is_na(v)) {
                        throw std::invalid_argument("global_idx contains NA");
                    }
                    if (v < 0.0) {
                        throw std::invalid_argument(
                            "global_idx must be non-negative");
                    }
                    double iv;
                    if (std::modf(v, &iv) != 0.0) {
                        throw std::invalid_argument(
                            "global_idx must contain integers");
                    }
                    global_indices.push_back(static_cast<uint64_t>(iv));
                }
            }
        } else {
            global_indices.resize(combinations.nrows());
            std::iota(global_indices.begin(), global_indices.end(), 0);
        }

        std::unordered_map<std::string, int> column_index_map;
        column_index_map.reserve(column_names_cpp.size());
        for (size_t c = 0; c < column_names_cpp.size(); ++c) {
            column_index_map[column_names_cpp[c]] = static_cast<int>(c);
        }
        
        std::ofstream output_stream;
        if (append) {
            output_stream.open(output_file, std::ios::out | std::ios::app);
        } else {
            output_stream.open(output_file, std::ios::out | std::ios::trunc);
        }
        if (!output_stream.is_open()) {
            throw std::runtime_error("Failed to open output file: " + output_file);
        }
        
        if (!append) {
            std::string header =
                "Combination,ACME_Mean,ACME_2.5%,ACME_97.5%,ACME_p-value,"
                "ADE_Mean,ADE_2.5%,ADE_97.5%,ADE_p-value,"
                "Total_Effect_Mean,Total_Effect_2.5%,Total_Effect_97.5%,Total_"
                "Effect_p-value\n";
            output_stream << header;
        }
        
        const size_t num_combinations = combinations.nrows();
        const size_t chunk_size = 1024;

        for (size_t chunk_begin = 0; chunk_begin < num_combinations;
             chunk_begin += chunk_size) {
            size_t chunk_end = std::min(chunk_begin + chunk_size, num_combinations);
            std::vector<std::string> chunk_results(chunk_end - chunk_begin);

            std::exception_ptr error = nullptr;
            std::mutex error_mutex;

            MediationWorker worker(data_eigen, column_index_map, nrep,
                                   exposure_vec_cpp, mediator_vec_cpp,
                                   outcome_vec_cpp, global_indices, pert,
                                   mediator_family, outcome_family,
                                   replace_outcome, base_seed, chunk_begin,
                                   chunk_results, &error, &error_mutex);

            RcppParallel::parallelFor(chunk_begin, chunk_end, worker);

            if (error) {
                std::rethrow_exception(error);
            }

            for (const auto& line : chunk_results) {
                output_stream << line;
            }
        }
        output_stream.flush();
        
    } catch (const std::exception& e) {
        Rcpp::stop("Error in mediation_analysis_cpp: %s", e.what());
    }
}

// [[Rcpp::export]]
double fastmed_test_p_value_cpp(NumericVector samples) {
    std::vector<double> vec = as<std::vector<double>>(samples);
    if (vec.empty()) {
        return 1.0;
    }
    const double est = mean_cpp(vec);
    return pval_mediate(vec, est);
}

// [[Rcpp::export]]
List fastmed_test_calculate_statistics_cpp(NumericVector samples) {
    std::vector<double> vec = as<std::vector<double>>(samples);
    StatisticsSummary stats = calculate_statistics(vec);
    return List::create(_["mean"] = stats.mean,
                        _["percentile_2_5"] = stats.percentile_2_5,
                        _["percentile_97_5"] = stats.percentile_97_5,
                        _["p_value"] = stats.p_value);
}

// [[Rcpp::export]]
IntegerMatrix fastmed_test_two_bootstrap_samples(int n,
                                                 uint64_t base_seed,
                                                 uint64_t global_combination_idx,
                                                 uint64_t rep_idx) {
    if (n <= 0) {
        throw std::invalid_argument("n must be positive");
    }

    std::mt19937_64 gen(derive_seed(base_seed, global_combination_idx, rep_idx));
    std::uniform_int_distribution<int> dis(0, n - 1);

    IntegerMatrix draws(n, 2);
    for (int i = 0; i < n; ++i) {
        draws(i, 0) = dis(gen);
    }
    for (int i = 0; i < n; ++i) {
        draws(i, 1) = dis(gen);
    }
    return draws;
}

// [[Rcpp::export]]
List fastmed_test_ols_sigmas_cpp(NumericVector exposure,
                                NumericVector mediator,
                                NumericVector outcome) {
    const auto n = static_cast<int>(exposure.size());
    if (mediator.size() != n || outcome.size() != n) {
        throw std::invalid_argument(
            "Mismatch in lengths between exposure, mediator, and outcome");
    }
    if (n <= 0) {
        throw std::invalid_argument("Vectors must be non-empty");
    }

    VectorXd exposure_eigen = as<VectorXd>(exposure);
    VectorXd mediator_eigen = as<VectorXd>(mediator);
    VectorXd outcome_eigen = as<VectorXd>(outcome);

    MatrixXd X_med(n, 2);
    X_med.col(0).setOnes();
    X_med.col(1) = exposure_eigen;
    VectorXd beta_med = linear_regression(X_med, mediator_eigen);

    MatrixXd X_out(n, 3);
    X_out.col(0).setOnes();
    X_out.col(1) = mediator_eigen;
    X_out.col(2) = exposure_eigen;
    VectorXd beta_out = linear_regression(X_out, outcome_eigen);

    const int p_med = static_cast<int>(X_med.cols());
    const int p_out = static_cast<int>(X_out.cols());
    if (n <= p_med) {
        throw std::invalid_argument("Insufficient observations for mediator model");
    }
    if (n <= p_out) {
        throw std::invalid_argument("Insufficient observations for outcome model");
    }

    VectorXd resid_med = mediator_eigen - X_med * beta_med;
    double sigma_med =
        std::sqrt(resid_med.squaredNorm() / static_cast<double>(n - p_med));

    VectorXd resid_out = outcome_eigen - X_out * beta_out;
    double sigma_out =
        std::sqrt(resid_out.squaredNorm() / static_cast<double>(n - p_out));

    return List::create(_["sigma_med"] = sigma_med,
                        _["sigma_out"] = sigma_out,
                        _["df_med"] = n - p_med,
                        _["df_out"] = n - p_out);
}

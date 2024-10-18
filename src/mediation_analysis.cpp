// mediation_analysis.cpp
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppParallel.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
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

const double PERCENTILE_2_5 = 2.5;
const double PERCENTILE_97_5 = 97.5;
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

double std_dev_cpp(const std::vector<double>& data) {
    if (data.size() < 2) {
        throw std::runtime_error(
                "Cannot calculate standard deviation with fewer than 2 elements");
    }
    
    double m = mean_cpp(data);
    double accum = std::accumulate(
        data.begin(), data.end(), 0.0,
        [m](double sum, double val) { return sum + (val - m) * (val - m); });
    return std::sqrt(accum / (data.size() - 1));
}

double percentile_cpp(const std::vector<double>& data, double perc) {
    if (data.empty()) {
        throw std::runtime_error("Cannot calculate percentile of empty vector");
    }
    std::vector<double> sorted_data = data;
    size_t n = sorted_data.size();
    size_t k = static_cast<size_t>((n - 1) * perc / 100.0);
    std::nth_element(sorted_data.begin(), sorted_data.begin() + k,
                     sorted_data.end());
    return sorted_data[k];
}

double p_value_cpp(const std::vector<double>& data, double estimate) {
    if (data.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double count = static_cast<double>(std::count_if(
        data.begin(), data.end(),
        [estimate](double v) { return std::abs(v) >= std::abs(estimate); }));
    return count / data.size();
}

StatisticsSummary calculate_statistics(const std::vector<double>& samples) {
    return {mean_cpp(samples), percentile_cpp(samples, PERCENTILE_2_5),
            percentile_cpp(samples, PERCENTILE_97_5),
            p_value_cpp(samples, 0.0)};
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
    const std::vector<std::string>& column_names;
    const int nrep;
    const std::vector<std::string>& exposure_vec;
    const std::vector<std::string>& mediator_vec;
    const std::vector<std::string>& outcome_vec;
    
    std::unordered_map<std::string, int> column_index_map;
    std::string pert_method;
    std::shared_ptr<BufferedFileWriter> writer;
    
public:
    MediationWorker(const MatrixXd& data_,
                    const std::vector<std::string>& column_names_, int nrep_,
                    const std::vector<std::string>& exposure_vec_,
                    const std::vector<std::string>& mediator_vec_,
                    const std::vector<std::string>& outcome_vec_,
                    const std::string& pert_method_,
                    std::shared_ptr<BufferedFileWriter> writer_)
        : data(data_),
          column_names(column_names_),
          nrep(nrep_),
          exposure_vec(exposure_vec_),
          mediator_vec(mediator_vec_),
          outcome_vec(outcome_vec_),
          pert_method(pert_method_),
          writer(writer_) {
        for (size_t c = 0; c < column_names.size(); ++c) {
            column_index_map[column_names[c]] = c;
        }
    }
    
    void operator()(std::size_t begin, std::size_t end) {
        try {
            std::mt19937 gen(std::random_device{}());
            std::uniform_int_distribution<> dis(0, data.rows() - 1);
            
            for (std::size_t idx = begin; idx < end; ++idx) {
                if (idx >= exposure_vec.size()) break;
                
                std::string result = process_combination(idx, gen, dis);
                writer->write(result);
            }
        } catch (const std::exception& e) {
            Rcpp::stop(std::string("Error in worker: ") + e.what());
        }
    }
    
private:
    std::string process_combination(std::size_t idx, std::mt19937& gen,
                                    std::uniform_int_distribution<>& dis) {
        std::string exposure_col = exposure_vec[idx];
        std::string mediator_col = mediator_vec[idx];
        std::string outcome_col = outcome_vec[idx];
        
        int exp_idx = get_column_index(exposure_col);
        int med_idx = get_column_index(mediator_col);
        int out_idx = get_column_index(outcome_col);
        
        VectorXd exposure = data.col(exp_idx);
        VectorXd mediator = data.col(med_idx);
        VectorXd outcome = data.col(out_idx);
        
        // Fit initial models
        MatrixXd X_med(exposure.size(), 2);
        X_med.col(0).setOnes();
        X_med.col(1) = exposure;
        VectorXd beta_med = linear_regression(X_med, mediator);
        
        MatrixXd X_out(exposure.size(), 3);
        X_out.col(0).setOnes();
        X_out.col(1) = mediator;
        X_out.col(2) = exposure;
        VectorXd beta_out = linear_regression(X_out, outcome);
        
        std::vector<BootstrapResult> bootstrap_results;
        
        if (pert_method == "asymptotic") {
            bootstrap_results = perform_bootstrap_asymptotic(
                X_med, X_out, beta_med, beta_out, mediator, outcome, gen);
        } else if (pert_method == "bootstrap") {
            bootstrap_results = perform_bootstrap_resample(exposure, mediator,
                                                           outcome, gen, dis);
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
            const MatrixXd& X_med, const MatrixXd& X_out, const VectorXd& beta_med,
            const VectorXd& beta_out, const VectorXd& mediator,
            const VectorXd& outcome, std::mt19937& gen) {
        std::vector<BootstrapResult> results;
        results.reserve(nrep);
        
        VectorXd resid_med = mediator - X_med * beta_med;
        double sigma_med = std_dev_cpp(std::vector<double>(
            resid_med.data(), resid_med.data() + resid_med.size()));
        
        VectorXd resid_out = outcome - X_out * beta_out;
        double sigma_out = std_dev_cpp(std::vector<double>(
            resid_out.data(), resid_out.data() + resid_out.size()));
        
        MatrixXd cov_beta_med =
            sigma_med * sigma_med * (X_med.transpose() * X_med).inverse();
        MatrixXd cov_beta_out =
            sigma_out * sigma_out * (X_out.transpose() * X_out).inverse();
        
        Eigen::LLT<MatrixXd> chol_med(cov_beta_med);
        Eigen::LLT<MatrixXd> chol_out(cov_beta_out);
        
        std::normal_distribution<> dist(0.0, 1.0);
        std::normal_distribution<> dist_med_error(0.0, sigma_med);
        std::normal_distribution<> dist_out_error(0.0, sigma_out);
        
        for (int rep = 0; rep < nrep; ++rep) {
            VectorXd beta_med_boot =
                beta_med + chol_med.matrixL() *
                VectorXd::NullaryExpr(beta_med.size(), [&]() {
                    return dist(gen);
                });
            VectorXd beta_out_boot =
                beta_out + chol_out.matrixL() *
                VectorXd::NullaryExpr(beta_out.size(), [&]() {
                    return dist(gen);
                });
            
            double m0 =
                beta_med_boot[0] + beta_med_boot[1] * X0 + dist_med_error(gen);
            double m1 =
                beta_med_boot[0] + beta_med_boot[1] * X1 + dist_med_error(gen);
            
            double y00 = beta_out_boot[0] + beta_out_boot[1] * m0 +
                beta_out_boot[2] * X0 + dist_out_error(gen);
            double y10 = beta_out_boot[0] + beta_out_boot[1] * m1 +
                beta_out_boot[2] * X0 + dist_out_error(gen);
            double y01 = beta_out_boot[0] + beta_out_boot[1] * m0 +
                beta_out_boot[2] * X1 + dist_out_error(gen);
            double y11 = beta_out_boot[0] + beta_out_boot[1] * m1 +
                beta_out_boot[2] * X1 + dist_out_error(gen);
            
            results.push_back({
                y10 - y00,  // indirect_effect_0
                y11 - y01,  // indirect_effect_1
                y01 - y00,  // direct_effect_0
                y11 - y10,  // direct_effect_1
                y11 - y00   // total_effect
            });
        }
        
        return results;
    }
    
    std::vector<BootstrapResult> perform_bootstrap_resample(
            const VectorXd& exposure, const VectorXd& mediator,
            const VectorXd& outcome, std::mt19937& gen,
            std::uniform_int_distribution<>& dis) {
        std::vector<BootstrapResult> results;
        results.reserve(nrep);
        
        auto X_med_boot = std::make_unique<MatrixXd>(exposure.size(), 2);
        auto X_out_boot = std::make_unique<MatrixXd>(exposure.size(), 3);
        
        for (int rep = 0; rep < nrep; ++rep) {
            VectorXd boot_exposure(exposure.size());
            VectorXd boot_mediator(mediator.size());
            VectorXd boot_outcome(outcome.size());
            
            for (int m = 0; m < exposure.size(); ++m) {
                int sample_idx = dis(gen);
                boot_exposure[m] = exposure[sample_idx];
                boot_mediator[m] = mediator[sample_idx];
                boot_outcome[m] = outcome[sample_idx];
            }
            
            X_med_boot->col(0).setOnes();
            X_med_boot->col(1) = boot_exposure;
            
            X_out_boot->col(0).setOnes();
            X_out_boot->col(1) = boot_mediator;
            X_out_boot->col(2) = boot_exposure;
            
            try {
                VectorXd beta_med_boot =
                    linear_regression(*X_med_boot, boot_mediator);
                VectorXd beta_out_boot =
                    linear_regression(*X_out_boot, boot_outcome);
                
                double m0 = beta_med_boot[0] + beta_med_boot[1] * X0;
                double m1 = beta_med_boot[0] + beta_med_boot[1] * X1;
                
                double y00 = beta_out_boot[0] + beta_out_boot[1] * m0 +
                    beta_out_boot[2] * X0;
                double y10 = beta_out_boot[0] + beta_out_boot[1] * m1 +
                    beta_out_boot[2] * X0;
                double y01 = beta_out_boot[0] + beta_out_boot[1] * m0 +
                    beta_out_boot[2] * X1;
                double y11 = beta_out_boot[0] + beta_out_boot[1] * m1 +
                    beta_out_boot[2] * X1;
                
                results.push_back({
                    y10 - y00,  // indirect_effect_0
                    y11 - y01,  // indirect_effect_1
                    y01 - y00,  // direct_effect_0
                    y11 - y10,  // direct_effect_1
                    y11 - y00   // total_effect
                });
            } catch (const std::exception& e) {
                Rcpp::warning(
                    "Linear regression failed: %s Skipping this bootstrap "
                    "replicate.",
                    e.what());
                --rep;
            }
        }
        
        return results;
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
        result_stream << std::fixed << std::setprecision(6);
        result_stream << combination << ",";
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
                            std::string pert = "asymptotic") {
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
        
        auto writer = std::make_shared<BufferedFileWriter>(output_file, 100);
        
        std::string header =
            "Combination,ACME_Mean,ACME_2.5%,ACME_97.5%,ACME_p-value,"
            "ADE_Mean,ADE_2.5%,ADE_97.5%,ADE_p-value,"
            "Total_Effect_Mean,Total_Effect_2.5%,Total_Effect_97.5%,Total_"
            "Effect_p-value\n";
        writer->write(header);
        
        MediationWorker worker(data_eigen, column_names_cpp, nrep,
                               exposure_vec_cpp, mediator_vec_cpp,
                               outcome_vec_cpp, pert, writer);
        
        RcppParallel::parallelFor(0, combinations.nrows(), worker);
        
        writer->flush();
        
    } catch (const std::exception& e) {
        Rcpp::stop("Error in mediation_analysis_cpp: %s", e.what());
    }
}

// mediation_analysis.cpp
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppParallel.h>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <stdexcept>

// [[Rcpp::depends(RcppEigen, RcppParallel)]]

using namespace Rcpp;
using namespace RcppParallel;
using namespace Eigen;

// Helper functions
double mean_cpp(const std::vector<double>& data) {
    if (data.empty()) {
        throw std::runtime_error("Cannot compute mean of empty vector");
    }
    return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

double std_dev_cpp(const std::vector<double>& data) {
    if (data.size() < 2) {
        throw std::runtime_error("Need at least two elements to compute standard deviation");
    }
    double m = mean_cpp(data);
    double accum = 0.0;
    for (const auto& val : data) {
        accum += (val - m) * (val - m);
    }
    return std::sqrt(accum / (data.size() - 1));
}

double percentile_cpp(std::vector<double> data, double perc) {
    if(data.empty()) {
        throw std::runtime_error("Cannot compute percentile of empty vector");
    }
    std::sort(data.begin(), data.end());
    double idx = perc / 100.0 * (data.size() - 1);
    size_t lower = static_cast<size_t>(std::floor(idx));
    size_t upper = static_cast<size_t>(std::ceil(idx));
    if(lower == upper) {
        return data[lower];
    }
    double weight = idx - lower;
    return data[lower] * (1.0 - weight) + data[upper] * weight;
}

double p_value_cpp(const std::vector<double>& data, double estimate) {
    if (data.empty()) {
        throw std::runtime_error("Cannot compute p-value of empty vector");
    }
    size_t count = std::count_if(data.begin(), data.end(), [&](double v) { return v <= estimate; });
    return static_cast<double>(count) / data.size();
}

// Linear regression using Eigen
VectorXd linear_regression(const MatrixXd& X, const VectorXd& y) {
    if (X.rows() != y.rows()) {
        throw std::runtime_error("Mismatch in number of rows between X and y");
    }
    // Using QR decomposition for better numerical stability
    return X.colPivHouseholderQr().solve(y);
}

// Mutex for thread-safe file writing
std::mutex file_mutex;

// Worker class for parallel processing
struct MediationWorker : public Worker {
    // Inputs
    const RMatrix<double> data;
    const std::vector<std::string> column_names;
    const std::vector<std::string> exposure_cols;
    const std::vector<std::string> mediator_cols;
    const std::vector<std::string> outcome_cols;
    const int nrep;
    const Rcpp::DataFrame combinations;
    const std::string output_file;
    const Rcpp::DataFrame combinations;

    // Constructor
    MediationWorker(const NumericMatrix data_,
                    const CharacterVector column_names_,
                    const CharacterVector exposure_cols_,
                    const CharacterVector mediator_cols_,
                    const CharacterVector outcome_cols_,
                    int nrep_,
                    std::string output_file_,
                    Rcpp::DataFrame combinations_)
        : data(data_),
          column_names(Rcpp::as<std::vector<std::string>>(column_names_)),
          exposure_cols(Rcpp::as<std::vector<std::string>>(exposure_cols_)),
          mediator_cols(Rcpp::as<std::vector<std::string>>(mediator_cols_)),
          outcome_cols(Rcpp::as<std::vector<std::string>>(outcome_cols_)),
          nrep(nrep_),
          output_file(output_file_),
          combinations(combinations_) {}
    
    // Function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        try {
            std::ofstream outfile(output_file, std::ios::app);
            if (!outfile.is_open()) {
                throw std::runtime_error("Unable to open output file");
            }
            
            // Random number generator
            std::mt19937 gen(std::random_device{}());
            std::uniform_int_distribution<> dis(0, data.nrow() - 1);
            
            // Get column vectors from combinations DataFrame
            Rcpp::CharacterVector exposure_vec = combinations["exposure"];
            Rcpp::CharacterVector mediator_vec = combinations["mediator"];
            Rcpp::CharacterVector outcome_vec = combinations["outcome"];
            
            // Local buffer to accumulate results
            std::vector<std::string> local_buffer;
            const size_t buffer_size = 100;
            
            
            for (std::size_t idx = begin; idx < end; ++idx) {
                std::string exposure_col = Rcpp::as<std::string>(exposure_vec[idx]);
                std::string mediator_col = Rcpp::as<std::string>(mediator_vec[idx]);
                std::string outcome_col = Rcpp::as<std::string>(outcome_vec[idx]);
                
                // Find column indices
                int exp_idx = -1, med_idx = -1, out_idx = -1;
                for (int c = 0; c < column_names.size(); ++c) {
                    if (column_names[c] == exposure_col) exp_idx = c;
                    if (column_names[c] == mediator_col) med_idx = c;
                    if (column_names[c] == outcome_col) out_idx = c;
                }
                
                if (exp_idx == -1 || med_idx == -1 || out_idx == -1) {
                    throw std::runtime_error("Column not found");
                }
                
                // Extract data for current combination
                VectorXd exposure(data.nrow()), mediator(data.nrow()), outcome(data.nrow());
                for (int row = 0; row < data.nrow(); ++row) {
                    exposure(row) = data(row, exp_idx);
                    mediator(row) = data(row, med_idx);
                    outcome(row) = data(row, out_idx);
                }
                
                // Prepare design matrices
                MatrixXd X_med(data.nrow(), 2);
                X_med.col(0).setOnes();
                X_med.col(1) = exposure;
                
                MatrixXd X_out(data.nrow(), 3);
                X_out.col(0).setOnes();
                X_out.col(1) = mediator;
                X_out.col(2) = exposure;
                
                // Perform initial regressions
                VectorXd beta_med = linear_regression(X_med, mediator);
                VectorXd beta_out = linear_regression(X_out, outcome);
                
                // Bootstrap
                std::vector<double> ind0_samples, ind1_samples, dir0_samples, dir1_samples, tot_samples;
                ind0_samples.reserve(nrep);
                ind1_samples.reserve(nrep);
                dir0_samples.reserve(nrep);
                dir1_samples.reserve(nrep);
                tot_samples.reserve(nrep);
                
                for (int rep = 0; rep < nrep; ++rep) {
                    // Sample with replacement
                    VectorXd boot_exposure(data.nrow()), boot_mediator(data.nrow()), boot_outcome(data.nrow());
                    for (int m = 0; m < data.nrow(); ++m) {
                        int sample_idx = dis(gen);
                        boot_exposure[m] = exposure[sample_idx];
                        boot_mediator[m] = mediator[sample_idx];
                        boot_outcome[m] = outcome[sample_idx];
                    }
                    
                    // Regressions on bootstrap sample
                    VectorXd beta_med_boot = linear_regression(X_med, boot_mediator);
                    VectorXd beta_out_boot = linear_regression(X_out, boot_outcome);
                    
                    // Calculate effects
                    double m0 = beta_med_boot[0];
                    double m1 = beta_med_boot[0] + beta_med_boot[1];
                    double y00 = beta_out_boot[0] + beta_out_boot[1] * m0;
                    double y10 = beta_out_boot[0] + beta_out_boot[1] * m1;
                    double y01 = y00 + beta_out_boot[2];
                    double y11 = y10 + beta_out_boot[2];
                    
                    ind0_samples.push_back(y10 - y00);
                    ind1_samples.push_back(y11 - y01);
                    dir0_samples.push_back(y01 - y00);
                    dir1_samples.push_back(y11 - y10);
                    tot_samples.push_back(y11 - y00);
                }
                
                // Calculate statistics
                std::string combination = exposure_col + "_" + mediator_col + "_" + outcome_col;
                
                // Prepare result string
                std::stringstream result;
                result << std::fixed << std::setprecision(6);
                result << combination << ","
                       << mean_cpp(ind1_samples) << "," << std_dev_cpp(ind1_samples) << ","
                       << percentile_cpp(ind1_samples, 2.5) << "," << percentile_cpp(ind1_samples, 97.5) << ","
                       << p_value_cpp(ind1_samples, 0.0) << ","
                       << mean_cpp(ind0_samples) << "," << std_dev_cpp(ind0_samples) << ","
                       << percentile_cpp(ind0_samples, 2.5) << "," << percentile_cpp(ind0_samples, 97.5) << ","
                       << p_value_cpp(ind0_samples, 0.0) << ","
                       << mean_cpp(dir1_samples) << "," << std_dev_cpp(dir1_samples) << ","
                       << percentile_cpp(dir1_samples, 2.5) << "," << percentile_cpp(dir1_samples, 97.5) << ","
                       << p_value_cpp(dir1_samples, 0.0) << ","
                       << mean_cpp(dir0_samples) << "," << std_dev_cpp(dir0_samples) << ","
                       << percentile_cpp(dir0_samples, 2.5) << "," << percentile_cpp(dir0_samples, 97.5) << ","
                       << p_value_cpp(dir0_samples, 0.0) << ","
                       << mean_cpp(tot_samples) << "," << std_dev_cpp(tot_samples) << ","
                       << percentile_cpp(tot_samples, 2.5) << "," << percentile_cpp(tot_samples, 97.5) << ","
                       << p_value_cpp(tot_samples, 0.0) << "\n";
                
                // Add result to local buffer
                local_buffer.push_back(result.str());
                
                // If buffer is full, write to file
                if (local_buffer.size() >= buffer_size) {
                    // Lock mutex and write buffer to file
                    {
                        std::lock_guard<std::mutex> lock(file_mutex);
                        std::ofstream outfile(output_file, std::ios::app);
                        if (!outfile.is_open()) {
                            throw std::runtime_error("Unable to open output file");
                        }
                        for (const auto& line : local_buffer) {
                            outfile << line;
                        }
                        outfile.close();
                    }
                    // Clear local buffer
                    local_buffer.clear();
                }
            }
            
            // After processing all combinations, write any remaining results in the buffer
            if (!local_buffer.empty()) {
                std::lock_guard<std::mutex> lock(file_mutex);
                std::ofstream outfile(output_file, std::ios::app);
                if (!outfile.is_open()) {
                    throw std::runtime_error("Unable to open output file");
                }
                for (const auto& line : local_buffer) {
                    outfile << line;
                }
                outfile.close();
            }
        } catch (const std::exception& e) {
            Rcpp::stop(std::string("Error in worker: ") + e.what());
        }
    }
};

// [[Rcpp::export]]
void mediation_analysis_cpp(NumericMatrix data,
                            CharacterVector column_names,
                            CharacterVector exposure_cols,
                            CharacterVector mediator_cols,
                            CharacterVector outcome_cols,
                            DataFrame combinations,
                            int nrep,
                            std::string output_file) {
    try {
        MediationWorker worker(data, column_names, exposure_cols, mediator_cols, outcome_cols, nrep, output_file, combinations);
        parallelFor(0, combinations.nrows(), worker);
    } catch (const std::exception& e) {
        Rcpp::stop(std::string("Error in mediation_analysis_cpp: ") + e.what());
    }
}

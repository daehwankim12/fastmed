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
#include <sstream>
#include <cmath>
#include <unordered_map>

// [[Rcpp::depends(RcppEigen, RcppParallel)]]

using namespace Rcpp;
using namespace RcppParallel;
using namespace Eigen;

double mean_cpp(const std::vector<double>& data) {
    if (data.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

double std_dev_cpp(const std::vector<double>& data) {
    if (data.size() < 2) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    double m = mean_cpp(data);
    double accum = std::accumulate(data.begin(), data.end(), 0.0,
                                   [m](double sum, double val) { return sum + (val - m) * (val - m); });
    return std::sqrt(accum / (data.size() - 1));
}

double percentile_cpp(const std::vector<double>& data, double perc) {
    if (data.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::vector<double> sorted_data = data;
    std::sort(sorted_data.begin(), sorted_data.end());
    double idx = perc / 100.0 * (sorted_data.size() - 1);
    size_t lower = static_cast<size_t>(std::floor(idx));
    size_t upper = static_cast<size_t>(std::ceil(idx));
    if (lower == upper) {
        return sorted_data[lower];
    }
    double weight = idx - lower;
    return sorted_data[lower] * (1.0 - weight) + sorted_data[upper] * weight;
}

double p_value_cpp(const std::vector<double>& data, double estimate) {
    if (data.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return static_cast<double>(std::count_if(data.begin(), data.end(), 
                                             [estimate](double v) { return v <= estimate; })) / data.size();
}

VectorXd linear_regression(const MatrixXd& X, const VectorXd& y) {
    if (X.rows() != y.rows()) {
        throw std::runtime_error("Mismatch in number of rows between X and y");
    }
    Eigen::LDLT<MatrixXd> ldlt(X.transpose() * X);
    if (ldlt.rcond() < 1e-10) {
        throw std::runtime_error("Matrix is (near) singular");
    }
    return ldlt.solve(X.transpose() * y);
}

class BufferedFileWriter {
private:
    std::ofstream file;
    std::vector<std::string> buffer;
    std::mutex mutex;
    const size_t max_buffer_size;
    
public:
    BufferedFileWriter(const std::string& filename, size_t buffer_size = 1000)
        : file(filename, std::ios::out | std::ios::trunc), max_buffer_size(buffer_size) {
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
    double ind0, ind1, dir0, dir1, tot;
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
    BufferedFileWriter& writer;
    
public:
    MediationWorker(const MatrixXd& data_,
                    const std::vector<std::string>& column_names_,
                    int nrep_,
                    const std::vector<std::string>& exposure_vec_,
                    const std::vector<std::string>& mediator_vec_,
                    const std::vector<std::string>& outcome_vec_,
                    BufferedFileWriter& writer_)
        : data(data_), column_names(column_names_), nrep(nrep_),
          exposure_vec(exposure_vec_), mediator_vec(mediator_vec_), outcome_vec(outcome_vec_),
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
                if (idx >= exposure_vec.size()) break;  // Ensure we don't go out of bounds
                
                std::string result = process_combination(idx, gen, dis);
                writer.write(result);
            }
        } catch (const std::exception& e) {
            Rcpp::stop(std::string("Error in worker: ") + e.what());
        }
    }
    
private:
    std::string process_combination(std::size_t idx, std::mt19937& gen, std::uniform_int_distribution<>& dis) {
        std::string exposure_col = exposure_vec[idx];
        std::string mediator_col = mediator_vec[idx];
        std::string outcome_col = outcome_vec[idx];
        
        int exp_idx = get_column_index(exposure_col);
        int med_idx = get_column_index(mediator_col);
        int out_idx = get_column_index(outcome_col);
        
        VectorXd exposure = data.col(exp_idx);
        VectorXd mediator = data.col(med_idx);
        VectorXd outcome = data.col(out_idx);
        
        std::vector<BootstrapResult> bootstrap_results = perform_bootstrap(exposure, mediator, outcome, gen, dis);
        
        return format_results(exposure_col, mediator_col, outcome_col, bootstrap_results);
    }
    
    int get_column_index(const std::string& col_name) {
        auto it = column_index_map.find(col_name);
        if (it == column_index_map.end()) {
            throw std::runtime_error("Column not found: " + col_name);
        }
        return it->second;
    }
    
    std::vector<BootstrapResult> perform_bootstrap(const VectorXd& exposure, const VectorXd& mediator, const VectorXd& outcome,
                                                   std::mt19937& gen, std::uniform_int_distribution<>& dis) {
        std::vector<BootstrapResult> results;
        results.reserve(nrep);
        
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
            
            MatrixXd X_med_boot(exposure.size(), 2);
            X_med_boot.col(0).setOnes();
            X_med_boot.col(1) = boot_exposure;
            
            MatrixXd X_out_boot(exposure.size(), 3);
            X_out_boot.col(0).setOnes();
            X_out_boot.col(1) = boot_mediator;
            X_out_boot.col(2) = boot_exposure;
            
            try {
                VectorXd beta_med_boot = linear_regression(X_med_boot, boot_mediator);
                VectorXd beta_out_boot = linear_regression(X_out_boot, boot_outcome);
                
                double m0 = beta_med_boot[0];
                double m1 = beta_med_boot[0] + beta_med_boot[1];
                double y00 = beta_out_boot[0] + beta_out_boot[1] * m0;
                double y10 = beta_out_boot[0] + beta_out_boot[1] * m1;
                double y01 = y00 + beta_out_boot[2];
                double y11 = y10 + beta_out_boot[2];
                
                results.push_back({y10 - y00, y11 - y01, y01 - y00, y11 - y10, y11 - y00});
            } catch (const std::runtime_error& e) {
                Rcpp::warning("Linear regression failed: %s Skipping this bootstrap replicate.", e.what());
                --rep; // Retry this replicate
            }
        }
        
        return results;
    }
    
    std::string format_results(const std::string& exposure_col, const std::string& mediator_col, const std::string& outcome_col,
                               const std::vector<BootstrapResult>& results) {
        std::string combination = exposure_col + "_" + mediator_col + "_" + outcome_col;
        
        std::vector<double> acme_samples, ade_samples, total_samples;
        acme_samples.reserve(nrep);
        ade_samples.reserve(nrep);
        total_samples.reserve(nrep);
        
        for (const auto& result : results) {
            acme_samples.push_back(result.ind0);  // ACME
            ade_samples.push_back(result.dir0);   // ADE
            total_samples.push_back(result.tot);  // Total effect
        }
        
        std::stringstream result_stream;
        result_stream << std::fixed << std::setprecision(6);
        result_stream << combination << ",";
        write_statistics(result_stream, acme_samples);
        result_stream << ",";
        write_statistics(result_stream, ade_samples);
        result_stream << ",";
        write_statistics(result_stream, total_samples);
        result_stream << "\n";
        
        return result_stream.str();
    }
    
    void write_statistics(std::stringstream& stream, const std::vector<double>& samples) {
        stream << mean_cpp(samples) << "," << std_dev_cpp(samples) << ","
               << percentile_cpp(samples, 2.5) << "," << percentile_cpp(samples, 97.5) << ","
               << p_value_cpp(samples, 0.0);
    }
};

// [[Rcpp::export]]
void mediation_analysis_cpp(NumericMatrix data,
                            CharacterVector column_names,
                            DataFrame combinations,
                            int nrep,
                            std::string output_file) {
    try {
        if (data.nrow() == 0 || data.ncol() == 0) {
            throw std::invalid_argument("Data matrix is empty");
        }
        if (column_names.size() != data.ncol()) {
            throw std::invalid_argument("Column names size does not match data columns");
        }
        if (combinations.nrows() == 0) {
            throw std::invalid_argument("Combinations dataframe is empty");
        }
        if (nrep <= 0) {
            throw std::invalid_argument("Number of bootstrap replicates must be positive");
        }
        
        MatrixXd data_eigen = as<MatrixXd>(data);
        std::vector<std::string> column_names_cpp = as<std::vector<std::string>>(column_names);
        
        std::vector<std::string> exposure_vec_cpp = as<std::vector<std::string>>(combinations["exposure"]);
        std::vector<std::string> mediator_vec_cpp = as<std::vector<std::string>>(combinations["mediator"]);
        std::vector<std::string> outcome_vec_cpp = as<std::vector<std::string>>(combinations["outcome"]);
        
        BufferedFileWriter writer(output_file, 1000);
        
        // 열 이름을 직접 파일에 쓰기
        std::string header = "Combination,ACME_Mean,ACME_SD,ACME_2.5%,ACME_97.5%,ACME_p-value,"
        "ADE_Mean,ADE_SD,ADE_2.5%,ADE_97.5%,ADE_p-value,"
        "Total_Effect_Mean,Total_Effect_SD,Total_Effect_2.5%,Total_Effect_97.5%,Total_Effect_p-value\n";
        writer.write(header);
        writer.flush();
        
        MediationWorker worker(data_eigen, column_names_cpp, nrep,
                               exposure_vec_cpp, mediator_vec_cpp, outcome_vec_cpp, writer);
        
        parallelFor(0, combinations.nrows(), worker);
        
        writer.flush();
        
    } catch (const std::exception& e) {
        Rcpp::stop("Error in mediation_analysis_cpp: %s", e.what());
    }
}

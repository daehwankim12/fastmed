#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

struct StatisticsSummary {
    double mean;
    double percentile_2_5;
    double percentile_97_5;
    double p_value;
};

inline double mean_cpp(const std::vector<double>& data) {
    if (data.empty()) {
        throw std::runtime_error("Cannot calculate mean of empty vector");
    }
    return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

inline double quantile_type7_sorted(const std::vector<double>& sorted, double prob) {
    if (sorted.empty()) {
        throw std::runtime_error("quantile_type7_sorted: empty");
    }
    if (!(prob >= 0.0 && prob <= 1.0)) {
        throw std::invalid_argument("quantile_type7_sorted: prob must be in [0,1]");
    }

    const size_t n = sorted.size();
    if (n == 1) {
        return sorted[0];
    }

    const double h = (static_cast<double>(n) - 1.0) * prob + 1.0;  // 1-indexed
    const double hf = std::floor(h);
    size_t j = static_cast<size_t>(hf);  // 1..n
    const double g = h - hf;

    if (j <= 1) {
        return sorted[0];
    }
    if (j >= n) {
        return sorted[n - 1];
    }

    const size_t idx = j - 1;  // 0-index
    return (1.0 - g) * sorted[idx] + g * sorted[idx + 1];
}

inline StatisticsSummary calculate_statistics_inplace(std::vector<double>& samples) {
    if (samples.empty()) {
        throw std::runtime_error("calculate_statistics_inplace: empty");
    }

    double sum = 0.0;
    size_t pos = 0;
    size_t neg = 0;
    for (double v : samples) {
        sum += v;
        if (v > 0.0) {
            ++pos;
        } else if (v < 0.0) {
            ++neg;
        }
    }

    const double est = sum / static_cast<double>(samples.size());
    const double alpha = 0.05;

    std::sort(samples.begin(), samples.end());

    double p = 1.0;
    if (est != 0.0) {
        p = 2.0 * static_cast<double>(std::min(pos, neg)) /
            static_cast<double>(samples.size());
        if (p > 1.0) {
            p = 1.0;
        }
    }

    return {est,
            quantile_type7_sorted(samples, alpha / 2.0),
            quantile_type7_sorted(samples, 1.0 - alpha / 2.0),
            p};
}

// mediation::pval-style sign test:
// if estimate == 0 => 1; else p = 2*min(#pos,#neg)/N
inline double pval_mediate(const std::vector<double>& sims, double estimate) {
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


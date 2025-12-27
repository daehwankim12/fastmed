#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>

inline uint64_t derive_seed(uint64_t base_seed,
                            uint64_t global_combination_idx,
                            uint64_t rep_idx) {
    uint64_t x = base_seed ^ (global_combination_idx << 32) ^ rep_idx;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

class SafePoissonSampler {
public:
    explicit SafePoissonSampler(double lambda) {
        if (!std::isfinite(lambda) || lambda < 0.0) {
            throw std::runtime_error(
                "SafePoissonSampler: lambda must be finite and non-negative");
        }

        if (lambda <= max_exact_lambda()) {
            mode_ = Mode::Exact;
            poisson_.param(std::poisson_distribution<int>::param_type(lambda));
            return;
        }

        mode_ = Mode::NormalApprox;
        const double sd = std::sqrt(lambda);
        normal_.param(std::normal_distribution<double>::param_type(lambda, sd));
    }

    int operator()(std::mt19937_64& rng) {
        if (mode_ == Mode::Exact) {
            return poisson_(rng);
        }

        const double x = normal_(rng);
        if (!std::isfinite(x)) {
            throw std::runtime_error("SafePoissonSampler: draw is not finite");
        }
        if (x <= 0.0) {
            return 0;
        }
        const double max_int = static_cast<double>(std::numeric_limits<int>::max());
        if (x >= max_int) {
            return std::numeric_limits<int>::max();
        }
        return static_cast<int>(std::llround(x));
    }

private:
    enum class Mode { Exact, NormalApprox };

    static constexpr double max_exact_lambda() {
        return 1e6;
    }

    Mode mode_ = Mode::Exact;
    std::poisson_distribution<int> poisson_;
    std::normal_distribution<double> normal_;
};

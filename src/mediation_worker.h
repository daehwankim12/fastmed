#pragma once

#include <RcppEigen.h>
#include <RcppParallel.h>

#include "fastmed_loop_order.h"
#include "glm_fit.h"
#include "match_mediation_rng.h"

#include <cstdint>
#include <random>
#include <string>
#include <vector>

struct BootstrapResult {
    double indirect_effect_0, indirect_effect_1;
    double direct_effect_0, direct_effect_1;
    double total_effect;
};

class MediationWorker : public RcppParallel::Worker {
  private:
    const Eigen::Map<const Eigen::MatrixXd>& data;
    const Eigen::VectorXd& weights;
    const std::vector<std::string>& column_names;
    const int nrep;
    const std::vector<int>& exposure_col_idx;
    const std::vector<int>& mediator_col_idx;
    const std::vector<int>& outcome_col_idx;
    const std::vector<GlmFamily>& mediator_fams;
    const std::vector<GlmFamily>& outcome_fams;
    const std::vector<char>& mediator_fams_ok;
    const std::vector<char>& outcome_fams_ok;
    std::string pert_method;
    const bool replace_outcome;
    const bool legacy_output_schema;
    const bool excel_safe_csv;
    const double treat_value;
    const double control_value;
    const uint64_t base_seed;
    const bool match_mediation;
    const LoopOrder loop_order;
    const bool fail_fast;
    const bool include_failure_reason;
    const size_t chunk_begin;
    std::vector<std::string>& output_lines;
    const MatchMediationRngBlock* match_rng;

  public:
    MediationWorker(const Eigen::Map<const Eigen::MatrixXd>& data_,
                    const Eigen::VectorXd& weights_,
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
                    bool excel_safe_csv_,
                    double treat_value_,
                    double control_value_,
                    uint64_t base_seed_,
                    bool match_mediation_,
                    LoopOrder loop_order_,
                    bool fail_fast_,
                    bool include_failure_reason_,
                    size_t chunk_begin_,
                    std::vector<std::string>& output_lines_,
                    const MatchMediationRngBlock* match_rng_);

    void operator()(std::size_t begin, std::size_t end);

    std::string process_combination_serial(std::size_t idx,
                                           Eigen::MatrixXd& X_med,
                                           Eigen::MatrixXd& X_out,
                                           Eigen::VectorXd& off_m,
                                           Eigen::VectorXd& off_y,
                                           Eigen::VectorXd& outcome_obs_buf,
                                           Eigen::VectorXd& weights_obs_buf);

  private:
    std::string format_na_row(const std::string& exposure_col,
                              const std::string& mediator_col,
                              const std::string& outcome_col,
                              const std::string& failure_reason);

    std::string handle_failure(const std::string& exposure_col,
                               const std::string& mediator_col,
                               const std::string& outcome_col,
                               const std::string& failure_reason);

    std::string process_combination(std::size_t idx,
                                    Eigen::MatrixXd& X_med,
                                    Eigen::MatrixXd& X_out,
                                    Eigen::VectorXd& off_m,
                                    Eigen::VectorXd& off_y,
                                    Eigen::VectorXd& outcome_obs_buf,
                                    Eigen::VectorXd& weights_obs_buf);

    std::vector<BootstrapResult> perform_bootstrap_asymptotic(
        const GlmFit& fit_m,
        const GlmFit& fit_y,
        GlmFamily fam_m,
        GlmFamily fam_y,
        int n,
        uint64_t global_combination_idx,
        const Eigen::Ref<const Eigen::VectorXd>& exposure_obs,
        const Eigen::Ref<const Eigen::VectorXd>& outcome_obs,
        const Eigen::Ref<const Eigen::VectorXd>& weights_obs,
        bool replace_outcome_);

    std::vector<BootstrapResult> perform_bootstrap_asymptotic_mediation_rng(
        const GlmFit& fit_m,
        const GlmFit& fit_y,
        GlmFamily fam_m,
        GlmFamily fam_y,
        int n,
        const Eigen::Ref<const Eigen::VectorXd>& exposure_obs,
        const Eigen::Ref<const Eigen::VectorXd>& outcome_obs,
        const Eigen::Ref<const Eigen::VectorXd>& weights_obs,
        bool replace_outcome_,
        const MatchMediationRngSlice& rng_slice);

    std::vector<BootstrapResult> perform_bootstrap_asymptotic_mediation_poisson_serial(
        const GlmFit& fit_m,
        const GlmFit& fit_y,
        GlmFamily fam_y,
        int n,
        const Eigen::Ref<const Eigen::VectorXd>& exposure_obs,
        const Eigen::Ref<const Eigen::VectorXd>& outcome_obs,
        const Eigen::Ref<const Eigen::VectorXd>& weights_obs,
        bool replace_outcome_);

    std::vector<BootstrapResult> perform_bootstrap_resample(
        const Eigen::Ref<const Eigen::VectorXd>& exposure,
        const Eigen::Ref<const Eigen::VectorXd>& mediator,
        const Eigen::Ref<const Eigen::VectorXd>& outcome,
        GlmFamily fam_m,
        GlmFamily fam_y,
        int n,
        uint64_t global_combination_idx,
        const Eigen::Ref<const Eigen::VectorXd>& weights_obs,
        bool replace_outcome_);

    std::vector<BootstrapResult> perform_bootstrap_resample_mediation_rng(
        const GlmFit& fit_m,
        const GlmFit& fit_y,
        const Eigen::Ref<const Eigen::VectorXd>& exposure,
        const Eigen::Ref<const Eigen::VectorXd>& mediator,
        const Eigen::Ref<const Eigen::VectorXd>& outcome,
        GlmFamily fam_m,
        GlmFamily fam_y,
        int n,
        const Eigen::Ref<const Eigen::VectorXd>& weights_obs,
        bool replace_outcome_,
        BootstrapResult* t0_out,
        const MatchMediationRngSlice& rng_slice);

    std::vector<BootstrapResult> perform_bootstrap_resample_mediation_poisson_serial(
        const GlmFit& fit_m,
        const GlmFit& fit_y,
        const Eigen::Ref<const Eigen::VectorXd>& exposure,
        const Eigen::Ref<const Eigen::VectorXd>& mediator,
        const Eigen::Ref<const Eigen::VectorXd>& outcome,
        GlmFamily fam_y,
        int n,
        const Eigen::Ref<const Eigen::VectorXd>& weights_obs,
        bool replace_outcome_,
        BootstrapResult* t0_out);

    BootstrapResult simulate_effect_draw(
        const Eigen::VectorXd& bm,
        const Eigen::VectorXd& by,
        GlmFamily fam_m,
        GlmFamily fam_y,
        double sigma2_m,
        std::mt19937_64& rng,
        const Eigen::Ref<const Eigen::VectorXd>& exposure_obs,
        const Eigen::Ref<const Eigen::VectorXd>& outcome_obs,
        const Eigen::Ref<const Eigen::VectorXd>& weights_obs,
        bool replace_outcome_);

    std::string format_results(const std::string& exposure_col,
                               const std::string& mediator_col,
                               const std::string& outcome_col,
                               const std::vector<BootstrapResult>& results,
                               const BootstrapResult* t0);
};

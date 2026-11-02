################################################################################
# Module: 05_pipeline_utils.R
# Purpose: Shared helpers for experiment pipelines
################################################################################

source("01_interval_utils.R")
source("01_conformal_utils.R")
source("02_models.R")
source("02_dgp_registry.R")


resolve_propensity_scores <- function(X, T, ps_fit_idx, scenario,
                                      ps_mode = "estimated_ps",
                                      true_ps = NULL,
                                      ...) {
    if (ps_mode == "true_ps") {
        ps_all <- if (!is.null(true_ps)) true_ps else get_true_ps(scenario, X, ...)
        return(list(ps_all = ps_all, model = NULL, mode = ps_mode, fit_idx = ps_fit_idx))
    }

    if (ps_mode == "estimated_ps") {
        ps_model <- fit_propensity_score(X[ps_fit_idx, , drop = FALSE], T[ps_fit_idx])
        ps_all <- predict_propensity_score(ps_model, X)
        return(list(ps_all = ps_all, model = ps_model, mode = ps_mode, fit_idx = ps_fit_idx))
    }

    if (ps_mode == "misspecified_ps") {
        if (!scenario %in% c("S3", "S3_1", "S3_2")) {
            stop("misspecified_ps is only defined for the S3 family.")
        }
        ps_model <- fit_misspecified_ps_S3(X[ps_fit_idx, , drop = FALSE], T[ps_fit_idx])
        ps_all <- predict_misspecified_ps_S3(ps_model, X)
        return(list(ps_all = ps_all, model = ps_model, mode = ps_mode, fit_idx = ps_fit_idx))
    }

    stop("ps_mode must be true_ps, estimated_ps, or misspecified_ps.")
}


split_fold1_for_stage1 <- function(fold1_idx, T, fit_prop = 0.75, seed = NULL) {
    if (fit_prop <= 0 || fit_prop >= 1) {
        stop("fit_prop must be in (0, 1).")
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }

    sample_fit <- function(idx, arm_label) {
        n_arm <- length(idx)
        if (n_arm < 2) {
            stop(sprintf("Need at least two observations in treatment arm %s for Stage 1 split.", arm_label))
        }
        n_fit <- min(max(1, floor(n_arm * fit_prop)), n_arm - 1)
        sample(idx, n_fit)
    }

    treated_idx <- fold1_idx[T[fold1_idx] == 1]
    control_idx <- fold1_idx[T[fold1_idx] == 0]
    fit_treated_idx <- sample_fit(treated_idx, 1)
    fit_control_idx <- sample_fit(control_idx, 0)
    fit_idx <- c(fit_treated_idx, fit_control_idx)
    cf_idx <- setdiff(fold1_idx, fit_idx)

    list(
        fit_idx = fit_idx,
        cf_idx = cf_idx,
        fit_treated_idx = fit_treated_idx,
        fit_control_idx = fit_control_idx,
        cf_treated_idx = setdiff(treated_idx, fit_treated_idx),
        cf_control_idx = setdiff(control_idx, fit_control_idx)
    )
}


stage1_weight_diagnostics <- function(weights_1_to_0, weights_0_to_1,
                                      weight_upper = 20,
                                      weight_lower = 0.05) {
    list(
        source_1_target_0 = weight_tail_summary(weights_1_to_0, weight_upper, weight_lower),
        source_0_target_1 = weight_tail_summary(weights_0_to_1, weight_upper, weight_lower)
    )
}


build_ite_unions_from_counterfactuals <- function(T_vec, Y_vec,
                                                  intervals_Y0,
                                                  intervals_Y1) {
    n <- length(T_vec)
    intervals_ITE <- vector("list", n)
    for (i in seq_len(n)) {
        if (T_vec[i] == 1) {
            intervals_cf <- normalize_interval_union(intervals_Y0[[i]])
            intervals_ITE[[i]] <- normalize_interval_union(cbind(
                Y_vec[i] - intervals_cf[, 2],
                Y_vec[i] - intervals_cf[, 1]
            ))
        } else {
            intervals_cf <- normalize_interval_union(intervals_Y1[[i]])
            intervals_ITE[[i]] <- normalize_interval_union(cbind(
                intervals_cf[, 1] - Y_vec[i],
                intervals_cf[, 2] - Y_vec[i]
            ))
        }
    }
    intervals_ITE
}


prediction_sets_to_dataframe <- function(intervals_list) {
    data.frame(
        mean_length = mean(sapply(intervals_list, union_total_length)),
        mean_n_intervals = mean(sapply(intervals_list, nrow))
    )
}

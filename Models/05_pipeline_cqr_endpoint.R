################################################################################
# Module: 05_pipeline_cqr_endpoint.R
# Purpose: Lei-style endpoint CQR baseline for induced ITE intervals
################################################################################

source("05_pipeline_utils.R")
source("03_outcome_cqr.R")
source("04_ite_cqr_endpoint.R")


run_cqr_endpoint <- function(X, T, Y,
                             fold1_idx,
                             fold2_idx,
                             test_idx,
                             scenario,
                             ps_mode = "estimated_ps",
                             alpha = 0.1,
                             train_prop_stage1 = 0.75,
                             train_prop_stage2 = 0.75,
                             n_post = 1000,
                             quantile_stage2 = 0.5,
                             weight_upper = 20,
                             weight_lower = 0.05,
                             true_ps = NULL,
                             verbose = TRUE,
                             seed = NULL,
                             ...) {
    alpha_stage <- alpha / 2
    q_lower <- alpha_stage / 2
    q_upper <- 1 - alpha_stage / 2
    seed_Y1 <- if (!is.null(seed)) seed * 1000 + 43 else NULL
    seed_Y0 <- if (!is.null(seed)) seed * 1000 + 44 else NULL
    seed_stage2 <- if (!is.null(seed)) seed * 1000 + 45 else NULL
    seed_stage1_split <- if (!is.null(seed)) seed * 1000 + 46 else NULL

    stage1_split <- split_fold1_for_stage1(
        fold1_idx, T, fit_prop = train_prop_stage1, seed = seed_stage1_split
    )

    ps_info <- resolve_propensity_scores(X, T, stage1_split$fit_idx, scenario, ps_mode, true_ps, ...)
    ps_all <- ps_info$ps_all
    weights_1_to_0 <- compute_conformal_weights(ps_all, 1, 0)
    weights_0_to_1 <- compute_conformal_weights(ps_all, 0, 1)

    if (verbose) cat("=== CQR_endpoint ===\n")
    fit_Y1 <- fit_outcome_cqr(
        X[stage1_split$fit_treated_idx, , drop = FALSE], Y[stage1_split$fit_treated_idx],
        X_cal = X[stage1_split$cf_treated_idx, , drop = FALSE],
        Y_cal = Y[stage1_split$cf_treated_idx],
        weights_cal = weights_1_to_0[stage1_split$cf_treated_idx],
        lower_quantile = q_lower,
        upper_quantile = q_upper,
        n_post = n_post,
        seed = seed_Y1
    )
    fit_Y0 <- fit_outcome_cqr(
        X[stage1_split$fit_control_idx, , drop = FALSE], Y[stage1_split$fit_control_idx],
        X_cal = X[stage1_split$cf_control_idx, , drop = FALSE],
        Y_cal = Y[stage1_split$cf_control_idx],
        weights_cal = weights_0_to_1[stage1_split$cf_control_idx],
        lower_quantile = q_lower,
        upper_quantile = q_upper,
        n_post = n_post,
        seed = seed_Y0
    )

    X_fold2 <- X[fold2_idx, , drop = FALSE]
    T_fold2 <- T[fold2_idx]
    Y_fold2 <- Y[fold2_idx]
    intervals_ITE_fold2 <- data.frame(lower = NA_real_, upper = NA_real_)[rep(1, length(fold2_idx)), ]

    treated_mask <- T_fold2 == 1
    if (any(treated_mask)) {
        intervals_Y0 <- predict_outcome_cqr(
            fit_Y0, X_fold2[treated_mask, , drop = FALSE],
            alpha = alpha_stage,
            weights_test = weights_0_to_1[fold2_idx][treated_mask],
            weight_upper = weight_upper,
            weight_lower = weight_lower
        )
        intervals_ITE_fold2$lower[treated_mask] <- Y_fold2[treated_mask] - intervals_Y0$upper
        intervals_ITE_fold2$upper[treated_mask] <- Y_fold2[treated_mask] - intervals_Y0$lower
    }
    control_mask <- T_fold2 == 0
    if (any(control_mask)) {
        intervals_Y1 <- predict_outcome_cqr(
            fit_Y1, X_fold2[control_mask, , drop = FALSE],
            alpha = alpha_stage,
            weights_test = weights_1_to_0[fold2_idx][control_mask],
            weight_upper = weight_upper,
            weight_lower = weight_lower
        )
        intervals_ITE_fold2$lower[control_mask] <- intervals_Y1$lower - Y_fold2[control_mask]
        intervals_ITE_fold2$upper[control_mask] <- intervals_Y1$upper - Y_fold2[control_mask]
    }

    fit_ITE <- fit_ite_cqr_endpoint(
        X_fold2,
        as.matrix(intervals_ITE_fold2),
        train_prop = train_prop_stage2,
        quantile = quantile_stage2,
        seed = seed_stage2
    )
    intervals_ITE_test <- predict_ite_cqr_endpoint(fit_ITE, X[test_idx, , drop = FALSE], alpha_stage)
    widths <- intervals_ITE_test$upper - intervals_ITE_test$lower

    structure(
        list(
            intervals_ITE_test = intervals_ITE_test,
            intervals_ITE_fold2 = intervals_ITE_fold2,
            fit_Y0 = fit_Y0,
            fit_Y1 = fit_Y1,
            fit_ITE = fit_ITE,
            ps_all = ps_all,
            weights_stage1 = list(source_1_target_0 = weights_1_to_0, source_0_target_1 = weights_0_to_1),
            diagnostics = list(
                stage1 = list(
                    ps_fit_source = if (ps_mode == "true_ps") "true_ps" else "I_fit",
                    ps_fit_n = length(stage1_split$fit_idx),
                    fit_n = length(stage1_split$fit_idx),
                    cal_n = length(stage1_split$cf_idx),
                    fit_n_treated = length(stage1_split$fit_treated_idx),
                    fit_n_control = length(stage1_split$fit_control_idx),
                    cal_n_treated = length(stage1_split$cf_treated_idx),
                    cal_n_control = length(stage1_split$cf_control_idx),
                    avg_width = mean(intervals_ITE_fold2$upper - intervals_ITE_fold2$lower),
                    avg_n_intervals = 1,
                    weight_diagnostics = stage1_weight_diagnostics(weights_1_to_0, weights_0_to_1)
                ),
                stage2_endpoint = list(
                    quantile = quantile_stage2,
                    score_mean = fit_ITE$score_mean
                ),
                test = list(avg_width = mean(widths), sd_width = stats::sd(widths), widths = widths)
            ),
            params = list(alpha = alpha, alpha_stage = alpha_stage, ps_mode = ps_mode)
        ),
        class = "cqr_endpoint_result"
    )
}

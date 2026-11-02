################################################################################
# Module: 05_pipeline_aps_classification.R
# Purpose: Main APS_classification method with multiple Stage 1 modes
################################################################################

source("05_pipeline_utils.R")
source("03_outcome_aps.R")
source("04_ite_aps_classification.R")


run_aps_classification <- function(X, T, Y,
                                   fold1_idx,
                                   fold2_idx,
                                   test_idx,
                                   scenario,
                                   stage1_mode = "estimated_ps",
                                   stage1_model = "dp",
                                   alpha = 0.1,
                                   train_prop_stage1 = 0.75,
                                   train_prop_stage2 = 0.75,
                                   stage1_n_grid = 500,
                                   stage1_grid_resolution = 0.05,
                                   stage1_max_grid_size = 2000,
                                   stage2_n_grid = 500,
                                   stage2_grid_resolution = 0.05,
                                   stage2_max_grid_size = 2000,
                                   n_samples = 100,
                                   pos_ratio = 0.5,
                                   strategy = "region_uniform",
                                   boundary_ratio = 0.5,
                                   stage2_method = "mixture",
                                   gamma_score = 0,
                                   dp_n_fit_max = 10000,
                                   dp_n_iter = 250,
                                   dp_burn = 125,
                                   dp_thin = 5,
                                   dp_alpha = 5,
                                   dp_kappa0 = 0.05,
                                   dp_nu0_extra = 3,
                                   dp_prior_scale = 1,
                                   dp_jitter = 1e-7,
                                   mdn_n_components = 5,
                                   mdn_hidden = c(64L, 64L),
                                   mdn_epochs = 500,
                                   mdn_batch_size = 256,
                                   mdn_lr = 1e-3,
                                   mdn_weight_decay = 1e-4,
                                   mdn_patience = 30,
                                   mdn_val_prop = 0.2,
                                   mdn_min_sigma = 0.02,
                                   mdn_grid_tail_prob = 5e-4,
                                   mdn_grid_sd_mult = 6,
                                   nsf_n_bins = 12,
                                   nsf_hidden = c(128L, 128L),
                                   nsf_epochs = 500,
                                   nsf_batch_size = 256,
                                   nsf_lr = 1e-3,
                                   nsf_weight_decay = 1e-4,
                                   nsf_patience = 30,
                                   nsf_val_prop = 0.2,
                                   nsf_tail_bound = 5,
                                   nsf_min_bin_width = 1e-3,
                                   nsf_min_bin_height = 1e-3,
                                   nsf_min_derivative = 1e-3,
                                   nsf_grid_tail_prob = 5e-4,
                                   nsf_grid_sd_mult = 6,
                                   weight_upper = 20,
                                   weight_lower = 0.05,
                                   true_ps = NULL,
                                   verbose = TRUE,
                                   seed = NULL,
                                   ...,
                                   drbart_nburn = 5000,
                                   drbart_nsim = 5000,
                                   drbart_nthin = 1,
                                   drbart_m_mean = 200,
                                   drbart_m_var = 100,
                                   drbart_variance = "ux",
                                   drbart_tree_alpha = 0.95,
                                   drbart_tree_beta = 2,
                                   drbart_lambda = 1,
                                   drbart_nu = 2,
                                   drbart_kfac = 2,
                                   drbart_phi0 = 1,
                                   drbart_grid_tail_prob = 5e-4,
                                   drbart_predict_batch_size = 100,
                                   mdn_device = "cpu") {
    if (!stage1_mode %in% c("true_ps", "estimated_ps", "misspecified_ps", "oracle")) {
        stop("stage1_mode must be true_ps, estimated_ps, misspecified_ps, or oracle.")
    }
    if (!stage1_model %in% c("dp", "bart", "mdn", "nsf", "drbart")) {
        stop("stage1_model must be dp, bart, mdn, nsf, or drbart.")
    }
    if (!stage2_method %in% c("mixture", "min")) {
        stop("stage2_method must be mixture or min.")
    }
    if (stage1_mode == "misspecified_ps" &&
        !scenario %in% c("S3", "S3_1", "S3_2")) {
        stop("misspecified_ps is only available for the S3 family.")
    }
    stage1_grid <- normalize_outcome_grid_params(
        stage1_n_grid,
        stage1_grid_resolution,
        stage1_max_grid_size
    )
    stage2_grid <- stage2_grid_params(
        stage2_n_grid,
        stage2_grid_resolution,
        stage2_max_grid_size
    )

    alpha_stage <- alpha / 2
    seed_Y1 <- if (!is.null(seed)) seed * 1000 + 43 else NULL
    seed_Y0 <- if (!is.null(seed)) seed * 1000 + 44 else NULL
    seed_stage2 <- if (!is.null(seed)) seed * 1000 + 45 else NULL
    seed_stage1_split <- if (!is.null(seed)) seed * 1000 + 46 else NULL

    stage1_split <- split_fold1_for_stage1(
        fold1_idx, T, fit_prop = train_prop_stage1, seed = seed_stage1_split
    )

    ps_mode <- if (stage1_mode == "oracle") "true_ps" else stage1_mode
    ps_info <- resolve_propensity_scores(X, T, stage1_split$fit_idx, scenario, ps_mode, true_ps, ...)
    ps_all <- ps_info$ps_all
    weights_1_to_0 <- compute_conformal_weights(ps_all, 1, 0)
    weights_0_to_1 <- compute_conformal_weights(ps_all, 0, 1)

    gamma_stage2 <- if (stage1_mode == "estimated_ps") gamma_score else 0
    if (verbose) {
        cat("=== APS_classification ===\n")
        cat(sprintf("Stage 1 mode: %s\n", stage1_mode))
        cat(sprintf("Stage 1 outcome model: %s\n", stage1_model))
        cat(sprintf("Stage 2 method: %s\n", stage2_method))
        cat(sprintf("Stage 2 gamma_score used: %.3f\n", gamma_stage2))
        cat(sprintf(
            "Stage 1 grid: n_grid=%d, grid_resolution=%.4f, max_grid_size=%d\n",
            stage1_grid$n_grid,
            stage1_grid$grid_resolution,
            stage1_grid$max_grid_size
        ))
        cat(sprintf(
            "Stage 2 grid: n_grid=%d, grid_resolution=%.4f, max_grid_size=%d\n",
            stage2_grid$n_grid,
            stage2_grid$grid_resolution,
            stage2_grid$max_grid_size
        ))
    }

    if (stage1_mode == "oracle") {
        fit_Y1 <- fit_outcome_aps_oracle(
            X[stage1_split$fit_treated_idx, , drop = FALSE], Y[stage1_split$fit_treated_idx],
            scenario = scenario,
            arm = 1,
            X_cal = X[stage1_split$cf_treated_idx, , drop = FALSE],
            Y_cal = Y[stage1_split$cf_treated_idx],
            weights_cal = weights_1_to_0[stage1_split$cf_treated_idx],
            n_grid = stage1_grid$n_grid,
            grid_resolution = stage1_grid$grid_resolution,
            max_grid_size = stage1_grid$max_grid_size,
            seed = seed_Y1,
            ...
        )
        fit_Y0 <- fit_outcome_aps_oracle(
            X[stage1_split$fit_control_idx, , drop = FALSE], Y[stage1_split$fit_control_idx],
            scenario = scenario,
            arm = 0,
            X_cal = X[stage1_split$cf_control_idx, , drop = FALSE],
            Y_cal = Y[stage1_split$cf_control_idx],
            weights_cal = weights_0_to_1[stage1_split$cf_control_idx],
            n_grid = stage1_grid$n_grid,
            grid_resolution = stage1_grid$grid_resolution,
            max_grid_size = stage1_grid$max_grid_size,
            seed = seed_Y0,
            ...
        )
        predict_stage1 <- predict_outcome_aps_oracle
    } else {
        fit_Y1 <- fit_outcome_aps(
            X[stage1_split$fit_treated_idx, , drop = FALSE], Y[stage1_split$fit_treated_idx],
            X_cal = X[stage1_split$cf_treated_idx, , drop = FALSE],
            Y_cal = Y[stage1_split$cf_treated_idx],
            weights_cal = weights_1_to_0[stage1_split$cf_treated_idx],
            stage1_model = stage1_model,
            n_grid = stage1_grid$n_grid,
            grid_resolution = stage1_grid$grid_resolution,
            max_grid_size = stage1_grid$max_grid_size,
            dp_n_fit_max = dp_n_fit_max,
            dp_n_iter = dp_n_iter,
            dp_burn = dp_burn,
            dp_thin = dp_thin,
            dp_alpha = dp_alpha,
            dp_kappa0 = dp_kappa0,
            dp_nu0_extra = dp_nu0_extra,
            dp_prior_scale = dp_prior_scale,
            dp_jitter = dp_jitter,
            mdn_n_components = mdn_n_components,
            mdn_hidden = mdn_hidden,
            mdn_epochs = mdn_epochs,
            mdn_batch_size = mdn_batch_size,
            mdn_lr = mdn_lr,
            mdn_weight_decay = mdn_weight_decay,
            mdn_patience = mdn_patience,
            mdn_val_prop = mdn_val_prop,
            mdn_min_sigma = mdn_min_sigma,
            mdn_grid_tail_prob = mdn_grid_tail_prob,
            mdn_grid_sd_mult = mdn_grid_sd_mult,
            mdn_device = mdn_device,
            nsf_n_bins = nsf_n_bins,
            nsf_hidden = nsf_hidden,
            nsf_epochs = nsf_epochs,
            nsf_batch_size = nsf_batch_size,
            nsf_lr = nsf_lr,
            nsf_weight_decay = nsf_weight_decay,
            nsf_patience = nsf_patience,
            nsf_val_prop = nsf_val_prop,
            nsf_tail_bound = nsf_tail_bound,
            nsf_min_bin_width = nsf_min_bin_width,
            nsf_min_bin_height = nsf_min_bin_height,
            nsf_min_derivative = nsf_min_derivative,
            nsf_grid_tail_prob = nsf_grid_tail_prob,
            nsf_grid_sd_mult = nsf_grid_sd_mult,
            drbart_nburn = drbart_nburn,
            drbart_nsim = drbart_nsim,
            drbart_nthin = drbart_nthin,
            drbart_m_mean = drbart_m_mean,
            drbart_m_var = drbart_m_var,
            drbart_variance = drbart_variance,
            drbart_tree_alpha = drbart_tree_alpha,
            drbart_tree_beta = drbart_tree_beta,
            drbart_lambda = drbart_lambda,
            drbart_nu = drbart_nu,
            drbart_kfac = drbart_kfac,
            drbart_phi0 = drbart_phi0,
            drbart_grid_tail_prob = drbart_grid_tail_prob,
            drbart_predict_batch_size = drbart_predict_batch_size,
            seed = seed_Y1,
            verbose = verbose
        )
        fit_Y0 <- fit_outcome_aps(
            X[stage1_split$fit_control_idx, , drop = FALSE], Y[stage1_split$fit_control_idx],
            X_cal = X[stage1_split$cf_control_idx, , drop = FALSE],
            Y_cal = Y[stage1_split$cf_control_idx],
            weights_cal = weights_0_to_1[stage1_split$cf_control_idx],
            stage1_model = stage1_model,
            n_grid = stage1_grid$n_grid,
            grid_resolution = stage1_grid$grid_resolution,
            max_grid_size = stage1_grid$max_grid_size,
            dp_n_fit_max = dp_n_fit_max,
            dp_n_iter = dp_n_iter,
            dp_burn = dp_burn,
            dp_thin = dp_thin,
            dp_alpha = dp_alpha,
            dp_kappa0 = dp_kappa0,
            dp_nu0_extra = dp_nu0_extra,
            dp_prior_scale = dp_prior_scale,
            dp_jitter = dp_jitter,
            mdn_n_components = mdn_n_components,
            mdn_hidden = mdn_hidden,
            mdn_epochs = mdn_epochs,
            mdn_batch_size = mdn_batch_size,
            mdn_lr = mdn_lr,
            mdn_weight_decay = mdn_weight_decay,
            mdn_patience = mdn_patience,
            mdn_val_prop = mdn_val_prop,
            mdn_min_sigma = mdn_min_sigma,
            mdn_grid_tail_prob = mdn_grid_tail_prob,
            mdn_grid_sd_mult = mdn_grid_sd_mult,
            mdn_device = mdn_device,
            nsf_n_bins = nsf_n_bins,
            nsf_hidden = nsf_hidden,
            nsf_epochs = nsf_epochs,
            nsf_batch_size = nsf_batch_size,
            nsf_lr = nsf_lr,
            nsf_weight_decay = nsf_weight_decay,
            nsf_patience = nsf_patience,
            nsf_val_prop = nsf_val_prop,
            nsf_tail_bound = nsf_tail_bound,
            nsf_min_bin_width = nsf_min_bin_width,
            nsf_min_bin_height = nsf_min_bin_height,
            nsf_min_derivative = nsf_min_derivative,
            nsf_grid_tail_prob = nsf_grid_tail_prob,
            nsf_grid_sd_mult = nsf_grid_sd_mult,
            drbart_nburn = drbart_nburn,
            drbart_nsim = drbart_nsim,
            drbart_nthin = drbart_nthin,
            drbart_m_mean = drbart_m_mean,
            drbart_m_var = drbart_m_var,
            drbart_variance = drbart_variance,
            drbart_tree_alpha = drbart_tree_alpha,
            drbart_tree_beta = drbart_tree_beta,
            drbart_lambda = drbart_lambda,
            drbart_nu = drbart_nu,
            drbart_kfac = drbart_kfac,
            drbart_phi0 = drbart_phi0,
            drbart_grid_tail_prob = drbart_grid_tail_prob,
            drbart_predict_batch_size = drbart_predict_batch_size,
            seed = seed_Y0,
            verbose = verbose
        )
        predict_stage1 <- predict_outcome_aps
    }

    X_fold2 <- X[fold2_idx, , drop = FALSE]
    T_fold2 <- T[fold2_idx]
    Y_fold2 <- Y[fold2_idx]
    n_fold2 <- length(fold2_idx)
    intervals_Y0 <- vector("list", n_fold2)
    intervals_Y1 <- vector("list", n_fold2)

    treated_mask <- T_fold2 == 1
    if (any(treated_mask)) {
        intervals_Y0[treated_mask] <- predict_stage1(
            fit_Y0,
            X_fold2[treated_mask, , drop = FALSE],
            alpha = alpha_stage,
            weights_test = weights_0_to_1[fold2_idx][treated_mask],
            weight_upper = weight_upper,
            weight_lower = weight_lower,
            n_grid = stage1_grid$n_grid,
            grid_resolution = stage1_grid$grid_resolution,
            max_grid_size = stage1_grid$max_grid_size,
            ...
        )
    }

    control_mask <- T_fold2 == 0
    if (any(control_mask)) {
        intervals_Y1[control_mask] <- predict_stage1(
            fit_Y1,
            X_fold2[control_mask, , drop = FALSE],
            alpha = alpha_stage,
            weights_test = weights_1_to_0[fold2_idx][control_mask],
            weight_upper = weight_upper,
            weight_lower = weight_lower,
            n_grid = stage1_grid$n_grid,
            grid_resolution = stage1_grid$grid_resolution,
            max_grid_size = stage1_grid$max_grid_size,
            ...
        )
    }

    intervals_ITE_fold2 <- build_ite_unions_from_counterfactuals(
        T_fold2, Y_fold2, intervals_Y0, intervals_Y1
    )
    stage1_widths <- sapply(intervals_ITE_fold2, union_total_length)
    stage1_n_intervals <- sapply(intervals_ITE_fold2, nrow)

    if (stage2_method == "mixture") {
        fit_ITE <- fit_ite_aps_classification(
            X_fold2,
            intervals_ITE_fold2,
            train_prop = train_prop_stage2,
            n_samples = n_samples,
            pos_ratio = pos_ratio,
            strategy = strategy,
            boundary_ratio = boundary_ratio,
            gamma_score = gamma_stage2,
            n_grid = stage2_grid$n_grid,
            grid_resolution = stage2_grid$grid_resolution,
            max_grid_size = stage2_grid$max_grid_size,
            verbose = verbose,
            seed = seed_stage2
        )
        intervals_ITE_test <- predict_ite_aps_classification(
            fit_ITE,
            X[test_idx, , drop = FALSE],
            alpha = alpha_stage
        )
    } else {
        fit_ITE <- fit_ite_aps_classification_min(
            X_fold2,
            T_fold2,
            intervals_ITE_fold2,
            train_prop = train_prop_stage2,
            n_samples = n_samples,
            pos_ratio = pos_ratio,
            strategy = strategy,
            boundary_ratio = boundary_ratio,
            gamma_score = gamma_stage2,
            n_grid = stage2_grid$n_grid,
            grid_resolution = stage2_grid$grid_resolution,
            max_grid_size = stage2_grid$max_grid_size,
            verbose = verbose,
            seed = seed_stage2
        )
        intervals_ITE_test <- predict_ite_aps_classification_min(
            fit_ITE,
            X[test_idx, , drop = FALSE],
            alpha = alpha_stage
        )
    }
    widths_test <- sapply(intervals_ITE_test, union_total_length)
    n_intervals_test <- sapply(intervals_ITE_test, nrow)

    dp_diag <- function(fit) {
        if (!is.null(fit$dp_diagnostics)) {
            fit$dp_diagnostics
        } else {
            list(
                k_mean_last_half = NA_real_,
                k_min_last_half = NA_real_,
                k_max_last_half = NA_real_,
                posterior_draws = NA_real_,
                n_fit_used = NA_real_
            )
        }
    }
    diag_Y1 <- dp_diag(fit_Y1)
    diag_Y0 <- dp_diag(fit_Y0)

    mdn_diag <- function(fit) {
        if (!is.null(fit$mdn_diagnostics)) {
            fit$mdn_diagnostics
        } else {
            list(
                n_fit_used = NA_real_,
                n_components = NA_real_,
                hidden = NA_character_,
                epochs_trained = NA_real_,
                best_val_loss = NA_real_,
                final_train_loss = NA_real_,
                final_val_loss = NA_real_
            )
        }
    }
    mdn_diag_Y1 <- mdn_diag(fit_Y1)
    mdn_diag_Y0 <- mdn_diag(fit_Y0)

    nsf_diag <- function(fit) {
        if (!is.null(fit$nsf_diagnostics)) {
            fit$nsf_diagnostics
        } else {
            list(
                n_fit_used = NA_real_,
                n_bins = NA_real_,
                hidden = NA_character_,
                tail_bound = NA_real_,
                epochs_trained = NA_real_,
                best_val_loss = NA_real_,
                final_train_loss = NA_real_,
                final_val_loss = NA_real_
            )
        }
    }
    nsf_diag_Y1 <- nsf_diag(fit_Y1)
    nsf_diag_Y0 <- nsf_diag(fit_Y0)

    weights_val_Y1 <- if (!is.null(fit_Y1$weights_val)) fit_Y1$weights_val else numeric(0)
    weights_val_Y0 <- if (!is.null(fit_Y0$weights_val)) fit_Y0$weights_val else numeric(0)

    result <- structure(
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
                    mode = stage1_mode,
                    model = if (stage1_mode == "oracle") "oracle" else stage1_model,
                    ps_fit_source = if (ps_mode == "true_ps") "true_ps" else "I_fit",
                    ps_fit_n = length(stage1_split$fit_idx),
                    fit_n = length(stage1_split$fit_idx),
                    cal_n = length(stage1_split$cf_idx),
                    fit_n_treated = length(stage1_split$fit_treated_idx),
                    fit_n_control = length(stage1_split$fit_control_idx),
                    cal_n_treated = length(stage1_split$cf_treated_idx),
                    cal_n_control = length(stage1_split$cf_control_idx),
                    avg_width = mean(stage1_widths),
                    sd_width = stats::sd(stage1_widths),
                    avg_n_intervals = mean(stage1_n_intervals),
                    stage1_widths = stage1_widths,
                    stage1_n_intervals = stage1_n_intervals,
                    weight_diagnostics = stage1_weight_diagnostics(weights_1_to_0, weights_0_to_1),
                    calibration_ess_Y1 = effective_sample_size(weights_val_Y1),
                    calibration_ess_Y0 = effective_sample_size(weights_val_Y0),
                    dp_k_mean_Y1 = diag_Y1$k_mean_last_half,
                    dp_k_mean_Y0 = diag_Y0$k_mean_last_half,
                    dp_k_min_Y1 = diag_Y1$k_min_last_half,
                    dp_k_min_Y0 = diag_Y0$k_min_last_half,
                    dp_k_max_Y1 = diag_Y1$k_max_last_half,
                    dp_k_max_Y0 = diag_Y0$k_max_last_half,
                    mdn_n_components_Y1 = mdn_diag_Y1$n_components,
                    mdn_n_components_Y0 = mdn_diag_Y0$n_components,
                    mdn_epochs_trained_Y1 = mdn_diag_Y1$epochs_trained,
                    mdn_epochs_trained_Y0 = mdn_diag_Y0$epochs_trained,
                    mdn_best_val_loss_Y1 = mdn_diag_Y1$best_val_loss,
                    mdn_best_val_loss_Y0 = mdn_diag_Y0$best_val_loss,
                    mdn_final_train_loss_Y1 = mdn_diag_Y1$final_train_loss,
                    mdn_final_train_loss_Y0 = mdn_diag_Y0$final_train_loss,
                    mdn_final_val_loss_Y1 = mdn_diag_Y1$final_val_loss,
                    mdn_final_val_loss_Y0 = mdn_diag_Y0$final_val_loss,
                    nsf_n_bins_Y1 = nsf_diag_Y1$n_bins,
                    nsf_n_bins_Y0 = nsf_diag_Y0$n_bins,
                    nsf_epochs_trained_Y1 = nsf_diag_Y1$epochs_trained,
                    nsf_epochs_trained_Y0 = nsf_diag_Y0$epochs_trained,
                    nsf_best_val_loss_Y1 = nsf_diag_Y1$best_val_loss,
                    nsf_best_val_loss_Y0 = nsf_diag_Y0$best_val_loss,
                    nsf_final_train_loss_Y1 = nsf_diag_Y1$final_train_loss,
                    nsf_final_train_loss_Y0 = nsf_diag_Y0$final_train_loss,
                    nsf_final_val_loss_Y1 = nsf_diag_Y1$final_val_loss,
                    nsf_final_val_loss_Y0 = nsf_diag_Y0$final_val_loss,
                    nsf_tail_bound_Y1 = nsf_diag_Y1$tail_bound,
                    nsf_tail_bound_Y0 = nsf_diag_Y0$tail_bound,
                    grid_params_Y1 = fit_Y1$grid_params,
                    grid_params_Y0 = fit_Y0$grid_params
                ),
                stage2_classification = list(
                    method = stage2_method,
                    n_samples = n_samples,
                    pos_ratio = pos_ratio,
                    strategy = strategy,
                    boundary_ratio = boundary_ratio,
                    gamma_score = gamma_stage2,
                    requested_gamma_score = gamma_score,
                    region = fit_ITE$region,
                    grid_params = fit_ITE$grid_params,
                    grid_diagnostics = fit_ITE$grid_diagnostics,
                    score_mean = fit_ITE$score_mean
                ),
                test = list(
                    avg_width = mean(widths_test),
                    sd_width = stats::sd(widths_test),
                    avg_n_intervals = mean(n_intervals_test),
                    sd_n_intervals = stats::sd(n_intervals_test),
                    widths = widths_test,
                    n_intervals = n_intervals_test
                )
            ),
            params = list(
                alpha = alpha,
                alpha_stage = alpha_stage,
                stage1_mode = stage1_mode,
                stage1_model = if (stage1_mode == "oracle") "oracle" else stage1_model,
                ps_mode = ps_mode,
                stage2_method = stage2_method,
                gamma_score = gamma_stage2,
                requested_gamma_score = gamma_score,
                train_prop_stage1 = train_prop_stage1,
                train_prop_stage2 = train_prop_stage2,
                stage1_n_grid = stage1_grid$n_grid,
                stage1_grid_resolution = stage1_grid$grid_resolution,
                stage1_max_grid_size = stage1_grid$max_grid_size,
                stage2_n_grid = stage2_grid$n_grid,
                stage2_grid_resolution = stage2_grid$grid_resolution,
                stage2_max_grid_size = stage2_grid$max_grid_size,
                n_samples = n_samples,
                pos_ratio = pos_ratio,
                strategy = strategy,
                boundary_ratio = boundary_ratio,
                mdn_n_components = mdn_n_components,
                mdn_hidden = paste(normalize_mdn_hidden(mdn_hidden), collapse = ","),
                mdn_epochs = mdn_epochs,
                mdn_batch_size = mdn_batch_size,
                mdn_lr = mdn_lr,
                mdn_weight_decay = mdn_weight_decay,
                mdn_patience = mdn_patience,
                mdn_val_prop = mdn_val_prop,
                mdn_min_sigma = mdn_min_sigma,
                mdn_grid_tail_prob = mdn_grid_tail_prob,
                mdn_grid_sd_mult = mdn_grid_sd_mult,
                nsf_n_bins = nsf_n_bins,
                nsf_hidden = paste(normalize_nsf_hidden(nsf_hidden), collapse = ","),
                nsf_epochs = nsf_epochs,
                nsf_batch_size = nsf_batch_size,
                nsf_lr = nsf_lr,
                nsf_weight_decay = nsf_weight_decay,
                nsf_patience = nsf_patience,
                nsf_val_prop = nsf_val_prop,
                nsf_tail_bound = nsf_tail_bound,
                nsf_min_bin_width = nsf_min_bin_width,
                nsf_min_bin_height = nsf_min_bin_height,
                nsf_min_derivative = nsf_min_derivative,
                nsf_grid_tail_prob = nsf_grid_tail_prob,
                nsf_grid_sd_mult = nsf_grid_sd_mult
            )
        ),
        class = "aps_classification_result"
    )
    if (stage1_mode != "oracle" && stage1_model == "mdn") {
        result$params$mdn_device <- mdn_device
        result$params$mdn <- c(fit_Y1$density_model$params, list(
            device_requested = mdn_device,
            torch_version = as.character(utils::packageVersion("torch"))
        ))
        result$diagnostics$stage1$mdn_device_Y1 <- fit_Y1$density_model$device
        result$diagnostics$stage1$mdn_device_Y0 <- fit_Y0$density_model$device
    }
    if (stage1_mode != "oracle" && stage1_model == "drbart") {
        result$params$drbart <- fit_Y1$density_model$params
        result$diagnostics$stage1$drbart_Y1 <- fit_Y1$density_model$diagnostics
        result$diagnostics$stage1$drbart_Y0 <- fit_Y0$density_model$diagnostics
    }
    result
}

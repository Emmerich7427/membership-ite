#!/usr/bin/env Rscript
################################################################################
# Script: 02_run_aps_classification.R
# Purpose: Run APS_classification experiments
################################################################################

get_script_dir <- function() {
    file_arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
    normalizePath(getwd())
}

SCRIPT_DIR <- get_script_dir()
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))

parse_int_csv <- function(x) {
    pieces <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
    pieces <- pieces[nzchar(pieces)]
    if (length(pieces) == 0L) {
        return(integer(0))
    }
    out <- as.integer(pieces)
    if (any(is.na(out)) || any(out < 1L)) {
        stop("Expected a comma-separated list of positive integers.")
    }
    out
}

parse_args <- function(args) {
    config <- list(
        scenario = "S3", rho = 0, alpha = 0.1, n_simulations = 1,
        run_single_seed = FALSE, target_seed = NULL,
        stage1_mode = "estimated_ps", stage1_model = "dp", gamma_score = 0,
        separate_stage1_model_dir = FALSE,
        stage2_method = "mixture",
        stage1_n_grid = 500, stage1_grid_resolution = 0.05, stage1_max_grid_size = 2000,
        stage2_n_grid = 1000, stage2_grid_resolution = 0.02, stage2_max_grid_size = 5000,
        n_samples = 100, pos_ratio = 0.5, strategy = "region_uniform", boundary_ratio = 0.5,
        dp_n_fit_max = 10000, dp_n_iter = 250, dp_burn = 125, dp_thin = 5,
        dp_alpha = 5, dp_kappa0 = 0.05, dp_nu0_extra = 3, dp_prior_scale = 1,
        mdn_n_components = 5, mdn_hidden = c(64L, 64L),
        mdn_epochs = 500, mdn_batch_size = 256,
        mdn_lr = 1e-3, mdn_weight_decay = 1e-4, mdn_patience = 30,
        mdn_val_prop = 0.2, mdn_min_sigma = 0.02,
        mdn_grid_tail_prob = 5e-4, mdn_grid_sd_mult = 6,
        mdn_device = "cpu",
        nsf_n_bins = 12, nsf_hidden = c(128L, 128L),
        nsf_epochs = 500, nsf_batch_size = 256,
        nsf_lr = 1e-3, nsf_weight_decay = 1e-4, nsf_patience = 30,
        nsf_val_prop = 0.2, nsf_tail_bound = 5,
        nsf_min_bin_width = 1e-3, nsf_min_bin_height = 1e-3,
        nsf_min_derivative = 1e-3,
        nsf_grid_tail_prob = 5e-4, nsf_grid_sd_mult = 6,
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
        overlap_strength = 1.0
    )
    i <- 1
    while (i <= length(args)) {
        arg <- args[i]
        if (arg == "--scenario" && i < length(args)) {
            config$scenario <- args[i + 1]; i <- i + 2
        } else if (arg == "--output-dir" && i < length(args)) {
            config$output_dir <- normalizePath(args[i + 1], winslash = "/", mustWork = FALSE)
            i <- i + 2
        } else if (arg == "--rho" && i < length(args)) {
            config$rho <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--alpha" && i < length(args)) {
            config$alpha <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--n-sim" && i < length(args)) {
            config$n_simulations <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--seed" && i < length(args)) {
            config$run_single_seed <- TRUE; config$target_seed <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--stage1-mode" && i < length(args)) {
            config$stage1_mode <- args[i + 1]; i <- i + 2
        } else if (arg == "--stage1-model" && i < length(args)) {
            config$stage1_model <- args[i + 1]; i <- i + 2
        } else if (arg == "--separate-stage1-model-dir") {
            config$separate_stage1_model_dir <- TRUE; i <- i + 1
        } else if (arg == "--stage2-method" && i < length(args)) {
            config$stage2_method <- args[i + 1]; i <- i + 2
        } else if (arg == "--gamma_score" && i < length(args)) {
            config$gamma_score <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--stage1_n_grid" && i < length(args)) {
            config$stage1_n_grid <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--stage1_grid_resolution" && i < length(args)) {
            config$stage1_grid_resolution <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--stage1_max_grid_size" && i < length(args)) {
            config$stage1_max_grid_size <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--stage2_n_grid" && i < length(args)) {
            config$stage2_n_grid <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--stage2_grid_resolution" && i < length(args)) {
            config$stage2_grid_resolution <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--stage2_max_grid_size" && i < length(args)) {
            config$stage2_max_grid_size <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--n_samples" && i < length(args)) {
            config$n_samples <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--strategy" && i < length(args)) {
            config$strategy <- args[i + 1]; i <- i + 2
        } else if (arg == "--overlap-strength" && i < length(args)) {
            config$overlap_strength <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--dp_n_iter" && i < length(args)) {
            config$dp_n_iter <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--dp_burn" && i < length(args)) {
            config$dp_burn <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--dp_thin" && i < length(args)) {
            config$dp_thin <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--mdn_device" && i < length(args)) {
            config$mdn_device <- args[i + 1]
            i <- i + 2
        } else if (arg == "--mdn_n_components" && i < length(args)) {
            config$mdn_n_components <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--mdn_hidden" && i < length(args)) {
            config$mdn_hidden <- parse_int_csv(args[i + 1]); i <- i + 2
        } else if (arg == "--mdn_epochs" && i < length(args)) {
            config$mdn_epochs <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--mdn_batch_size" && i < length(args)) {
            config$mdn_batch_size <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--mdn_lr" && i < length(args)) {
            config$mdn_lr <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--mdn_weight_decay" && i < length(args)) {
            config$mdn_weight_decay <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--mdn_patience" && i < length(args)) {
            config$mdn_patience <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--mdn_val_prop" && i < length(args)) {
            config$mdn_val_prop <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--mdn_min_sigma" && i < length(args)) {
            config$mdn_min_sigma <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--mdn_grid_tail_prob" && i < length(args)) {
            config$mdn_grid_tail_prob <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--mdn_grid_sd_mult" && i < length(args)) {
            config$mdn_grid_sd_mult <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--nsf_n_bins" && i < length(args)) {
            config$nsf_n_bins <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--nsf_hidden" && i < length(args)) {
            config$nsf_hidden <- parse_int_csv(args[i + 1]); i <- i + 2
        } else if (arg == "--nsf_epochs" && i < length(args)) {
            config$nsf_epochs <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--nsf_batch_size" && i < length(args)) {
            config$nsf_batch_size <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--nsf_lr" && i < length(args)) {
            config$nsf_lr <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--nsf_weight_decay" && i < length(args)) {
            config$nsf_weight_decay <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--nsf_patience" && i < length(args)) {
            config$nsf_patience <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--nsf_val_prop" && i < length(args)) {
            config$nsf_val_prop <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--nsf_tail_bound" && i < length(args)) {
            config$nsf_tail_bound <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--nsf_min_bin_width" && i < length(args)) {
            config$nsf_min_bin_width <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--nsf_min_bin_height" && i < length(args)) {
            config$nsf_min_bin_height <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--nsf_min_derivative" && i < length(args)) {
            config$nsf_min_derivative <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--nsf_grid_tail_prob" && i < length(args)) {
            config$nsf_grid_tail_prob <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--drbart_nburn" && i < length(args)) {
            config$drbart_nburn <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--drbart_nsim" && i < length(args)) {
            config$drbart_nsim <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--drbart_nthin" && i < length(args)) {
            config$drbart_nthin <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--drbart_m_mean" && i < length(args)) {
            config$drbart_m_mean <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--drbart_m_var" && i < length(args)) {
            config$drbart_m_var <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--drbart_variance" && i < length(args)) {
            config$drbart_variance <- args[i + 1]
            i <- i + 2
        } else if (arg == "--drbart_tree_alpha" && i < length(args)) {
            config$drbart_tree_alpha <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--drbart_tree_beta" && i < length(args)) {
            config$drbart_tree_beta <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--drbart_lambda" && i < length(args)) {
            config$drbart_lambda <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--drbart_nu" && i < length(args)) {
            config$drbart_nu <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--drbart_kfac" && i < length(args)) {
            config$drbart_kfac <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--drbart_phi0" && i < length(args)) {
            config$drbart_phi0 <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--drbart_grid_tail_prob" && i < length(args)) {
            config$drbart_grid_tail_prob <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--drbart_predict_batch_size" && i < length(args)) {
            config$drbart_predict_batch_size <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--nsf_grid_sd_mult" && i < length(args)) {
            config$nsf_grid_sd_mult <- as.numeric(args[i + 1]); i <- i + 2
        } else {
            stop(sprintf("Unknown or incomplete argument: %s", arg))
        }
    }
    if (config$stage1_model == "mdn" &&
        !config$mdn_device %in% c("cpu", "cuda", "auto")) {
        stop("--mdn_device must be cpu, cuda, or auto.")
    }
    config
}

format_rho_label <- function(rho) paste0("rho_", format(rho, trim = TRUE, scientific = FALSE))
format_gamma_label <- function(gamma_score) paste0("gamma_", format(gamma_score, trim = TRUE, scientific = FALSE))
format_method_label <- function(stage2_method) {
    if (stage2_method == "mixture") return("APS_classification")
    if (stage2_method == "min") return("APS_classification_min")
    stop("stage2_method must be mixture or min.")
}

build_output_dir <- function(config) {
    if (!is.null(config$output_dir) && nzchar(config$output_dir)) return(config$output_dir)
    method_label <- format_method_label(config$stage2_method)
    parts <- c(PROJECT_ROOT, "Results", config$scenario, format_rho_label(config$rho),
               method_label, config$stage1_mode)
    use_model_dir <- config$stage1_mode != "oracle" &&
        (isTRUE(config$separate_stage1_model_dir) || config$stage1_model %in% c("mdn", "nsf", "drbart"))
    if (use_model_dir) {
        parts <- c(parts, config$stage1_model)
    }
    if (config$stage1_mode == "estimated_ps") {
        parts <- c(parts, format_gamma_label(config$gamma_score))
    }
    do.call(file.path, as.list(parts))
}

load_data_simulation <- function(path) {
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    if (!exists("data_simul", envir = env, inherits = FALSE)) {
        stop(sprintf("data_simul not found in %s", path))
    }
    get("data_simul", envir = env)
}

stage1_test_diagnostics_aps <- function(results, data_simul, config, n_train) {
    alpha_stage <- config$alpha / 2
    X_test <- data_simul$Xtest
    T_test <- data_simul$Ttest
    predictor <- if (config$stage1_mode == "oracle") predict_outcome_aps_oracle else predict_outcome_aps

    out <- list(Y1_cov = NA_real_, Y1_len = NA_real_, Y0_cov = NA_real_, Y0_len = NA_real_)
    idx_T0 <- which(T_test == 0)
    if (length(idx_T0) > 0) {
        sets_Y1 <- predictor(
            results$fit_Y1, X_test[idx_T0, , drop = FALSE], alpha = alpha_stage,
            weights_test = results$weights_stage1$source_1_target_0[n_train + idx_T0]
        )
        cov <- mapply(point_in_union, data_simul$Y1test[idx_T0], sets_Y1)
        out$Y1_cov <- mean(cov)
        out$Y1_len <- mean(sapply(sets_Y1, union_total_length))
    }

    idx_T1 <- which(T_test == 1)
    if (length(idx_T1) > 0) {
        sets_Y0 <- predictor(
            results$fit_Y0, X_test[idx_T1, , drop = FALSE], alpha = alpha_stage,
            weights_test = results$weights_stage1$source_0_target_1[n_train + idx_T1]
        )
        cov <- mapply(point_in_union, data_simul$Y0test[idx_T1], sets_Y0)
        out$Y0_cov <- mean(cov)
        out$Y0_len <- mean(sapply(sets_Y0, union_total_length))
    }
    out
}

run_one <- function(seed, config, data_dir) {
    data_simul <- load_data_simulation(file.path(data_dir, sprintf("synthetic_%03d.RData", seed)))
    n_train <- nrow(data_simul$Xtrain)
    n_test <- nrow(data_simul$Xtest)

    folds <- split_data_by_treatment(
        data_simul$Xtrain, data_simul$Ytrain, data_simul$Ttrain,
        fold_prop = 0.5, seed = seed * 1000 + 42
    )
    X_all <- rbind(data_simul$Xtrain, data_simul$Xtest)
    T_all <- c(data_simul$Ttrain, data_simul$Ttest)
    Y_all <- c(data_simul$Ytrain, data_simul$Ytest)
    true_ps <- c(data_simul$pstrain, data_simul$pstest)

    results <- run_aps_classification(
        X_all, T_all, Y_all,
        fold1_idx = folds$fold1_idx,
        fold2_idx = folds$fold2_idx,
        test_idx = (n_train + 1):(n_train + n_test),
        scenario = config$scenario,
        stage1_mode = config$stage1_mode,
        stage1_model = config$stage1_model,
        alpha = config$alpha,
        stage1_n_grid = config$stage1_n_grid,
        stage1_grid_resolution = config$stage1_grid_resolution,
        stage1_max_grid_size = config$stage1_max_grid_size,
        stage2_n_grid = config$stage2_n_grid,
        stage2_grid_resolution = config$stage2_grid_resolution,
        stage2_max_grid_size = config$stage2_max_grid_size,
        n_samples = config$n_samples,
        pos_ratio = config$pos_ratio,
        strategy = config$strategy,
        boundary_ratio = config$boundary_ratio,
        stage2_method = config$stage2_method,
        gamma_score = config$gamma_score,
        dp_n_fit_max = config$dp_n_fit_max,
        dp_n_iter = config$dp_n_iter,
        dp_burn = config$dp_burn,
        dp_thin = config$dp_thin,
        dp_alpha = config$dp_alpha,
        dp_kappa0 = config$dp_kappa0,
        dp_nu0_extra = config$dp_nu0_extra,
        dp_prior_scale = config$dp_prior_scale,
        mdn_n_components = config$mdn_n_components,
        mdn_hidden = config$mdn_hidden,
        mdn_epochs = config$mdn_epochs,
        mdn_batch_size = config$mdn_batch_size,
        mdn_lr = config$mdn_lr,
        mdn_weight_decay = config$mdn_weight_decay,
        mdn_patience = config$mdn_patience,
        mdn_val_prop = config$mdn_val_prop,
        mdn_min_sigma = config$mdn_min_sigma,
        mdn_grid_tail_prob = config$mdn_grid_tail_prob,
        mdn_grid_sd_mult = config$mdn_grid_sd_mult,
        mdn_device = config$mdn_device,
        nsf_n_bins = config$nsf_n_bins,
        nsf_hidden = config$nsf_hidden,
        nsf_epochs = config$nsf_epochs,
        nsf_batch_size = config$nsf_batch_size,
        nsf_lr = config$nsf_lr,
        nsf_weight_decay = config$nsf_weight_decay,
        nsf_patience = config$nsf_patience,
        nsf_val_prop = config$nsf_val_prop,
        nsf_tail_bound = config$nsf_tail_bound,
        nsf_min_bin_width = config$nsf_min_bin_width,
        nsf_min_bin_height = config$nsf_min_bin_height,
        nsf_min_derivative = config$nsf_min_derivative,
        nsf_grid_tail_prob = config$nsf_grid_tail_prob,
        nsf_grid_sd_mult = config$nsf_grid_sd_mult,
        drbart_nburn = config$drbart_nburn,
        drbart_nsim = config$drbart_nsim,
        drbart_nthin = config$drbart_nthin,
        drbart_m_mean = config$drbart_m_mean,
        drbart_m_var = config$drbart_m_var,
        drbart_variance = config$drbart_variance,
        drbart_tree_alpha = config$drbart_tree_alpha,
        drbart_tree_beta = config$drbart_tree_beta,
        drbart_lambda = config$drbart_lambda,
        drbart_nu = config$drbart_nu,
        drbart_kfac = config$drbart_kfac,
        drbart_phi0 = config$drbart_phi0,
        drbart_grid_tail_prob = config$drbart_grid_tail_prob,
        drbart_predict_batch_size = config$drbart_predict_batch_size,
        true_ps = true_ps,
        verbose = TRUE,
        seed = seed,
        overlap_strength = config$overlap_strength
    )

    pred_summary <- summarize_prediction_sets(data_simul$tautest, results$intervals_ITE_test)
    stage1_test <- stage1_test_diagnostics_aps(results, data_simul, config, n_train)
    wdiag <- results$diagnostics$stage1$weight_diagnostics

    result_row <- data.frame(
        seed = seed,
        scenario = config$scenario,
        rho = config$rho,
        method = format_method_label(config$stage2_method),
        stage1_mode = config$stage1_mode,
        stage1_model = results$params$stage1_model,
        stage2_method = results$params$stage2_method,
        alpha = config$alpha,
        stage1_n_grid = results$params$stage1_n_grid,
        stage1_grid_resolution = results$params$stage1_grid_resolution,
        stage1_max_grid_size = results$params$stage1_max_grid_size,
        stage2_n_grid = results$params$stage2_n_grid,
        stage2_grid_resolution = results$params$stage2_grid_resolution,
        stage2_max_grid_size = results$params$stage2_max_grid_size,
        gamma_score = results$params$gamma_score,
        requested_gamma_score = config$gamma_score,
        coverage = pred_summary$coverage,
        mean_length = pred_summary$mean_length,
        se_length = pred_summary$se_length,
        min_length = pred_summary$min_length,
        max_length = pred_summary$max_length,
        mean_n_intervals = pred_summary$mean_n_intervals,
        stage1_mean_length = results$diagnostics$stage1$avg_width,
        stage1_se_length = stats::sd(results$diagnostics$stage1$stage1_widths) /
            sqrt(length(results$diagnostics$stage1$stage1_widths)),
        stage1_coverage_Y1 = stage1_test$Y1_cov,
        stage1_length_Y1 = stage1_test$Y1_len,
        stage1_coverage_Y0 = stage1_test$Y0_cov,
        stage1_length_Y0 = stage1_test$Y0_len,
        ess_Y1 = results$diagnostics$stage1$calibration_ess_Y1,
        ess_Y0 = results$diagnostics$stage1$calibration_ess_Y0,
        weight_1_to_0_max = wdiag$source_1_target_0$max,
        weight_1_to_0_q99 = wdiag$source_1_target_0$q99,
        weight_0_to_1_max = wdiag$source_0_target_1$max,
        weight_0_to_1_q99 = wdiag$source_0_target_1$q99,
        dp_k_mean_Y1 = results$diagnostics$stage1$dp_k_mean_Y1,
        dp_k_mean_Y0 = results$diagnostics$stage1$dp_k_mean_Y0,
        mdn_n_components_Y1 = results$diagnostics$stage1$mdn_n_components_Y1,
        mdn_n_components_Y0 = results$diagnostics$stage1$mdn_n_components_Y0,
        mdn_epochs_trained_Y1 = results$diagnostics$stage1$mdn_epochs_trained_Y1,
        mdn_epochs_trained_Y0 = results$diagnostics$stage1$mdn_epochs_trained_Y0,
        mdn_best_val_loss_Y1 = results$diagnostics$stage1$mdn_best_val_loss_Y1,
        mdn_best_val_loss_Y0 = results$diagnostics$stage1$mdn_best_val_loss_Y0,
        mdn_final_train_loss_Y1 = results$diagnostics$stage1$mdn_final_train_loss_Y1,
        mdn_final_train_loss_Y0 = results$diagnostics$stage1$mdn_final_train_loss_Y0,
        mdn_final_val_loss_Y1 = results$diagnostics$stage1$mdn_final_val_loss_Y1,
        mdn_final_val_loss_Y0 = results$diagnostics$stage1$mdn_final_val_loss_Y0,
        nsf_n_bins_Y1 = results$diagnostics$stage1$nsf_n_bins_Y1,
        nsf_n_bins_Y0 = results$diagnostics$stage1$nsf_n_bins_Y0,
        nsf_epochs_trained_Y1 = results$diagnostics$stage1$nsf_epochs_trained_Y1,
        nsf_epochs_trained_Y0 = results$diagnostics$stage1$nsf_epochs_trained_Y0,
        nsf_best_val_loss_Y1 = results$diagnostics$stage1$nsf_best_val_loss_Y1,
        nsf_best_val_loss_Y0 = results$diagnostics$stage1$nsf_best_val_loss_Y0,
        nsf_final_train_loss_Y1 = results$diagnostics$stage1$nsf_final_train_loss_Y1,
        nsf_final_train_loss_Y0 = results$diagnostics$stage1$nsf_final_train_loss_Y0,
        nsf_final_val_loss_Y1 = results$diagnostics$stage1$nsf_final_val_loss_Y1,
        nsf_final_val_loss_Y0 = results$diagnostics$stage1$nsf_final_val_loss_Y0,
        nsf_tail_bound_Y1 = results$diagnostics$stage1$nsf_tail_bound_Y1,
        nsf_tail_bound_Y0 = results$diagnostics$stage1$nsf_tail_bound_Y0
    )
    if (config$stage1_mode != "oracle" && config$stage1_model == "mdn") {
        attr(result_row, "mdn") <- list(
            params = results$params$mdn,
            device_Y1 = results$diagnostics$stage1$mdn_device_Y1,
            device_Y0 = results$diagnostics$stage1$mdn_device_Y0
        )
    }
    if (config$stage1_mode != "oracle" && config$stage1_model == "drbart") {
        attr(result_row, "drbart") <- list(
            params = results$params$drbart,
            Y1 = results$diagnostics$stage1$drbart_Y1,
            Y0 = results$diagnostics$stage1$drbart_Y0
        )
    }
    result_row
}

bind_rows_fill <- function(rows) {
    cols <- unique(unlist(lapply(rows, names)))
    rows <- lapply(rows, function(x) {
        missing <- setdiff(cols, names(x))
        for (col in missing) x[[col]] <- NA
        x[cols]
    })
    do.call(rbind, rows)
}

config <- parse_args(commandArgs(trailingOnly = TRUE))

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(file.path(PROJECT_ROOT, "Models"))
source("05_pipeline_aps_classification.R")

data_dir <- file.path(PROJECT_ROOT, "Data", config$scenario, format_rho_label(config$rho))
output_dir <- build_output_dir(config)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

seeds <- if (config$run_single_seed) config$target_seed else 0:(config$n_simulations - 1)
all_results <- vector("list", length(seeds))
for (i in seq_along(seeds)) {
    seed <- seeds[i]
    results_df <- tryCatch(
        run_one(seed, config, data_dir),
        error = function(e) {
            cat(sprintf("[Seed %03d] ERROR: %s\n", seed, e$message))
            data.frame(seed = seed, scenario = config$scenario, rho = config$rho,
                       method = format_method_label(config$stage2_method),
                       stage1_mode = config$stage1_mode,
                       stage1_model = config$stage1_model, alpha = config$alpha,
                       stage2_method = config$stage2_method,
                       stage1_n_grid = config$stage1_n_grid,
                       stage1_grid_resolution = config$stage1_grid_resolution,
                       stage1_max_grid_size = config$stage1_max_grid_size,
                       stage2_n_grid = config$stage2_n_grid,
                       stage2_grid_resolution = config$stage2_grid_resolution,
                       stage2_max_grid_size = config$stage2_max_grid_size,
                       coverage = NA_real_, mean_length = NA_real_)
        }
    )
    all_results[[i]] <- results_df
    if(nrow(results_df)!=1L || any(!is.finite(results_df$coverage)) || any(!is.finite(results_df$mean_length)))
        stop(sprintf("Seed %03d failed; no result file written.",seed))
    save(results_df, file = file.path(output_dir, sprintf("result_%03d.RData", seed)))
}

if (!config$run_single_seed) {
    combined_results <- bind_rows_fill(all_results)
    save(combined_results, file = file.path(output_dir, "combined_results.RData"))
    write.csv(combined_results, file.path(output_dir, "combined_results.csv"), row.names = FALSE)
}

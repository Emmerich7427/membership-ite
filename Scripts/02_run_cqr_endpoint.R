#!/usr/bin/env Rscript
################################################################################
# Script: 02_run_cqr_endpoint.R
# Purpose: Run CQR_endpoint baseline experiments
################################################################################

get_script_dir <- function() {
    file_arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
    normalizePath(getwd())
}

SCRIPT_DIR <- get_script_dir()
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))

parse_args <- function(args) {
    config <- list(
        scenario = "S3", rho = 0, alpha = 0.1, n_simulations = 1,
        run_single_seed = FALSE, target_seed = NULL,
        ps_mode = "estimated_ps", quantile_stage2 = 0.5,
        overlap_strength = 1.0
    )
    i <- 1
    while (i <= length(args)) {
        arg <- args[i]
        if (arg == "--scenario" && i < length(args)) {
            config$scenario <- args[i + 1]; i <- i + 2
        } else if (arg == "--rho" && i < length(args)) {
            config$rho <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--alpha" && i < length(args)) {
            config$alpha <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--n-sim" && i < length(args)) {
            config$n_simulations <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--seed" && i < length(args)) {
            config$run_single_seed <- TRUE; config$target_seed <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--ps-mode" && i < length(args)) {
            config$ps_mode <- args[i + 1]; i <- i + 2
        } else if (arg == "--overlap-strength" && i < length(args)) {
            config$overlap_strength <- as.numeric(args[i + 1]); i <- i + 2
        } else {
            stop(sprintf("Unknown or incomplete argument: %s", arg))
        }
    }
    config
}

format_rho_label <- function(rho) paste0("rho_", format(rho, trim = TRUE, scientific = FALSE))
build_output_dir <- function(config) {
    file.path(PROJECT_ROOT, "Results", config$scenario, format_rho_label(config$rho),
              "CQR_endpoint", config$ps_mode)
}

load_data_simulation <- function(path) {
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    get("data_simul", envir = env)
}

stage1_test_diagnostics_cqr <- function(results, data_simul, config, n_train) {
    alpha_stage <- config$alpha / 2
    X_test <- data_simul$Xtest
    T_test <- data_simul$Ttest
    out <- list(Y1_cov = NA_real_, Y1_len = NA_real_, Y0_cov = NA_real_, Y0_len = NA_real_)

    idx_T0 <- which(T_test == 0)
    if (length(idx_T0) > 0) {
        sets_Y1 <- predict_outcome_cqr(
            results$fit_Y1, X_test[idx_T0, , drop = FALSE], alpha = alpha_stage,
            weights_test = results$weights_stage1$source_1_target_0[n_train + idx_T0]
        )
        out$Y1_cov <- mean(data_simul$Y1test[idx_T0] >= sets_Y1$lower &
                           data_simul$Y1test[idx_T0] <= sets_Y1$upper)
        out$Y1_len <- mean(sets_Y1$upper - sets_Y1$lower)
    }

    idx_T1 <- which(T_test == 1)
    if (length(idx_T1) > 0) {
        sets_Y0 <- predict_outcome_cqr(
            results$fit_Y0, X_test[idx_T1, , drop = FALSE], alpha = alpha_stage,
            weights_test = results$weights_stage1$source_0_target_1[n_train + idx_T1]
        )
        out$Y0_cov <- mean(data_simul$Y0test[idx_T1] >= sets_Y0$lower &
                           data_simul$Y0test[idx_T1] <= sets_Y0$upper)
        out$Y0_len <- mean(sets_Y0$upper - sets_Y0$lower)
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

    results <- run_cqr_endpoint(
        X_all, T_all, Y_all,
        fold1_idx = folds$fold1_idx,
        fold2_idx = folds$fold2_idx,
        test_idx = (n_train + 1):(n_train + n_test),
        scenario = config$scenario,
        ps_mode = config$ps_mode,
        alpha = config$alpha,
        quantile_stage2 = config$quantile_stage2,
        true_ps = c(data_simul$pstrain, data_simul$pstest),
        verbose = TRUE,
        seed = seed,
        overlap_strength = config$overlap_strength
    )

    pred_summary <- summarize_prediction_sets(data_simul$tautest, results$intervals_ITE_test)
    stage1_test <- stage1_test_diagnostics_cqr(results, data_simul, config, n_train)
    wdiag <- results$diagnostics$stage1$weight_diagnostics
    data.frame(
        seed = seed,
        scenario = config$scenario,
        rho = config$rho,
        method = "CQR_endpoint",
        ps_mode = config$ps_mode,
        alpha = config$alpha,
        coverage = pred_summary$coverage,
        mean_length = pred_summary$mean_length,
        se_length = pred_summary$se_length,
        min_length = pred_summary$min_length,
        max_length = pred_summary$max_length,
        mean_n_intervals = pred_summary$mean_n_intervals,
        stage1_mean_length = results$diagnostics$stage1$avg_width,
        stage1_coverage_Y1 = stage1_test$Y1_cov,
        stage1_length_Y1 = stage1_test$Y1_len,
        stage1_coverage_Y0 = stage1_test$Y0_cov,
        stage1_length_Y0 = stage1_test$Y0_len,
        weight_1_to_0_max = wdiag$source_1_target_0$max,
        weight_1_to_0_q99 = wdiag$source_1_target_0$q99,
        weight_0_to_1_max = wdiag$source_0_target_1$max,
        weight_0_to_1_q99 = wdiag$source_0_target_1$q99
    )
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
source("05_pipeline_cqr_endpoint.R")

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
                       method = "CQR_endpoint", ps_mode = config$ps_mode,
                       alpha = config$alpha, coverage = NA_real_, mean_length = NA_real_)
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

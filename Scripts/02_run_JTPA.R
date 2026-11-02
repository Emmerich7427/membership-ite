#!/usr/bin/env Rscript
# Usage: Rscript Scripts/02_run_JTPA.R membership|cqr SEED
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), ".."))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !args[1] %in% c("membership", "cqr"))
    stop("Usage: Rscript Scripts/02_run_JTPA.R membership|cqr SEED (1..10)")
seed <- as.numeric(args[2])
stopifnot(length(seed) == 1L, is.finite(seed), seed %in% 1:10)
method <- if (args[1] == "membership") "APS_classification" else "CQR_endpoint"
dest <- file.path(root, "Results/JTPA", method)
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
input <- new.env()
load(file.path(root, "Data/JTPA", sprintf("real_%03d.RData", seed)), input)
d <- input$data_simul

setwd(file.path(root, "Models"))
source(if (method == "APS_classification") "05_pipeline_aps_classification.R" else "05_pipeline_cqr_endpoint.R")
# Use the exact prediction and quantile functions used in the final JTPA fits.
source("03_outcome_JTPA.R")
source("06_JTPA_diagnostics.R")
folds <- split_data_by_treatment(d$Xtrain, d$Ytrain, d$Ttrain,
                                 fold_prop = 0.5, seed = seed * 1000 + 42)
common <- list(X = rbind(d$Xtrain, d$Xtest), T = c(d$Ttrain, d$Ttest),
    Y = c(d$Ytrain, d$Ytest), fold1_idx = folds$fold1_idx,
    fold2_idx = folds$fold2_idx, test_idx = nrow(d$Xtrain) + seq_len(nrow(d$Xtest)),
    scenario = "JTPA", alpha = 0.1, true_ps = NULL, verbose = TRUE, seed = seed)
if (method == "APS_classification") {
    config <- list(stage1_mode = "estimated_ps", stage1_model = "dp",
        stage1_n_grid = 500, stage1_grid_resolution = 0.05, stage1_max_grid_size = 2000,
        stage2_n_grid = 1000, stage2_grid_resolution = 0.02, stage2_max_grid_size = 5000,
        n_samples = 100, pos_ratio = 0.5, strategy = "region_uniform", boundary_ratio = 0.5,
        stage2_method = "mixture", gamma_score = 0,
        dp_n_fit_max = 10000, dp_n_iter = 250L, dp_burn = 125L, dp_thin = 5L,
        dp_alpha = 5, dp_kappa0 = 0.05, dp_nu0_extra = 3, dp_prior_scale = 1)
    fit <- do.call(run_aps_classification, c(common, config))
} else {
    config <- list(ps_mode = "estimated_ps", quantile_stage2 = 0.5)
    fit <- do.call(run_cqr_endpoint, c(common, config))
}
C <- jrd_predict_outcomes(fit, d, method, seed, config)
result <- jrd_evaluate(input$test_ids, d$Ttest, d$Ytest, C$C0, C$C1, fit$intervals_ITE_test)
result$row$seed <- seed
result$row$method <- method
saveRDS(result, file.path(dest, sprintf("result_%03d.rds", seed)))
cat("Saved ", method, " seed ", seed, "\n", sep = "")

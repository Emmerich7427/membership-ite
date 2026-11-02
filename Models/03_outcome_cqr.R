################################################################################
# Module: 03_outcome_cqr.R
# Purpose: Outcome-level CQR prediction sets
################################################################################

source("01_conformal_utils.R")
source("02_models.R")


fit_outcome_cqr <- function(X, Y,
                            weights = NULL,
                            X_cal = NULL,
                            Y_cal = NULL,
                            weights_cal = NULL,
                            train_prop = 0.75,
                            lower_quantile = 0.025,
                            upper_quantile = 0.975,
                            n_post = 200,
                            seed = NULL) {
    if (lower_quantile <= 0 || upper_quantile >= 1 || lower_quantile >= upper_quantile) {
        stop("Require 0 < lower_quantile < upper_quantile < 1.")
    }

    n <- length(Y)
    has_external_cal <- !is.null(X_cal) || !is.null(Y_cal) || !is.null(weights_cal)
    if (has_external_cal) {
        if (is.null(X_cal) || is.null(Y_cal)) {
            stop("X_cal and Y_cal must both be supplied when using external calibration data.")
        }
        X_train <- X
        Y_train <- Y
        X_val <- X_cal
        Y_val <- Y_cal
        train_idx <- seq_len(n)
        val_idx <- seq_along(Y_cal)
        weights_val <- if (is.null(weights_cal)) rep(1, length(Y_cal)) else weights_cal
        if (length(weights_val) != length(Y_cal)) {
            stop("weights_cal must have the same length as Y_cal.")
        }
    } else {
        n_train <- floor(n * train_prop)
        if (!is.null(seed)) {
            set.seed(seed)
        }
        train_idx <- sample(n, n_train)
        val_idx <- setdiff(seq_len(n), train_idx)

        X_train <- X[train_idx, , drop = FALSE]
        Y_train <- Y[train_idx]
        X_val <- X[val_idx, , drop = FALSE]
        Y_val <- Y[val_idx]
        weights_val <- if (is.null(weights)) rep(1, length(val_idx)) else weights[val_idx]
    }

    model_lower <- fit_quantile_bart(X_train, Y_train, lower_quantile, n_post = n_post)
    model_upper <- fit_quantile_bart(X_train, Y_train, upper_quantile, n_post = n_post)

    pred_lower_val <- predict_quantile_bart(model_lower, X_val)
    pred_upper_val <- predict_quantile_bart(model_upper, X_val)
    scores_val <- pmax(pred_lower_val - Y_val, Y_val - pred_upper_val)

    structure(
        list(
            model_lower = model_lower,
            model_upper = model_upper,
            scores_val = scores_val,
            weights_val = weights_val,
            train_idx = train_idx,
            val_idx = val_idx,
            split_source = if (has_external_cal) "external" else "internal",
            lower_quantile = lower_quantile,
            upper_quantile = upper_quantile,
            score_range = range(scores_val),
            score_mean = mean(scores_val)
        ),
        class = "outcome_cqr_fit"
    )
}


predict_outcome_cqr <- function(fit_obj,
                                X_test,
                                alpha = 0.05,
                                weights_test = NULL,
                                weight_upper = 20,
                                weight_lower = 0.05) {
    n_test <- nrow(X_test)
    pred_lower <- predict_quantile_bart(fit_obj$model_lower, X_test)
    pred_upper <- predict_quantile_bart(fit_obj$model_upper, X_test)

    if (is.null(weights_test)) {
        weights_test <- rep(1, n_test)
    }

    avg_weight <- mean(c(fit_obj$weights_val, weights_test))
    weights_val_trunc <- truncate_weights(fit_obj$weights_val / avg_weight, weight_upper, weight_lower)
    weights_test_trunc <- truncate_weights(weights_test / avg_weight, weight_upper, weight_lower)
    total_weight_val <- sum(weights_val_trunc)
    weights_val_final <- weights_val_trunc / total_weight_val

    slacks <- numeric(n_test)
    for (i in seq_len(n_test)) {
        qt <- (1 + weights_test_trunc[i] / total_weight_val) * (1 - alpha)
        slacks[i] <- weighted_quantile(fit_obj$scores_val, weights_val_final, qt, use_inf = TRUE)
    }

    data.frame(lower = pred_lower - slacks, upper = pred_upper + slacks)
}

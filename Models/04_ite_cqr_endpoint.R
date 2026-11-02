################################################################################
# Module: 04_ite_cqr_endpoint.R
# Purpose: Endpoint CQR model for induced ITE intervals
################################################################################

source("01_conformal_utils.R")
source("02_models.R")


fit_ite_cqr_endpoint <- function(X,
                                 Y_interval,
                                 train_prop = 0.75,
                                 quantile = 0.5,
                                 seed = NULL) {
    n <- nrow(Y_interval)
    n_train <- floor(n * train_prop)
    if (!is.null(seed)) {
        set.seed(seed)
    }
    train_idx <- sample(n, n_train)
    val_idx <- setdiff(seq_len(n), train_idx)

    X_train <- X[train_idx, , drop = FALSE]
    X_val <- X[val_idx, , drop = FALSE]
    lower_train <- Y_interval[train_idx, 1]
    upper_train <- Y_interval[train_idx, 2]
    lower_val <- Y_interval[val_idx, 1]
    upper_val <- Y_interval[val_idx, 2]

    model_lower <- fit_quantile_rf(X_train, lower_train, quantile)
    model_upper <- fit_quantile_rf(X_train, upper_train, quantile)

    pred_lower_val <- predict_quantile_rf(model_lower, X_val)
    pred_upper_val <- predict_quantile_rf(model_upper, X_val)
    scores_val <- pmax(pred_lower_val - lower_val, upper_val - pred_upper_val, 0)

    structure(
        list(
            model_lower = model_lower,
            model_upper = model_upper,
            scores_val = scores_val,
            train_idx = train_idx,
            val_idx = val_idx,
            quantile = quantile,
            score_range = range(scores_val),
            score_mean = mean(scores_val)
        ),
        class = "ite_cqr_endpoint_fit"
    )
}


predict_ite_cqr_endpoint <- function(fit_obj, X_test, alpha = 0.05) {
    pred_lower <- predict_quantile_rf(fit_obj$model_lower, X_test)
    pred_upper <- predict_quantile_rf(fit_obj$model_upper, X_test)
    n_cal <- length(fit_obj$scores_val)
    weights_val <- rep(1 / n_cal, n_cal)
    q_level <- (1 + 1 / n_cal) * (1 - alpha)
    slack <- weighted_quantile(fit_obj$scores_val, weights_val, q_level, use_inf = TRUE)
    data.frame(lower = pred_lower - slack, upper = pred_upper + slack)
}

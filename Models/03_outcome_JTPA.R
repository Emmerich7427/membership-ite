# Functions used for the final JTPA fits; source after the shared pipelines.
weighted_quantile <- function(scores, weights, quantile_level, use_inf = FALSE) {
    if (length(scores) != length(weights)) {
        stop("scores and weights must have the same length.")
    }
    if (any(weights < 0) || sum(weights) <= 0) {
        stop("weights must be nonnegative with positive total mass.")
    }

    ord <- order(scores)
    sorted_scores <- scores[ord]
    sorted_weights <- weights[ord] / sum(weights)
    qt_capped <- pmin(quantile_level, 1)
    idx <- find_inds(cumsum(sorted_weights), qt_capped)
    result <- sorted_scores[pmin(idx, length(sorted_scores))]

    if (use_inf && quantile_level >= 1) {
        result <- Inf
    }
    result
}

predict_outcome_aps <- function(fit_obj,
                                X_test,
                                alpha = 0.05,
                                weights_test = NULL,
                                weight_upper = 20,
                                weight_lower = 0.05,
                                n_grid = NULL,
                                grid_resolution = NULL,
                                max_grid_size = NULL,
                                ...) {
    n_test <- nrow(X_test)
    if (is.null(weights_test)) {
        weights_test <- rep(1, n_test)
    }
    grid_params <- resolve_outcome_prediction_grid_params(
        fit_obj,
        n_grid = n_grid,
        grid_resolution = grid_resolution,
        max_grid_size = max_grid_size
    )

    avg_weight <- mean(c(fit_obj$weights_val, weights_test))
    weights_val_trunc <- truncate_weights(fit_obj$weights_val / avg_weight, weight_upper, weight_lower)
    weights_test_trunc <- truncate_weights(weights_test / avg_weight, weight_upper, weight_lower)
    total_weight_val <- sum(weights_val_trunc)
    weights_val_final <- weights_val_trunc / total_weight_val

    mdn_params_test <- NULL
    nsf_params_test <- NULL
    if (inherits(fit_obj$density_model, "mdn_density")) {
        mdn_params_test <- predict_mdn_params(fit_obj$density_model, X_test)
    } else if (inherits(fit_obj$density_model, "nsf_density")) {
        nsf_params_test <- predict_nsf_params(fit_obj$density_model, X_test)
    }

    intervals_list <- vector("list", n_test)
    for (i in seq_len(n_test)) {
        qt <- (1 + weights_test_trunc[i] / total_weight_val) * (1 - alpha)
        cutoff_i <- weighted_quantile(fit_obj$scores_val, weights_val_final, qt)
        if (!is.null(mdn_params_test)) {
            y_grid <- mdn_y_grid_from_params(
                fit_obj$density_model,
                mdn_params_test,
                i,
                grid_params$n_grid,
                grid_params$grid_resolution,
                grid_params$max_grid_size
            )
            density_grid <- predict_mdn_density_from_params(
                fit_obj$density_model,
                mdn_params_test,
                i,
                y_grid
            )$density_mean
        } else if (!is.null(nsf_params_test)) {
            y_grid <- nsf_y_grid_from_params(
                fit_obj$density_model,
                nsf_params_test,
                i,
                grid_params$n_grid,
                grid_params$grid_resolution,
                grid_params$max_grid_size
            )
            density_grid <- predict_nsf_density_from_params(
                fit_obj$density_model,
                nsf_params_test,
                i,
                y_grid
            )$density_mean
        } else {
            y_grid <- outcome_y_grid_for_x(
                fit_obj$density_model,
                X_test[i, , drop = FALSE],
                grid_params$n_grid,
                grid_params$grid_resolution,
                grid_params$max_grid_size
            )
            density_grid <- predict_outcome_density(
                fit_obj$density_model,
                X_test[i, , drop = FALSE],
                y_grid
            )$density_mean
        }
        scores_grid <- compute_aps_scores_grid(density_grid, y_grid)
        density_threshold <- aps_density_threshold(density_grid, scores_grid, cutoff_i)
        intervals_list[[i]] <- linear_superlevel_intervals(
            grid = y_grid,
            values = density_grid,
            threshold = density_threshold,
            left_value = density_grid[1],
            right_value = density_grid[length(density_grid)]
        )
    }
    intervals_list
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
        slacks[i] <- weighted_quantile(fit_obj$scores_val, weights_val_final, qt)
    }

    data.frame(lower = pred_lower - slacks, upper = pred_upper + slacks)
}

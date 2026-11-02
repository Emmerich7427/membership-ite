################################################################################
# Module: 04_ite_aps_classification.R
# Purpose: Probit BART classification for ITE interval unions
################################################################################

source("01_interval_utils.R")
source("01_conformal_utils.R")


construct_plausible_region <- function(interval_unions, margin_ratio = 0.2) {
    normalized <- lapply(interval_unions, normalize_interval_union)
    nonempty <- normalized[sapply(normalized, nrow) > 0]
    if (length(nonempty) == 0) {
        stop("Cannot construct plausible region from empty interval unions.")
    }

    all_lower <- unlist(lapply(nonempty, function(u) u[, 1]))
    all_upper <- unlist(lapply(nonempty, function(u) u[, 2]))
    if (any(!is.finite(c(all_lower, all_upper)))) {
        stop("Cannot construct finite plausible region from unbounded interval unions.")
    }

    tau_min <- min(all_lower)
    tau_max <- max(all_upper)
    margin <- margin_ratio * (tau_max - tau_min)
    c(tau_min - margin, tau_max + margin)
}


construct_complement <- function(region, union_sorted) {
    union_sorted <- normalize_interval_union(union_sorted)
    complement <- matrix(numeric(0), ncol = 2)

    if (nrow(union_sorted) == 0) {
        return(matrix(region, nrow = 1))
    }
    if (union_sorted[1, 1] > region[1]) {
        complement <- rbind(complement, c(region[1], union_sorted[1, 1]))
    }
    if (nrow(union_sorted) > 1) {
        for (k in seq_len(nrow(union_sorted) - 1)) {
            if (union_sorted[k + 1, 1] > union_sorted[k, 2]) {
                complement <- rbind(complement, c(union_sorted[k, 2], union_sorted[k + 1, 1]))
            }
        }
    }
    if (union_sorted[nrow(union_sorted), 2] < region[2]) {
        complement <- rbind(complement, c(union_sorted[nrow(union_sorted), 2], region[2]))
    }
    normalize_interval_union(complement)
}


intersect_interval_unions <- function(union_a, union_b) {
    union_a <- normalize_interval_union(union_a)
    union_b <- normalize_interval_union(union_b)
    intersections <- matrix(numeric(0), ncol = 2)
    if (nrow(union_a) == 0 || nrow(union_b) == 0) {
        return(empty_interval_union())
    }
    for (i in seq_len(nrow(union_a))) {
        for (j in seq_len(nrow(union_b))) {
            lower <- max(union_a[i, 1], union_b[j, 1])
            upper <- min(union_a[i, 2], union_b[j, 2])
            if (upper > lower) {
                intersections <- rbind(intersections, c(lower, upper))
            }
        }
    }
    normalize_interval_union(intersections)
}


sample_from_complement_uniform <- function(complement, n_samples) {
    if (n_samples <= 0) {
        return(numeric(0))
    }
    complement <- normalize_interval_union(complement)
    if (nrow(complement) == 0) {
        return(numeric(0))
    }
    lengths <- complement[, 2] - complement[, 1]
    if (any(!is.finite(lengths))) {
        stop("Cannot sample from an unbounded complement.")
    }
    counts <- rmultinom(1, n_samples, lengths / sum(lengths))[, 1]
    samples <- numeric(0)
    for (k in seq_len(nrow(complement))) {
        if (counts[k] > 0) {
            samples <- c(samples, runif(counts[k], complement[k, 1], complement[k, 2]))
        }
    }
    samples
}


construct_boundary_band <- function(union_sorted, region, bandwidth_ratio = 0.02) {
    boundaries <- c(union_sorted[, 1], union_sorted[, 2])
    bandwidth <- bandwidth_ratio * (region[2] - region[1])
    normalize_interval_union(cbind(
        pmax(boundaries - bandwidth, region[1]),
        pmin(boundaries + bandwidth, region[2])
    ))
}


sample_classification_points <- function(region, interval_union,
                                         n_samples = 100,
                                         pos_ratio = 0.5,
                                         strategy = "region_uniform",
                                         boundary_ratio = 0.5) {
    union_sorted <- normalize_interval_union(interval_union)

    if (strategy == "region_uniform") {
        tau <- runif(n_samples, region[1], region[2])
        labels <- as.integer(points_in_union(tau, union_sorted))
        return(list(tau = tau, labels = labels))
    }

    complement <- construct_complement(region, union_sorted)
    if (nrow(complement) == 0 || union_total_length(complement) <= 0) {
        stop("Cannot sample negative examples from an empty complement.")
    }

    n_pos <- min(max(1, round(n_samples * pos_ratio)), n_samples - 1)
    n_neg <- n_samples - n_pos

    if (strategy == "uniform") {
        tau_pos <- if (union_total_length(union_sorted) > 0) {
            sample_from_union(union_sorted, n_pos)
        } else {
            rep(mean(union_sorted), n_pos)
        }
        tau_neg <- sample_from_complement_uniform(complement, n_neg)
    } else if (strategy == "boundary") {
        boundary_band <- construct_boundary_band(union_sorted, region)
        pos_band <- intersect_interval_unions(union_sorted, boundary_band)
        neg_band <- intersect_interval_unions(complement, boundary_band)
        n_pos_boundary <- min(floor(n_pos * boundary_ratio), n_pos)
        n_neg_boundary <- min(floor(n_neg * boundary_ratio), n_neg)
        tau_pos <- c(
            if (union_total_length(pos_band) > 0) sample_from_union(pos_band, n_pos_boundary) else numeric(0),
            if (union_total_length(union_sorted) > 0) {
                sample_from_union(union_sorted, n_pos - n_pos_boundary)
            } else {
                rep(mean(union_sorted), n_pos - n_pos_boundary)
            }
        )
        tau_neg <- c(
            if (union_total_length(neg_band) > 0) sample_from_complement_uniform(neg_band, n_neg_boundary) else numeric(0),
            sample_from_complement_uniform(complement, n_neg - n_neg_boundary)
        )
    } else {
        stop("strategy must be 'uniform', 'boundary', or 'region_uniform'.")
    }

    list(tau = c(tau_pos, tau_neg), labels = c(rep(1, length(tau_pos)), rep(0, length(tau_neg))))
}


stage2_grid_params <- function(n_grid = 500,
                               grid_resolution = 0.05,
                               max_grid_size = 2000) {
    if (length(n_grid) != 1 || is.na(n_grid) || n_grid < 2) {
        stop("n_grid must be at least 2.")
    }
    if (length(grid_resolution) != 1 || is.na(grid_resolution) || grid_resolution <= 0) {
        stop("grid_resolution must be positive.")
    }
    if (length(max_grid_size) != 1 || is.na(max_grid_size) || max_grid_size < 2) {
        stop("max_grid_size must be at least 2.")
    }

    list(
        n_grid = as.integer(n_grid),
        grid_resolution = as.numeric(grid_resolution),
        max_grid_size = as.integer(max_grid_size)
    )
}


resolve_ite_prediction_grid <- function(fit_obj,
                                        n_grid = NULL,
                                        grid_resolution = NULL,
                                        max_grid_size = NULL) {
    default_params <- stage2_grid_params()
    fit_params <- if (is.null(fit_obj$grid_params)) {
        default_params
    } else {
        stage2_grid_params(
            fit_obj$grid_params$n_grid,
            fit_obj$grid_params$grid_resolution,
            fit_obj$grid_params$max_grid_size
        )
    }

    supplied <- list(
        n_grid = n_grid,
        grid_resolution = grid_resolution,
        max_grid_size = max_grid_size
    )
    resolved <- fit_params

    for (name in names(supplied)) {
        if (!is.null(supplied[[name]])) {
            candidate <- resolved
            candidate[[name]] <- supplied[[name]]
            candidate <- stage2_grid_params(
                candidate$n_grid,
                candidate$grid_resolution,
                candidate$max_grid_size
            )

            if (!is.null(fit_obj$grid_params)) {
                same <- if (name == "grid_resolution") {
                    isTRUE(all.equal(candidate[[name]], fit_params[[name]], tolerance = 1e-12))
                } else {
                    identical(as.integer(candidate[[name]]), as.integer(fit_params[[name]]))
                }
                if (!same) {
                    stop(sprintf(
                        "Prediction-time %s (%s) does not match fit-time value (%s).",
                        name,
                        as.character(candidate[[name]]),
                        as.character(fit_params[[name]])
                    ))
                }
            }

            resolved <- candidate
        }
    }

    if (!is.null(fit_obj$tau_grid)) {
        return(as.numeric(fit_obj$tau_grid))
    }

    make_regular_grid(
        fit_obj$region,
        n_grid = resolved$n_grid,
        grid_resolution = resolved$grid_resolution,
        max_grid_size = resolved$max_grid_size
    )
}


predict_membership_on_grid <- function(bart_fit,
                                       X_row,
                                       tau_grid,
                                       X_names) {
    tau_grid <- as.numeric(tau_grid)
    if (length(tau_grid) == 0) {
        return(numeric(0))
    }

    x_row <- as.numeric(X_row[1, ])
    if (length(x_row) != length(X_names)) {
        stop("X_row and X_names have incompatible dimensions.")
    }

    X_pred <- as.data.frame(matrix(rep(x_row, each = length(tau_grid)),
                                  nrow = length(tau_grid)))
    colnames(X_pred) <- X_names
    X_pred$tau <- tau_grid

    f_draws <- predict(bart_fit, newdata = X_pred)
    if (is.null(dim(f_draws))) {
        p_hat <- as.numeric(f_draws)
    } else {
        p_hat <- colMeans(f_draws)
    }

    if (length(p_hat) != length(tau_grid)) {
        stop("BART prediction output does not match the requested tau grid.")
    }
    if (any(is.na(p_hat))) {
        stop("BART prediction output contains NA or NaN values.")
    }

    as.numeric(p_hat)
}


fit_ite_aps_classification <- function(X,
                                       interval_unions,
                                       train_prop = 0.75,
                                       n_samples = 100,
                                       pos_ratio = 0.5,
                                       strategy = "region_uniform",
                                       boundary_ratio = 0.5,
                                       gamma_score = 0,
                                       n_grid = 500,
                                       grid_resolution = 0.05,
                                       max_grid_size = 2000,
                                       verbose = FALSE,
                                       seed = NULL) {
    if (!requireNamespace("dbarts", quietly = TRUE)) {
        stop("Package 'dbarts' is required for APS classification.")
    }
    if (gamma_score < 0 || gamma_score >= 1) {
        stop("gamma_score must be in [0, 1).")
    }
    grid_params <- stage2_grid_params(n_grid, grid_resolution, max_grid_size)

    n <- nrow(X)
    n_train <- floor(n * train_prop)
    if (!is.null(seed)) {
        set.seed(seed)
    }
    train_idx <- sample(n, n_train)
    val_idx <- setdiff(seq_len(n), train_idx)
    X_train <- X[train_idx, , drop = FALSE]
    X_val <- X[val_idx, , drop = FALSE]
    unions_train <- interval_unions[train_idx]
    unions_val <- interval_unions[val_idx]

    region <- construct_plausible_region(unions_train)
    tau_grid <- make_regular_grid(
        region,
        n_grid = grid_params$n_grid,
        grid_resolution = grid_params$grid_resolution,
        max_grid_size = grid_params$max_grid_size
    )
    if (verbose) {
        cat(sprintf("  Plausible ITE region: [%.3f, %.3f]\n", region[1], region[2]))
        cat(sprintf("  Canonical ITE grid size: %d\n", length(tau_grid)))
    }

    tau_all <- c()
    x_ids <- c()
    labels <- c()
    for (i in seq_len(n_train)) {
        sampled <- sample_classification_points(
            region, unions_train[[i]], n_samples, pos_ratio, strategy, boundary_ratio
        )
        tau_all <- c(tau_all, sampled$tau)
        x_ids <- c(x_ids, rep(i, length(sampled$labels)))
        labels <- c(labels, sampled$labels)
    }

    X_expanded <- as.data.frame(X_train[x_ids, , drop = FALSE])
    colnames(X_expanded) <- paste0("X", seq_len(ncol(X_expanded)))
    X_names <- colnames(X_expanded)
    X_expanded$tau <- tau_all

    bart_fit <- dbarts::bart(
        x.train = X_expanded,
        y.train = labels,
        keeptrees = TRUE,
        verbose = FALSE
    )

    scores_val <- numeric(length(val_idx))
    n_calibration_sets_outside_region <- 0L
    for (i in seq_along(val_idx)) {
        union_i <- normalize_interval_union(unions_val[[i]])
        if (nrow(union_i) == 0) {
            scores_val[i] <- 0
            next
        }

        if (any(union_i[, 1] < region[1] | union_i[, 2] > region[2])) {
            n_calibration_sets_outside_region <- n_calibration_sets_outside_region + 1L
        }

        p_hat <- predict_membership_on_grid(
            bart_fit,
            X_val[i, , drop = FALSE],
            tau_grid,
            X_names
        )
        membership_quantile <- linear_quantile_over_union(
            grid = tau_grid,
            values = p_hat,
            interval_union = union_i,
            gamma = gamma_score,
            region = region,
            outside_value = 0
        )
        scores_val[i] <- 1 - membership_quantile
    }

    structure(
        list(
            bart_fit = bart_fit,
            scores_val = scores_val,
            weights_val = rep(1, length(val_idx)),
            train_idx = train_idx,
            val_idx = val_idx,
            X_names = X_names,
            region = region,
            tau_grid = tau_grid,
            grid_params = grid_params,
            grid_diagnostics = list(
                n_grid_actual = length(tau_grid),
                grid_spacing = diff(tau_grid),
                max_grid_spacing = max(diff(tau_grid)),
                n_calibration_sets_outside_region = n_calibration_sets_outside_region,
                prop_calibration_sets_outside_region =
                    n_calibration_sets_outside_region / max(1, length(val_idx))
            ),
            classification_params = list(
                n_samples = n_samples,
                pos_ratio = pos_ratio,
                strategy = strategy,
                boundary_ratio = boundary_ratio,
                gamma_score = gamma_score
            ),
            score_range = range(scores_val),
            score_mean = mean(scores_val)
        ),
        class = "ite_aps_classification_fit"
    )
}


predict_ite_aps_classification <- function(fit_obj,
                                           X_test,
                                           alpha = 0.05,
                                           n_grid = NULL,
                                           grid_resolution = NULL,
                                           max_grid_size = NULL) {
    n_test <- nrow(X_test)
    region <- fit_obj$region
    n_cal <- length(fit_obj$scores_val)
    weights_val <- rep(1 / n_cal, n_cal)
    q_level <- (1 + 1 / n_cal) * (1 - alpha)
    cutoff <- weighted_quantile(fit_obj$scores_val, weights_val, q_level, use_inf = TRUE)

    tau_grid <- resolve_ite_prediction_grid(
        fit_obj,
        n_grid = n_grid,
        grid_resolution = grid_resolution,
        max_grid_size = max_grid_size
    )
    intervals_list <- vector("list", n_test)

    for (i in seq_len(n_test)) {
        p_hat <- predict_membership_on_grid(
            fit_obj$bart_fit,
            X_test[i, , drop = FALSE],
            tau_grid,
            fit_obj$X_names
        )
        intervals_list[[i]] <- if (cutoff == Inf || cutoff >= 1) {
            normalize_interval_union(matrix(c(-Inf, Inf), nrow = 1))
        } else {
            linear_superlevel_intervals(
                grid = tau_grid,
                values = p_hat,
                threshold = 1 - cutoff,
                left_value = 0,
                right_value = 0
            )
        }
    }

    intervals_list
}


fit_ite_aps_classification_min <- function(X,
                                           T,
                                           interval_unions,
                                           train_prop = 0.75,
                                           n_samples = 100,
                                           pos_ratio = 0.5,
                                           strategy = "region_uniform",
                                           boundary_ratio = 0.5,
                                           gamma_score = 0,
                                           n_grid = 500,
                                           grid_resolution = 0.05,
                                           max_grid_size = 2000,
                                           verbose = FALSE,
                                           seed = NULL) {
    if (!requireNamespace("dbarts", quietly = TRUE)) {
        stop("Package 'dbarts' is required for APS classification.")
    }
    if (gamma_score < 0 || gamma_score >= 1) {
        stop("gamma_score must be in [0, 1).")
    }
    grid_params <- stage2_grid_params(n_grid, grid_resolution, max_grid_size)

    n <- nrow(X)
    n_train <- floor(n * train_prop)
    if (!is.null(seed)) {
        set.seed(seed)
    }
    train_idx <- sample(n, n_train)
    val_idx <- setdiff(seq_len(n), train_idx)
    X_train <- X[train_idx, , drop = FALSE]
    X_val <- X[val_idx, , drop = FALSE]
    T_train <- T[train_idx]
    unions_train <- interval_unions[train_idx]
    unions_val <- interval_unions[val_idx]

    region <- construct_plausible_region(unions_train)
    tau_grid <- make_regular_grid(
        region,
        n_grid = grid_params$n_grid,
        grid_resolution = grid_params$grid_resolution,
        max_grid_size = grid_params$max_grid_size
    )
    if (verbose) {
        cat(sprintf("  Plausible ITE region: [%.3f, %.3f]\n", region[1], region[2]))
        cat(sprintf("  Canonical ITE grid size: %d\n", length(tau_grid)))
    }

    fit_arm <- function(arm) {
        arm_idx <- which(T_train == arm)
        if (length(arm_idx) == 0) {
            stop(sprintf("No Stage 2 training observations in arm %s.", arm))
        }

        tau_all <- c()
        x_ids <- c()
        labels <- c()
        for (i in arm_idx) {
            sampled <- sample_classification_points(
                region, unions_train[[i]], n_samples, pos_ratio, strategy, boundary_ratio
            )
            tau_all <- c(tau_all, sampled$tau)
            x_ids <- c(x_ids, rep(i, length(sampled$labels)))
            labels <- c(labels, sampled$labels)
        }

        X_expanded <- as.data.frame(X_train[x_ids, , drop = FALSE])
        colnames(X_expanded) <- paste0("X", seq_len(ncol(X_expanded)))
        X_names <- colnames(X_expanded)
        X_expanded$tau <- tau_all

        bart_fit <- dbarts::bart(
            x.train = X_expanded,
            y.train = labels,
            keeptrees = TRUE,
            verbose = FALSE
        )

        list(bart_fit = bart_fit, X_names = X_names, n_train = length(arm_idx))
    }

    fit_0 <- fit_arm(0)
    fit_1 <- fit_arm(1)
    X_names <- fit_0$X_names

    scores_val <- numeric(length(val_idx))
    n_calibration_sets_outside_region <- 0L
    for (i in seq_along(val_idx)) {
        union_i <- normalize_interval_union(unions_val[[i]])
        if (nrow(union_i) == 0) {
            scores_val[i] <- 0
            next
        }

        if (any(union_i[, 1] < region[1] | union_i[, 2] > region[2])) {
            n_calibration_sets_outside_region <- n_calibration_sets_outside_region + 1L
        }

        p_hat <- pmin(
            predict_membership_on_grid(
                fit_0$bart_fit,
                X_val[i, , drop = FALSE],
                tau_grid,
                X_names
            ),
            predict_membership_on_grid(
                fit_1$bart_fit,
                X_val[i, , drop = FALSE],
                tau_grid,
                X_names
            )
        )
        membership_quantile <- linear_quantile_over_union(
            grid = tau_grid,
            values = p_hat,
            interval_union = union_i,
            gamma = gamma_score,
            region = region,
            outside_value = 0
        )
        scores_val[i] <- 1 - membership_quantile
    }

    structure(
        list(
            bart_fit_0 = fit_0$bart_fit,
            bart_fit_1 = fit_1$bart_fit,
            scores_val = scores_val,
            weights_val = rep(1, length(val_idx)),
            train_idx = train_idx,
            val_idx = val_idx,
            X_names = X_names,
            region = region,
            tau_grid = tau_grid,
            grid_params = grid_params,
            grid_diagnostics = list(
                n_grid_actual = length(tau_grid),
                grid_spacing = diff(tau_grid),
                max_grid_spacing = max(diff(tau_grid)),
                n_calibration_sets_outside_region = n_calibration_sets_outside_region,
                prop_calibration_sets_outside_region =
                    n_calibration_sets_outside_region / max(1, length(val_idx))
            ),
            classification_params = list(
                n_samples = n_samples,
                pos_ratio = pos_ratio,
                strategy = strategy,
                boundary_ratio = boundary_ratio,
                gamma_score = gamma_score,
                stage2_method = "min"
            ),
            n_train_arm0 = fit_0$n_train,
            n_train_arm1 = fit_1$n_train,
            score_range = range(scores_val),
            score_mean = mean(scores_val)
        ),
        class = "ite_aps_classification_min_fit"
    )
}


predict_ite_aps_classification_min <- function(fit_obj,
                                               X_test,
                                               alpha = 0.05,
                                               n_grid = NULL,
                                               grid_resolution = NULL,
                                               max_grid_size = NULL) {
    n_test <- nrow(X_test)
    region <- fit_obj$region
    n_cal <- length(fit_obj$scores_val)
    weights_val <- rep(1 / n_cal, n_cal)
    q_level <- (1 + 1 / n_cal) * (1 - alpha)
    cutoff <- weighted_quantile(fit_obj$scores_val, weights_val, q_level, use_inf = TRUE)

    tau_grid <- resolve_ite_prediction_grid(
        fit_obj,
        n_grid = n_grid,
        grid_resolution = grid_resolution,
        max_grid_size = max_grid_size
    )
    intervals_list <- vector("list", n_test)

    for (i in seq_len(n_test)) {
        p_hat <- pmin(
            predict_membership_on_grid(
                fit_obj$bart_fit_0,
                X_test[i, , drop = FALSE],
                tau_grid,
                fit_obj$X_names
            ),
            predict_membership_on_grid(
                fit_obj$bart_fit_1,
                X_test[i, , drop = FALSE],
                tau_grid,
                fit_obj$X_names
            )
        )
        intervals_list[[i]] <- if (cutoff == Inf || cutoff >= 1) {
            normalize_interval_union(matrix(c(-Inf, Inf), nrow = 1))
        } else {
            linear_superlevel_intervals(
                grid = tau_grid,
                values = p_hat,
                threshold = 1 - cutoff,
                left_value = 0,
                right_value = 0
            )
        }
    }

    intervals_list
}

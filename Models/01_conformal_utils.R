################################################################################
# Module: 01_conformal_utils.R
# Purpose: Weighted conformal utilities and weight diagnostics
################################################################################

find_inds <- function(a, b) {
    n <- length(a)
    b <- b - 1e-12
    combined <- c(a, b)
    rank_combined <- rank(combined, ties.method = "first")
    rank_b <- rank(b, ties.method = "first")
    rank_b_in_combined <- rank_combined[-seq_len(n)]
    rank_b_in_combined - rank_b + 1
}


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

    if (use_inf && quantile_level > 1) {
        result <- Inf
    }
    result
}


truncate_weights <- function(weights, upper = 20, lower = 0.05) {
    pmin(pmax(weights, lower), upper)
}


compute_conformal_weights <- function(propensity_scores,
                                      train_treatment,
                                      test_treatment) {
    e_x <- pmin(pmax(propensity_scores, 1e-6), 1 - 1e-6)

    if (train_treatment == test_treatment) {
        rep(1, length(e_x))
    } else if (train_treatment == 1 && test_treatment == 0) {
        (1 - e_x) / e_x
    } else if (train_treatment == 0 && test_treatment == 1) {
        e_x / (1 - e_x)
    } else {
        stop("Treatment values must be 0 or 1.")
    }
}


effective_sample_size <- function(weights) {
    weights <- weights[is.finite(weights) & weights >= 0]
    if (length(weights) == 0 || sum(weights) <= 0) {
        return(NA_real_)
    }
    sum(weights)^2 / sum(weights^2)
}


weight_tail_summary <- function(weights, upper = 20, lower = 0.05) {
    weights <- weights[is.finite(weights)]
    if (length(weights) == 0) {
        return(list(
            mean = NA_real_, sd = NA_real_, min = NA_real_, max = NA_real_,
            q95 = NA_real_, q99 = NA_real_, ess = NA_real_,
            prop_above_upper = NA_real_, prop_below_lower = NA_real_
        ))
    }

    list(
        mean = mean(weights),
        sd = stats::sd(weights),
        min = min(weights),
        max = max(weights),
        q95 = as.numeric(stats::quantile(weights, 0.95, names = FALSE)),
        q99 = as.numeric(stats::quantile(weights, 0.99, names = FALSE)),
        ess = effective_sample_size(weights),
        prop_above_upper = mean(weights > upper),
        prop_below_lower = mean(weights < lower)
    )
}


split_data_by_treatment <- function(X, Y, T, fold_prop = 0.5, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    treated_idx <- which(T == 1)
    control_idx <- which(T == 0)
    n_treated_fold1 <- ceiling(length(treated_idx) * fold_prop)
    n_control_fold1 <- ceiling(length(control_idx) * fold_prop)

    fold1_idx <- c(
        sample(treated_idx, n_treated_fold1),
        sample(control_idx, n_control_fold1)
    )
    fold2_idx <- setdiff(seq_along(Y), fold1_idx)

    list(fold1_idx = fold1_idx, fold2_idx = fold2_idx)
}

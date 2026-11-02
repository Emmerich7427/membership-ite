# Held-out diagnostics used for the final JTPA results.
jrd_sets <- function(x) {
    if (is.data.frame(x) || is.matrix(x)) {
        x <- as.matrix(x)
        if (ncol(x) != 2L) stop("Expected two interval endpoints.")
        x <- lapply(seq_len(nrow(x)), function(i) x[i, , drop = FALSE])
    }
    if (!is.list(x) || any(vapply(x, is.null, logical(1))) ||
        any(vapply(x, anyNA, logical(1)))) stop("Invalid prediction sets.")
    lapply(x, normalize_interval_union)
}

jrd_subset <- function(inner, outer, tol = 1e-8) {
    if (anyNA(inner) || anyNA(outer)) stop("Missing endpoint.")
    inner <- normalize_interval_union(inner)
    outer <- normalize_interval_union(outer)
    if (!nrow(inner)) return(TRUE)
    if (!nrow(outer)) return(FALSE)
    all(vapply(seq_len(nrow(inner)), function(i)
        any(outer[, 1] <= inner[i, 1] + tol &
            outer[, 2] >= inner[i, 2] - tol), logical(1)))
}

jrd_contains <- function(set, y) {
    if (anyNA(set) || !is.finite(y)) stop("Invalid factual outcome or set.")
    nrow(set) > 0L && any(set[, 1] <= y & y <= set[, 2])
}

jrd_induced <- function(T, Y, C0, C1) {
    if (length(T) != length(Y) || length(C0) != length(Y) ||
        length(C1) != length(Y) || anyNA(T) || any(!T %in% c(0, 1)) ||
        any(!is.finite(Y))) stop("Invalid induced-set inputs.")
    lapply(seq_along(Y), function(i) {
        s <- if (T[i] == 1) C0[[i]] else C1[[i]]
        if (T[i] == 1) s <- cbind(Y[i] - s[, 2], Y[i] - s[, 1])
        else s <- s - Y[i]
        normalize_interval_union(s)
    })
}

jrd_evaluate <- function(ids, T, Y, C0, C1, C_delta) {
    n <- length(ids)
    if (anyDuplicated(ids) || anyNA(ids) || length(T) != n || length(Y) != n ||
        !all(c(0, 1) %in% T)) stop("Invalid test sample.")
    C0 <- jrd_sets(C0); C1 <- jrd_sets(C1); C_delta <- jrd_sets(C_delta)
    if (any(c(length(C0), length(C1), length(C_delta)) != n)) stop("Misaligned sets.")
    S <- jrd_induced(T, Y, C0, C1)
    factual <- vapply(seq_len(n), function(i)
        jrd_contains(if (T[i] == 1) C1[[i]] else C0[[i]], Y[i]), logical(1))
    included <- vapply(seq_len(n), function(i) jrd_subset(S[[i]], C_delta[[i]]), logical(1))
    # Keep exact-endpoint inclusion as a numerical sensitivity check as well.
    included_exact <- vapply(seq_len(n), function(i)
        jrd_subset(S[[i]], C_delta[[i]], tol = 0), logical(1))
    subjects <- data.frame(subject_id = as.character(ids), treatment = T, outcome = Y,
        factual_covered = factual, included = included, included_exact = included_exact,
        length = vapply(C_delta, union_total_length, numeric(1)),
        n_intervals = vapply(C_delta, nrow, integer(1)),
        induced_empty = vapply(S, nrow, integer(1)) == 0L,
        contains_zero = vapply(C_delta, jrd_contains, logical(1), y = 0))
    row <- data.frame(n_test = n, n_control = sum(T == 0), n_treated = sum(T == 1),
        covered_control = sum(factual[T == 0]), covered_treated = sum(factual[T == 1]),
        stage1_coverage_control = mean(factual[T == 0]),
        stage1_coverage_treated = mean(factual[T == 1]),
        included_count = sum(included), inclusion_rate = mean(included),
        inclusion_rate_exact = mean(included_exact),
        inclusion_tolerance_changes = sum(included != included_exact),
        mean_length = mean(subjects$length),
        prop_disconnected = mean(subjects$n_intervals > 1L),
        prop_empty = mean(subjects$n_intervals == 0L),
        prop_infinite = mean(is.infinite(subjects$length)),
        prop_induced_empty = mean(subjects$induced_empty),
        prop_contains_zero = mean(subjects$contains_zero), ite_coverage = NA_real_)
    list(row = row, subjects = subjects,
         predictions = list(C0 = C0, C1 = C1, S_obs = S, C_delta = C_delta))
}

jrd_predict_outcomes <- function(fit, data, method, seed, config) {
    # Called only AFTER the historical pipeline has finished all predictions.
    # Dedicated seeds make extra posterior-predictive draws reproducible, without
    # consuming any randomness used by the historical fitting/calibration calls.
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) saved_seed <- .Random.seed
    on.exit(if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv)
            else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    n <- nrow(data$Xtest)
    test_idx <- nrow(data$Xtrain) + seq_len(n)
    C <- list(C0 = vector("list", n), C1 = vector("list", n))
    audit <- list()
    for (a in 0:1) for (observed_arm in 0:1) {
        ii <- which(data$Ttest == observed_arm)
        if (!length(ii)) stop("Both test arms must be present.")
        outcome_fit <- fit[[paste0("fit_Y", a)]]
        w <- fit$weights_stage1[[if (a == 0) "source_0_target_1" else "source_1_target_0"]][test_idx[ii]]
        # Preserve the original weighted rule, including batch-specific scaling.
        # Batch by observed arm, as the original fold-2 predictions do.
        avg <- mean(c(outcome_fit$weights_val, w))
        wcal <- truncate_weights(outcome_fit$weights_val / avg, 20, 0.05)
        wt <- truncate_weights(w / avg, 20, 0.05)
        qt <- (1 + wt / sum(wcal)) * (1 - fit$params$alpha_stage)
        diag_seed <- as.integer(1000000 + seed * 10 + a * 2 + observed_arm)
        set.seed(diag_seed)
        args <- list(fit_obj = outcome_fit, X_test = data$Xtest[ii, , drop = FALSE],
            alpha = fit$params$alpha_stage, weights_test = w, weight_upper = 20, weight_lower = 0.05)
        if (method == "APS_classification") {
            args$n_grid <- config$stage1_n_grid
            args$grid_resolution <- config$stage1_grid_resolution
            args$max_grid_size <- config$stage1_max_grid_size
        }
        predictor <- if (method == "APS_classification") predict_outcome_aps else predict_outcome_cqr
        cat(sprintf("[DIAGNOSTIC] outcome arm=%d; observed arm=%d; n=%d\n", a, observed_arm, length(ii)))
        C[[paste0("C", a)]][ii] <- jrd_sets(do.call(predictor, args))
        audit[[length(audit) + 1L]] <- data.frame(outcome_arm = a, observed_arm = observed_arm,
            n = length(ii), prediction_seed = diag_seed, weight_scale = avg,
            max_quantile_level = max(qt), quantile_above_one = sum(qt > 1))
    }
    C$weight_audit <- do.call(rbind, audit)
    C
}

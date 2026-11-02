################################################################################
# Module: 02_dgp_S3.R
# Purpose: Moderate-dimensional observational DGP for robustness experiments
################################################################################

gen_X_S3 <- function(n, p = 6) {
    X <- matrix(runif(n * p), nrow = n, ncol = p)
    colnames(X) <- paste0("X", seq_len(p))
    X
}


expit <- function(x) 1 / (1 + exp(-x))


softmax_rows <- function(logit_mat) {
    shifted <- logit_mat - apply(logit_mat, 1, max)
    exp_shifted <- exp(shifted)
    exp_shifted / rowSums(exp_shifted)
}


trimodal_base_centers_S3 <- c(-3, 0, 3)
trimodal_min_prob_S3 <- 0.08


mu0_S3 <- function(X) {
    X <- as.matrix(X)
    0.40 +
        0.70 * (X[, 1] - 0.5) -
        0.55 * (X[, 2] - 0.5) +
        0.45 * (X[, 3] - 0.5) +
        0.35 * (X[, 4] - 0.5) -
        0.30 * (X[, 5] - 0.5) +
        0.25 * (X[, 6] - 0.5) +
        0.70 * (X[, 1] - 0.5) * (X[, 2] - 0.5)
}


tau_S3 <- function(X) {
    X <- as.matrix(X)
    0.7 +
        0.35 * (X[, 1] - 0.5) -
        0.30 * (X[, 3] - 0.5) +
        0.25 * (X[, 4] - 0.5) +
        0.20 * (X[, 5] - 0.5) -
        0.20 * (X[, 6] - 0.5)
}


sigma_S3 <- function(X) {
    X <- as.matrix(X)
    0.08 + 0.030 * X[, 1] + 0.025 * X[, 4] +
        0.020 * abs(X[, 5] - 0.5) + 0.015 * X[, 6]
}


ps_S3 <- function(X, overlap_strength = 1.0, ...) {
    X <- as.matrix(X)
    eta <- -0.20 +
        1.20 * (X[, 1] - 0.5) -
        1.10 * (X[, 2] - 0.5) +
        0.90 * (X[, 3] - 0.5) +
        0.60 * (X[, 1] - 0.5) * (X[, 2] - 0.5)
    pmin(pmax(expit(overlap_strength * eta), 0.07), 0.93)
}


center_locations0_S3 <- function(X) {
    X <- as.matrix(X)
    cbind(
        trimodal_base_centers_S3[1] -
            0.25 * (X[, 1] - 0.5) + 0.15 * (X[, 4] - 0.5),
        trimodal_base_centers_S3[2] +
            0.20 * (X[, 2] - 0.5) - 0.15 * (X[, 5] - 0.5),
        trimodal_base_centers_S3[3] +
            0.20 * (X[, 3] - 0.5) + 0.15 * (X[, 6] - 0.5)
    )
}


center_locations1_S3 <- function(X) {
    X <- as.matrix(X)
    cbind(
        -1.20 -
            0.20 * (X[, 2] - 0.5) + 0.15 * (X[, 5] - 0.5),
        0.05 +
            0.20 * (X[, 1] - 0.5) - 0.15 * (X[, 4] - 0.5),
        1.35 +
            0.20 * (X[, 3] - 0.5) + 0.15 * (X[, 6] - 0.5)
    )
}


mix_probs0_S3 <- function(X) {
    X <- as.matrix(X)
    logits <- cbind(
        0.50 * (X[, 1] - 0.5) - 0.30 * (X[, 4] - 0.5),
        -0.20 * (X[, 2] - 0.5) + 0.30 * (X[, 5] - 0.5),
        -0.40 * (X[, 3] - 0.5) + 0.30 * (X[, 6] - 0.5)
    )
    trimodal_min_prob_S3 + (1 - 3 * trimodal_min_prob_S3) * softmax_rows(logits)
}


mix_probs1_S3 <- function(X) {
    X <- as.matrix(X)
    logits <- cbind(
        -0.30 * (X[, 1] - 0.5) + 0.40 * (X[, 5] - 0.5),
        0.40 * (X[, 2] - 0.5) - 0.25 * (X[, 4] - 0.5),
        0.30 * (X[, 3] - 0.5) + 0.25 * (X[, 6] - 0.5)
    )
    trimodal_min_prob_S3 + (1 - 3 * trimodal_min_prob_S3) * softmax_rows(logits)
}


component_sds0_S3 <- function(X) {
    X <- as.matrix(X)
    cbind(
        0.050 + 0.010 * X[, 1],
        0.060 + 0.010 * X[, 4],
        0.050 + 0.010 * X[, 6]
    )
}


component_sds1_S3 <- function(X) {
    X <- as.matrix(X)
    cbind(
        0.050 + 0.010 * X[, 2],
        0.060 + 0.010 * X[, 5],
        0.050 + 0.010 * X[, 6]
    )
}


centered_mixture_stats_S3 <- function(probs, raw_centers, component_sds) {
    center_mean <- rowSums(probs * raw_centers)
    centered_locations <- raw_centers - center_mean
    variance <- rowSums(probs * (centered_locations^2 + component_sds^2))
    list(centered_locations = centered_locations, variance = variance)
}


draw_centered_mixture_offset_S3 <- function(probs, raw_centers, component_sds) {
    stats <- centered_mixture_stats_S3(probs, raw_centers, component_sds)
    cumulative_probs <- t(apply(probs, 1, cumsum))
    u <- runif(nrow(probs))
    component_id <- 1 + (u > cumulative_probs[, 1]) + (u > cumulative_probs[, 2])
    rows <- seq_len(nrow(probs))
    stats$centered_locations[cbind(rows, component_id)] +
        component_sds[cbind(rows, component_id)] * rnorm(nrow(probs))
}


generate_data_S3 <- function(ntrain, ntest, seed, rho = 0,
                             p = 6, overlap_strength = 1.0, ...) {
    set.seed(seed)
    stopifnot(rho >= -1, rho <= 1)
    if (p != 6) {
        stop("S3 is designed for p = 6.")
    }

    n <- ntrain + ntest
    X <- gen_X_S3(n, p = p)
    mu0 <- mu0_S3(X)
    tau <- tau_S3(X)
    mu1 <- mu0 + tau
    sigma <- sigma_S3(X)
    ps <- ps_S3(X, overlap_strength = overlap_strength)

    centers0 <- center_locations0_S3(X)
    centers1 <- center_locations1_S3(X)
    probs0 <- mix_probs0_S3(X)
    probs1 <- mix_probs1_S3(X)
    comp_sds0 <- component_sds0_S3(X)
    comp_sds1 <- component_sds1_S3(X)

    offset0 <- draw_centered_mixture_offset_S3(probs0, centers0, comp_sds0)
    offset1 <- draw_centered_mixture_offset_S3(probs1, centers1, comp_sds1)

    z0 <- rnorm(n)
    z1 <- rnorm(n)
    eps0 <- z0
    eps1 <- rho * z0 + sqrt(1 - rho^2) * z1

    Y0 <- mu0 + offset0 + sigma * eps0
    Y1 <- mu1 + offset1 + sigma * eps1
    ite <- Y1 - Y0

    offset_var0 <- centered_mixture_stats_S3(probs0, centers0, comp_sds0)$variance
    offset_var1 <- centered_mixture_stats_S3(probs1, centers1, comp_sds1)$variance
    std <- sqrt(offset_var0 + offset_var1 + 2 * (1 - rho) * sigma^2)

    build_data_simul(X, Y0, Y1, tau, ite, std, ps, ntrain, ntest)
}


generate_data_S3_repeated <- function(ntrain, ntest, seed, outcomes_per_X,
                                      rho = 0, p = 6,
                                      overlap_strength = 1.0, ...) {
    set.seed(seed)
    stopifnot(rho >= -1, rho <= 1)
    if (p != 6) {
        stop("S3 is designed for p = 6.")
    }
    if (length(outcomes_per_X) != 1L || is.na(outcomes_per_X) ||
        outcomes_per_X < 1L || outcomes_per_X != as.integer(outcomes_per_X)) {
        stop("outcomes_per_X must be a positive integer.")
    }
    if (length(ntrain) != 1L || is.na(ntrain) || ntrain < 1L ||
        ntrain != as.integer(ntrain)) {
        stop("ntrain must be a positive integer.")
    }
    if (length(ntest) != 1L || is.na(ntest) || ntest < 1L ||
        ntest != as.integer(ntest)) {
        stop("ntest must be a positive integer.")
    }
    if (ntrain %% outcomes_per_X != 0L) {
        stop(sprintf(
            "ntrain = %d must be divisible by outcomes_per_X = %d.",
            ntrain, outcomes_per_X
        ))
    }

    n_unique_train <- ntrain %/% outcomes_per_X
    Xtrain_unique <- gen_X_S3(n_unique_train, p = p)
    Xtrain <- Xtrain_unique[
        rep(seq_len(n_unique_train), each = outcomes_per_X),
        ,
        drop = FALSE
    ]
    Xtest <- gen_X_S3(ntest, p = p)
    X <- rbind(Xtrain, Xtest)
    n <- nrow(X)

    mu0 <- mu0_S3(X)
    tau <- tau_S3(X)
    mu1 <- mu0 + tau
    sigma <- sigma_S3(X)
    ps <- ps_S3(X, overlap_strength = overlap_strength)

    centers0 <- center_locations0_S3(X)
    centers1 <- center_locations1_S3(X)
    probs0 <- mix_probs0_S3(X)
    probs1 <- mix_probs1_S3(X)
    comp_sds0 <- component_sds0_S3(X)
    comp_sds1 <- component_sds1_S3(X)

    offset0 <- draw_centered_mixture_offset_S3(probs0, centers0, comp_sds0)
    offset1 <- draw_centered_mixture_offset_S3(probs1, centers1, comp_sds1)

    z0 <- rnorm(n)
    z1 <- rnorm(n)
    eps0 <- z0
    eps1 <- rho * z0 + sqrt(1 - rho^2) * z1

    Y0 <- mu0 + offset0 + sigma * eps0
    Y1 <- mu1 + offset1 + sigma * eps1
    ite <- Y1 - Y0

    offset_var0 <- centered_mixture_stats_S3(probs0, centers0, comp_sds0)$variance
    offset_var1 <- centered_mixture_stats_S3(probs1, centers1, comp_sds1)$variance
    std <- sqrt(offset_var0 + offset_var1 + 2 * (1 - rho) * sigma^2)

    build_data_simul(X, Y0, Y1, tau, ite, std, ps, ntrain, ntest)
}


true_outcome_density_S3 <- function(X, y_grid, arm = 0, ...) {
    X <- as.matrix(X)
    if (nrow(X) != 1) {
        stop("true_outcome_density_S3 expects one row of X.")
    }

    if (arm == 0) {
        probs <- mix_probs0_S3(X)
        centers <- center_locations0_S3(X)
        comp_sds <- component_sds0_S3(X)
        base_mean <- mu0_S3(X)
    } else {
        probs <- mix_probs1_S3(X)
        centers <- center_locations1_S3(X)
        comp_sds <- component_sds1_S3(X)
        base_mean <- mu0_S3(X) + tau_S3(X)
    }

    stats <- centered_mixture_stats_S3(probs, centers, comp_sds)
    component_sd <- sqrt(comp_sds[1, ]^2 + sigma_S3(X)^2)
    density <- rep(0, length(y_grid))
    for (k in seq_len(3)) {
        density <- density + probs[1, k] * stats::dnorm(
            y_grid,
            mean = base_mean + stats$centered_locations[1, k],
            sd = component_sd[k]
        )
    }
    density
}


oracle_y_grid_S3 <- function(X, arm = 0, n_grid = 500,
                             grid_resolution = 0.05,
                             max_grid_size = 2000, ...) {
    X <- as.matrix(X)
    if (arm == 0) {
        probs <- mix_probs0_S3(X)
        centers <- center_locations0_S3(X)
        comp_sds <- component_sds0_S3(X)
        base_mean <- mu0_S3(X)
    } else {
        probs <- mix_probs1_S3(X)
        centers <- center_locations1_S3(X)
        comp_sds <- component_sds1_S3(X)
        base_mean <- mu0_S3(X) + tau_S3(X)
    }
    stats <- centered_mixture_stats_S3(probs, centers, comp_sds)
    component_sd <- sqrt(comp_sds[1, ]^2 + sigma_S3(X)^2)
    means <- as.numeric(base_mean + stats$centered_locations[1, ])
    lower <- min(means - 5 * component_sd)
    upper <- max(means + 5 * component_sd)
    n_grid_i <- min(max(n_grid, ceiling((upper - lower) / grid_resolution)), max_grid_size)
    seq(lower, upper, length.out = n_grid_i)
}

################################################################################
# Module: 02_dgp_S2.R
# Purpose: S2 multimodal geometry DGP, migrated from the old S5 design
################################################################################

expit <- function(x) 1 / (1 + exp(-x))


softmax_rows <- function(logit_mat) {
    shifted <- logit_mat - apply(logit_mat, 1, max)
    exp_shifted <- exp(shifted)
    exp_shifted / rowSums(exp_shifted)
}


gen_X_S2 <- function(n) {
    X <- cbind(runif(n), runif(n))
    colnames(X) <- paste0("X", 1:2)
    X
}


mu0_S2 <- function(X) {
    X <- as.matrix(X)
    1.0 * sin(2 * pi * X[, 1]) + 0.7 * cos(pi * X[, 2]) +
        0.35 * sin(pi * X[, 1] * X[, 2])
}


tau_S2 <- function(X) {
    X <- as.matrix(X)
    0.7 + 0.6 * sin(pi * X[, 2]) + 0.5 * (X[, 1] - 0.5) -
        0.35 * cos(2 * pi * X[, 1]) + 0.25 * X[, 1] * X[, 2]
}


center_locations0_S2 <- function(X) {
    X <- as.matrix(X)
    anchor <- 0.20 * sin(2 * pi * X[, 1] * X[, 2]) - 0.15 * cos(pi * X[, 2])
    gap_left <- 0.65 + 1.75 * expit(
        1.8 * sin(2 * pi * X[, 1]) - 1.2 * cos(pi * X[, 2]) +
            0.8 * (X[, 1] - X[, 2])
    )
    gap_right <- 0.60 + 1.95 * expit(
        -1.6 * cos(2 * pi * X[, 2]) +
            1.1 * sin(pi * X[, 1] * X[, 2]) - 0.7 * (X[, 1] - 0.5)
    )
    middle_shift <- 0.25 * sin(pi * (X[, 1] + X[, 2]))
    cbind(anchor - gap_left, anchor + middle_shift, anchor + gap_right)
}


center_locations1_S2 <- function(X) {
    X <- as.matrix(X)
    anchor <- -0.10 * sin(pi * X[, 1]) + 0.20 * cos(2 * pi * X[, 1] * X[, 2])
    gap_left <- 0.55 + 2.05 * expit(
        -1.5 * sin(pi * X[, 2]) + 1.3 * cos(2 * pi * X[, 1]) +
            0.6 * (X[, 2] - 0.5)
    )
    gap_right <- 0.75 + 1.70 * expit(
        1.7 * sin(2 * pi * X[, 2]) - 0.9 * cos(pi * (X[, 1] - X[, 2])) +
            0.5 * (X[, 1] + X[, 2] - 1)
    )
    middle_shift <- -0.20 * cos(pi * (X[, 1] + 0.5 * X[, 2]))
    cbind(anchor - gap_left, anchor + middle_shift, anchor + gap_right)
}


mix_probs0_S2 <- function(X) {
    X <- as.matrix(X)
    logits <- cbind(
        1.1 * sin(2 * pi * X[, 1]) + 0.5 * cos(pi * X[, 2]),
        0.7 * cos(pi * (X[, 1] - X[, 2])) - 0.4 * (X[, 2] - 0.5),
        -2.3 + 4.2 * exp(-((X[, 1] - 0.25)^2 + (X[, 2] - 0.75)^2) / 0.05) +
            0.8 * sin(pi * X[, 2])
    )
    0.02 + (1 - 3 * 0.02) * softmax_rows(logits)
}


mix_probs1_S2 <- function(X) {
    X <- as.matrix(X)
    logits <- cbind(
        0.8 * cos(2 * pi * X[, 2]) - 0.4 * (X[, 1] - 0.5),
        0.8 * sin(pi * (X[, 1] + X[, 2])) + 0.3 * (X[, 2] - 0.5),
        -2.0 + 4.0 * exp(-((X[, 1] - 0.75)^2 + (X[, 2] - 0.30)^2) / 0.06) -
            0.7 * cos(pi * X[, 1])
    )
    0.02 + (1 - 3 * 0.02) * softmax_rows(logits)
}


component_sds0_S2 <- function(X) {
    X <- as.matrix(X)
    cbind(
        0.07 + 0.04 * X[, 1],
        0.10 + 0.05 * sin(pi * X[, 2])^2,
        0.08 + 0.07 * exp(-((X[, 1] - 0.25)^2 + (X[, 2] - 0.75)^2) / 0.08)
    )
}


component_sds1_S2 <- function(X) {
    X <- as.matrix(X)
    cbind(
        0.08 + 0.05 * X[, 2],
        0.09 + 0.04 * cos(pi * X[, 1])^2,
        0.07 + 0.08 * exp(-((X[, 1] - 0.75)^2 + (X[, 2] - 0.30)^2) / 0.09)
    )
}


sigma_S2 <- function(X) {
    X <- as.matrix(X)
    0.07 + 0.025 * X[, 1] + 0.025 * X[, 2]
}


ps_S2 <- function(X, ...) {
    X <- as.matrix(X)
    eta <- -0.25 + 1.0 * (X[, 1] - 0.5) - 0.9 * (X[, 2] - 0.5) +
        0.7 * sin(pi * X[, 1] * X[, 2])
    pmin(pmax(expit(eta), 0.15), 0.85)
}


centered_mixture_stats_S2 <- function(probs, raw_centers, component_sds) {
    center_mean <- rowSums(probs * raw_centers)
    centered_locations <- raw_centers - center_mean
    variance <- rowSums(probs * (centered_locations^2 + component_sds^2))
    list(centered_locations = centered_locations, variance = variance)
}


draw_centered_mixture_offset_S2 <- function(probs, raw_centers, component_sds) {
    stats <- centered_mixture_stats_S2(probs, raw_centers, component_sds)
    cumulative_probs <- t(apply(probs, 1, cumsum))
    u <- runif(nrow(probs))
    component_id <- 1 + (u > cumulative_probs[, 1]) + (u > cumulative_probs[, 2])
    rows <- seq_len(nrow(probs))
    stats$centered_locations[cbind(rows, component_id)] +
        component_sds[cbind(rows, component_id)] * rnorm(nrow(probs))
}


generate_data_S2 <- function(ntrain, ntest, seed, rho = 0, ...) {
    set.seed(seed)
    stopifnot(rho >= -1, rho <= 1)

    n <- ntrain + ntest
    X <- gen_X_S2(n)
    mu0 <- mu0_S2(X)
    tau <- tau_S2(X)
    mu1 <- mu0 + tau
    sigma <- sigma_S2(X)
    ps <- ps_S2(X)

    centers0 <- center_locations0_S2(X)
    centers1 <- center_locations1_S2(X)
    probs0 <- mix_probs0_S2(X)
    probs1 <- mix_probs1_S2(X)
    comp_sds0 <- component_sds0_S2(X)
    comp_sds1 <- component_sds1_S2(X)

    offset0 <- draw_centered_mixture_offset_S2(probs0, centers0, comp_sds0)
    offset1 <- draw_centered_mixture_offset_S2(probs1, centers1, comp_sds1)

    z0 <- rnorm(n)
    z1 <- rnorm(n)
    eps0 <- z0
    eps1 <- rho * z0 + sqrt(1 - rho^2) * z1

    Y0 <- mu0 + offset0 + sigma * eps0
    Y1 <- mu1 + offset1 + sigma * eps1
    ite <- Y1 - Y0

    offset_var0 <- centered_mixture_stats_S2(probs0, centers0, comp_sds0)$variance
    offset_var1 <- centered_mixture_stats_S2(probs1, centers1, comp_sds1)$variance
    std <- sqrt(offset_var0 + offset_var1 + 2 * (1 - rho) * sigma^2)

    build_data_simul(X, Y0, Y1, tau, ite, std, ps, ntrain, ntest)
}


true_outcome_density_S2 <- function(X, y_grid, arm = 0, ...) {
    X <- as.matrix(X)
    if (nrow(X) != 1) {
        stop("true_outcome_density_S2 expects one row of X.")
    }

    if (arm == 0) {
        probs <- mix_probs0_S2(X)
        centers <- center_locations0_S2(X)
        comp_sds <- component_sds0_S2(X)
        base_mean <- mu0_S2(X)
    } else {
        probs <- mix_probs1_S2(X)
        centers <- center_locations1_S2(X)
        comp_sds <- component_sds1_S2(X)
        base_mean <- mu0_S2(X) + tau_S2(X)
    }

    stats <- centered_mixture_stats_S2(probs, centers, comp_sds)
    component_sd <- sqrt(comp_sds[1, ]^2 + sigma_S2(X)^2)
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


oracle_y_grid_S2 <- function(X, arm = 0, n_grid = 500,
                             grid_resolution = 0.05,
                             max_grid_size = 2000, ...) {
    X <- as.matrix(X)
    if (arm == 0) {
        centers <- center_locations0_S2(X)
        probs <- mix_probs0_S2(X)
        comp_sds <- component_sds0_S2(X)
        base_mean <- mu0_S2(X)
    } else {
        centers <- center_locations1_S2(X)
        probs <- mix_probs1_S2(X)
        comp_sds <- component_sds1_S2(X)
        base_mean <- mu0_S2(X) + tau_S2(X)
    }
    stats <- centered_mixture_stats_S2(probs, centers, comp_sds)
    component_sd <- sqrt(comp_sds[1, ]^2 + sigma_S2(X)^2)
    means <- as.numeric(base_mean + stats$centered_locations[1, ])
    lower <- min(means - 5 * component_sd)
    upper <- max(means + 5 * component_sd)
    n_grid_i <- min(max(n_grid, ceiling((upper - lower) / grid_resolution)), max_grid_size)
    seq(lower, upper, length.out = n_grid_i)
}

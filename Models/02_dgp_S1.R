################################################################################
# Module: 02_dgp_S1.R
# Purpose: S1 randomized univariate DGP, migrated from the old S2 design
################################################################################

mu0_S1 <- function(X) {
    x <- as.matrix(X)[, 1]
    sin(pi * x) + x^2 / 3 + exp(-x^2)
}


tau_S1 <- function(X) {
    x <- as.matrix(X)[, 1]
    1 + sin(2 * x) + 0.5 * x * exp(-x^2 / 2)
}


sigma_S1 <- function(X) {
    x <- as.matrix(X)[, 1]
    0.5 + 0.5 * abs(x)
}


ps_S1 <- function(X, ...) {
    rep(0.5, nrow(as.matrix(X)))
}


generate_data_S1 <- function(ntrain, ntest, seed, rho = 0, ...) {
    set.seed(seed)
    stopifnot(rho >= -1, rho <= 1)

    n <- ntrain + ntest
    X <- matrix(rnorm(n), ncol = 1)
    colnames(X) <- "X1"

    mu0 <- mu0_S1(X)
    tau <- tau_S1(X)
    mu1 <- mu0 + tau
    sigma <- sigma_S1(X)
    ps <- ps_S1(X)

    z0 <- rnorm(n)
    z1 <- rnorm(n)
    eps0 <- z0
    eps1 <- rho * z0 + sqrt(1 - rho^2) * z1

    Y0 <- mu0 + sigma * eps0
    Y1 <- mu1 + sigma * eps1
    ite <- Y1 - Y0
    std <- sqrt(2 * (1 - rho)) * sigma

    build_data_simul(X, Y0, Y1, tau, ite, std, ps, ntrain, ntest)
}


true_outcome_density_S1 <- function(X, y_grid, arm = 0, ...) {
    X <- as.matrix(X)
    if (nrow(X) != 1) {
        stop("true_outcome_density_S1 expects one row of X.")
    }
    mean_y <- mu0_S1(X) + ifelse(arm == 1, tau_S1(X), 0)
    stats::dnorm(y_grid, mean = mean_y, sd = sigma_S1(X))
}


oracle_y_grid_S1 <- function(X, arm = 0, n_grid = 500,
                             grid_resolution = 0.05,
                             max_grid_size = 2000, ...) {
    X <- as.matrix(X)
    mean_y <- mu0_S1(X) + ifelse(arm == 1, tau_S1(X), 0)
    sd_y <- sigma_S1(X)
    lower <- mean_y - 5 * sd_y
    upper <- mean_y + 5 * sd_y
    n_grid_i <- min(max(n_grid, ceiling((upper - lower) / grid_resolution)), max_grid_size)
    seq(lower, upper, length.out = n_grid_i)
}

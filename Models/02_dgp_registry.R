################################################################################
# Module: 02_dgp_registry.R
# Purpose: Shared DGP interface for data generation and oracle quantities
################################################################################

source("02_dgp_S1.R")
source("02_dgp_S2.R")
source("02_dgp_S3.R")
source("02_dgp_S3_propensity_misspec.R")
source("02_dgp_S4.R")


build_data_simul <- function(X, Y0, Y1, CATE, ITE, std, ps, ntrain, ntest) {
    trainid <- seq_len(ntrain)
    testid <- (ntrain + 1):(ntrain + ntest)

    Xtrain <- X[trainid, , drop = FALSE]
    Xtest <- X[testid, , drop = FALSE]
    pstrain <- ps[trainid]
    pstest <- ps[testid]

    Ttrain <- as.numeric(runif(ntrain) <= pstrain)
    Ttest <- as.numeric(runif(ntest) <= pstest)
    Ytrain <- ifelse(Ttrain == 1, Y1[trainid], Y0[trainid])
    Ytest <- ifelse(Ttest == 1, Y1[testid], Y0[testid])

    list(
        Xtrain = Xtrain,
        Ttrain = Ttrain,
        Ytrain = Ytrain,
        Xtest = Xtest,
        Ttest = Ttest,
        Ytest = Ytest,
        Y0test = Y0[testid],
        Y1test = Y1[testid],
        tautest = ITE[testid],
        CATEtest = CATE[testid],
        stdtest = std[testid],
        pstrain = pstrain,
        pstest = pstest
    )
}


generate_dataset <- function(scenario, ntrain, ntest, seed, rho = 0, ...) {
    if (scenario == "S4") return(generate_data_S4(ntrain, ntest, seed, rho, ...))
    if (scenario == "S1") {
        return(generate_data_S1(ntrain, ntest, seed, rho, ...))
    }
    if (scenario == "S2") {
        return(generate_data_S2(ntrain, ntest, seed, rho, ...))
    }
    if (scenario == "S3") {
        return(generate_data_S3(ntrain, ntest, seed, rho, ...))
    }
    stop(sprintf("Unknown or unsupported scenario: %s", scenario))
}


get_true_ps <- function(scenario, X, ...) {
    if (scenario == "S4") return(ps_S4(X, ...))
    if (scenario == "S1") return(ps_S1(X, ...))
    if (scenario == "S2") return(ps_S2(X, ...))
    if (scenario == "S3") return(ps_S3(X, ...))
    stop(sprintf("True propensity is not implemented for scenario %s.", scenario))
}


get_true_cate <- function(scenario, X, ...) {
    if (scenario == "S4") return(tau_S4(X, ...))
    if (scenario == "S1") return(tau_S1(X))
    if (scenario == "S2") return(tau_S2(X))
    if (scenario == "S3") return(tau_S3(X))
    stop(sprintf("True CATE is not implemented for scenario %s.", scenario))
}


get_true_outcome_density <- function(scenario, X, y_grid, arm = 0, ...) {
    if (scenario == "S4") return(true_outcome_density_S4(X, y_grid, arm, ...))
    if (scenario == "S1") return(true_outcome_density_S1(X, y_grid, arm, ...))
    if (scenario == "S2") return(true_outcome_density_S2(X, y_grid, arm, ...))
    if (scenario == "S3") return(true_outcome_density_S3(X, y_grid, arm, ...))
    stop(sprintf("Oracle density is not implemented for scenario %s.", scenario))
}


get_oracle_y_grid <- function(scenario, X, arm = 0, n_grid = 500,
                              grid_resolution = 0.05,
                              max_grid_size = 2000, ...) {
    if (scenario == "S4") return(oracle_y_grid_S4(X, arm, n_grid, grid_resolution, max_grid_size, ...))
    if (scenario == "S1") {
        return(oracle_y_grid_S1(X, arm, n_grid, grid_resolution, max_grid_size, ...))
    }
    if (scenario == "S2") {
        return(oracle_y_grid_S2(X, arm, n_grid, grid_resolution, max_grid_size, ...))
    }
    if (scenario == "S3") {
        return(oracle_y_grid_S3(X, arm, n_grid, grid_resolution, max_grid_size, ...))
    }
    stop(sprintf("Oracle grid is not implemented for scenario %s.", scenario))
}

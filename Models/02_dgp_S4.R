################################################################################
# Module: 02_dgp_S4.R
# Purpose: S4 NLSM-based semi-synthetic DGP with X-dependent bimodal outcomes
################################################################################

# S4 preserves the S4 NLSM empirical DGP object and sampling protocol, but
# replaces the standard normal residuals with X-dependent centered bimodal
# shifts plus small within-mode Gaussian noise. Conditional on X, both Y(0)
# and Y(1) are bimodal with zero-mean added components, while Y(1) - Y(0)
# has three separated modes induced by the difference of two binary shifts.


.S4_cache <- new.env(parent = emptyenv())


default_dgp_path_S4 <- function() {
  getOption(
    "S4.NLSM_dgp_path",
    file.path("..", "Data", "NLSM_simul.RData")
  )
}


standardized_scores_S4 <- function(X, reference_X = X) {
    X <- as.matrix(X)
    reference_X <- as.matrix(reference_X)

    if (ncol(X) != ncol(reference_X)) {
        stop("X and reference_X must have the same number of columns.")
    }

    center <- colMeans(reference_X)
    scale <- apply(reference_X, 2, stats::sd)
    scale[!is.finite(scale) | scale == 0] <- 1

    Z <- sweep(sweep(X, 2, center, "-"), 2, scale, "/")
    z_col <- function(j) {
        if (ncol(Z) >= j) Z[, j] else rep(0, nrow(Z))
    }

    list(
        score_q = 0.8 * z_col(1) - 0.6 * z_col(2) + 0.4 * z_col(4),
        score_h = 0.5 * z_col(3) + 0.5 * z_col(5) - 0.3 * z_col(6)
    )
}


bimodal_params_S4 <- function(X, IQR0, IQR1,
                                reference_X = X,
                                q_min = 0.20,
                                q_max = 0.80,
                                gap_multiplier = 2.60,
                                gap_x_strength = 0.35,
                                within_ratio = 0.15,
                                min_within_sd = 1e-6) {
    if (length(q_min) != 1L || length(q_max) != 1L ||
        !is.finite(q_min) || !is.finite(q_max) ||
        q_min <= 0 || q_max >= 1 || q_min >= q_max) {
        stop("q_min and q_max must satisfy 0 < q_min < q_max < 1.")
    }

    if (length(gap_multiplier) != 1L || !is.finite(gap_multiplier) ||
        gap_multiplier <= 0) {
        stop("gap_multiplier must be a positive finite scalar.")
    }

    if (length(gap_x_strength) != 1L || !is.finite(gap_x_strength)) {
        stop("gap_x_strength must be a finite scalar.")
    }

    if (length(within_ratio) != 1L || !is.finite(within_ratio) ||
        within_ratio <= 0) {
        stop("within_ratio must be a positive finite scalar.")
    }

    X <- as.matrix(X)
    IQR0 <- as.numeric(IQR0)
    IQR1 <- as.numeric(IQR1)

    if (length(IQR0) != nrow(X) || length(IQR1) != nrow(X)) {
        stop("IQR0 and IQR1 must have one value per row of X.")
    }

    scores <- standardized_scores_S4(X, reference_X = reference_X)
    q <- q_min + (q_max - q_min) * stats::plogis(scores$score_q)

    sigma0 <- 0.5 * IQR0
    sigma1 <- 0.5 * IQR1
    base_scale <- sqrt((sigma0^2 + sigma1^2) / 2)

    fallback_scale <- stats::median(base_scale[is.finite(base_scale) &
                                                   base_scale > 0])
    if (!is.finite(fallback_scale)) {
        fallback_scale <- 1
    }
    base_scale[!is.finite(base_scale) | base_scale <= 0] <- fallback_scale

    h <- gap_multiplier * base_scale *
        exp(gap_x_strength * tanh(scores$score_h))
    within_sd <- pmax(within_ratio * h, min_within_sd)

    denom <- sqrt(q * (1 - q))
    shift_lower <- -h * q / denom
    shift_upper <- h * (1 - q) / denom

    list(
        q = q,
        h = h,
        within_sd = within_sd,
        shift_lower = shift_lower,
        shift_upper = shift_upper
    )
}


draw_bimodal_shift_S4 <- function(q, h) {
    q <- as.numeric(q)
    h <- as.numeric(h)

    if (length(q) != length(h)) {
        stop("q and h must have the same length.")
    }
    if (any(!is.finite(q)) || any(!is.finite(h)) ||
        any(q <= 0) || any(q >= 1) || any(h < 0)) {
        stop("q must lie in (0, 1) and h must be nonnegative.")
    }

    B <- as.numeric(stats::runif(length(q)) <= q)
    h * (B - q) / sqrt(q * (1 - q))
}


validate_dgp_S4 <- function(data) {
    required <- c("X", "EY0", "Etau", "IQR0", "IQR1", "ps")
    missing <- setdiff(required, names(data))
    if (length(missing) > 0) {
        stop(
            "S4 DGP object is missing required fields: ",
            paste(missing, collapse = ", ")
        )
    }

    X <- as.matrix(data$X)
    n <- nrow(X)

    if (n < 1L || ncol(X) < 1L) {
        stop("S4 DGP field `X` must be a non-empty matrix.")
    }

    expected_length <- c(
        EY0 = length(data$EY0),
        Etau = length(data$Etau),
        IQR0 = length(data$IQR0),
        IQR1 = length(data$IQR1),
        ps = length(data$ps)
    )

    if (any(expected_length != n)) {
        stop(
            "S4 DGP fields must all have one value per row of X. ",
            "Observed nrow(X) = ", n, "; lengths: ",
            paste(
                paste0(names(expected_length), "=", expected_length),
                collapse = ", "
            )
        )
    }

    if (anyNA(X) ||
        anyNA(data$EY0) ||
        anyNA(data$Etau) ||
        anyNA(data$IQR0) ||
        anyNA(data$IQR1) ||
        anyNA(data$ps)) {
        stop("S4 DGP object contains missing values.")
    }

    if (any(!is.finite(X)) ||
        any(!is.finite(data$EY0)) ||
        any(!is.finite(data$Etau)) ||
        any(!is.finite(data$IQR0)) ||
        any(!is.finite(data$IQR1)) ||
        any(!is.finite(data$ps))) {
        stop("S4 DGP object contains non-finite values.")
    }

    if (any(data$IQR0 < 0) || any(data$IQR1 < 0)) {
        stop("S4 DGP IQR values must be nonnegative.")
    }

    if (any(data$ps < 0) || any(data$ps > 1)) {
        stop("S4 DGP propensity scores must lie in [0, 1].")
    }

    data$X <- X
    data$EY0 <- as.numeric(data$EY0)
    data$Etau <- as.numeric(data$Etau)
    data$IQR0 <- as.numeric(data$IQR0)
    data$IQR1 <- as.numeric(data$IQR1)
    data$ps <- as.numeric(data$ps)

    data
}


load_dgp_S4 <- function(dgp_path = default_dgp_path_S4(),
                          use_cache = TRUE) {
    if (!file.exists(dgp_path)) {
        stop(
            "Cannot find the S4 NLSM DGP file: ", dgp_path, "\n",
            "Set options(S4.NLSM_dgp_path = '.../NLSM_simul.RData') ",
            "or pass dgp_path explicitly."
        )
    }

    cache_key <- normalizePath(dgp_path, winslash = "/", mustWork = TRUE)

    if (use_cache && exists(cache_key, envir = .S4_cache, inherits = FALSE)) {
        return(get(cache_key, envir = .S4_cache, inherits = FALSE))
    }

    env <- new.env(parent = emptyenv())
    loaded_names <- load(dgp_path, envir = env)

    if (!("data" %in% loaded_names)) {
        stop(
            "The S4 DGP file must contain an object named `data`. ",
            "Objects found: ", paste(loaded_names, collapse = ", ")
        )
    }

    data <- validate_dgp_S4(env$data)

    if (use_cache) {
        assign(cache_key, data, envir = .S4_cache)
    }

    data
}


clear_dgp_cache_S4 <- function() {
    rm(list = ls(envir = .S4_cache, all.names = TRUE), envir = .S4_cache)
    invisible(NULL)
}


generate_data_S4 <- function(ntrain, ntest, seed, rho = 0,
                               dgp_path = default_dgp_path_S4(),
                               q_min = 0.20,
                               q_max = 0.80,
                               gap_multiplier = 2.60,
                               gap_x_strength = 0.35,
                               within_ratio = 0.15,
                               ...) {
    if (length(rho) != 1L || !is.finite(rho) || rho != 0) {
        stop(
            "S4 follows the S4 NLSM sampling protocol and requires rho = 0."
        )
    }

    if (length(ntrain) != 1L || length(ntest) != 1L ||
        ntrain < 1 || ntest < 1 ||
        ntrain != as.integer(ntrain) ||
        ntest != as.integer(ntest)) {
        stop("ntrain and ntest must be positive integers.")
    }

    data <- load_dgp_S4(dgp_path)

    X <- data$X
    EY0 <- data$EY0
    Etau <- data$Etau
    IQR0 <- data$IQR0
    IQR1 <- data$IQR1

    # Retain the S4 cap exactly.
    ps <- pmin(data$ps, 1)

    n <- nrow(X)
    n_required <- ntrain + ntest

    if (n_required > n) {
        stop(
            "S4 samples without replacement from the fixed NLSM pool. ",
            "Requested ntrain + ntest = ", n_required,
            ", but the available pool has only ", n, " observations."
        )
    }

    # Keep the S4 seed convention.
    set.seed(seed * 2020)

    params <- bimodal_params_S4(
        X = X,
        IQR0 = IQR0,
        IQR1 = IQR1,
        reference_X = X,
        q_min = q_min,
        q_max = q_max,
        gap_multiplier = gap_multiplier,
        gap_x_strength = gap_x_strength,
        within_ratio = within_ratio
    )

    shift0 <- draw_bimodal_shift_S4(params$q, params$h)
    shift1 <- draw_bimodal_shift_S4(params$q, params$h)
    noise0 <- params$within_sd * stats::rnorm(n)
    noise1 <- params$within_sd * stats::rnorm(n)

    Y0 <- EY0 + shift0 + noise0
    Y1 <- EY0 + Etau + shift1 + noise1

    ITE <- Y1 - Y0
    std <- sqrt(2 * params$h^2 + 2 * params$within_sd^2)

    # Keep the S4 empirical-pool permutation protocol.
    ids <- sample(n, n, replace = FALSE)
    keep <- ids[seq_len(n_required)]

    X <- X[keep, , drop = FALSE]
    Y0 <- Y0[keep]
    Y1 <- Y1[keep]
    CATE <- Etau[keep]
    ITE <- ITE[keep]
    std <- std[keep]
    ps <- ps[keep]

    build_data_simul(
        X = X,
        Y0 = Y0,
        Y1 = Y1,
        CATE = CATE,
        ITE = ITE,
        std = std,
        ps = ps,
        ntrain = ntrain,
        ntest = ntest
    )
}


################################################################################
# Oracle quantities
################################################################################

row_key_S4 <- function(X) {
    X <- as.matrix(X)

    apply(
        X,
        1,
        function(x) {
            paste(sprintf("%.17g", as.numeric(x)), collapse = "\034")
        }
    )
}


lookup_dgp_rows_S4 <- function(X,
                                 dgp_path = default_dgp_path_S4()) {
    data <- load_dgp_S4(dgp_path)
    X <- as.matrix(X)

    if (ncol(X) != ncol(data$X)) {
        stop(
            "S4 X has ", ncol(X), " columns, but the NLSM DGP has ",
            ncol(data$X), " columns."
        )
    }

    dgp_keys <- attr(data, "row_keys", exact = TRUE)
    if (is.null(dgp_keys)) {
        dgp_keys <- row_key_S4(data$X)
    }

    query_keys <- row_key_S4(X)
    idx <- match(query_keys, dgp_keys)

    if (anyNA(idx)) {
        stop(
            "Some rows of X are not present in the fixed S4 NLSM covariate pool. ",
            "S4 oracle functions are defined by lookup on that empirical pool."
        )
    }

    idx
}


ps_S4 <- function(X, dgp_path = default_dgp_path_S4(), ...) {
    data <- load_dgp_S4(dgp_path)
    idx <- lookup_dgp_rows_S4(X, dgp_path)
    pmin(data$ps[idx], 1)
}


tau_S4 <- function(X, dgp_path = default_dgp_path_S4(), ...) {
    data <- load_dgp_S4(dgp_path)
    idx <- lookup_dgp_rows_S4(X, dgp_path)
    data$Etau[idx]
}


true_ite_density_S4 <- function(X, tau_grid,
                                  dgp_path = default_dgp_path_S4(),
                                  q_min = 0.20,
                                  q_max = 0.80,
                                  gap_multiplier = 2.60,
                                  gap_x_strength = 0.35,
                                  within_ratio = 0.15,
                                  ...) {
    X <- as.matrix(X)

    if (nrow(X) != 1L) {
        stop("true_ite_density_S4 expects one row of X.")
    }

    data <- load_dgp_S4(dgp_path)
    idx <- lookup_dgp_rows_S4(X, dgp_path)
    params <- bimodal_params_S4(
        X = data$X,
        IQR0 = data$IQR0,
        IQR1 = data$IQR1,
        reference_X = data$X,
        q_min = q_min,
        q_max = q_max,
        gap_multiplier = gap_multiplier,
        gap_x_strength = gap_x_strength,
        within_ratio = within_ratio
    )

    q <- params$q[idx]
    gap <- params$h[idx] / sqrt(q * (1 - q))
    component_means <- data$Etau[idx] + c(-gap, 0, gap)
    component_probs <- c(
        q * (1 - q),
        q^2 + (1 - q)^2,
        q * (1 - q)
    )
    component_sd <- sqrt(2) * params$within_sd[idx]

    density <- rep(0, length(tau_grid))
    for (k in seq_along(component_means)) {
        density <- density + component_probs[k] * stats::dnorm(
            tau_grid,
            mean = component_means[k],
            sd = component_sd
        )
    }
    density
}


oracle_tau_grid_S4 <- function(X, n_grid = 500,
                                 grid_resolution = 0.05,
                                 max_grid_size = 2000,
                                 dgp_path = default_dgp_path_S4(),
                                 q_min = 0.20,
                                 q_max = 0.80,
                                 gap_multiplier = 2.60,
                                 gap_x_strength = 0.35,
                                 within_ratio = 0.15,
                                 ...) {
    X <- as.matrix(X)

    if (nrow(X) != 1L) {
        stop("oracle_tau_grid_S4 expects one row of X.")
    }

    data <- load_dgp_S4(dgp_path)
    idx <- lookup_dgp_rows_S4(X, dgp_path)
    params <- bimodal_params_S4(
        X = data$X,
        IQR0 = data$IQR0,
        IQR1 = data$IQR1,
        reference_X = data$X,
        q_min = q_min,
        q_max = q_max,
        gap_multiplier = gap_multiplier,
        gap_x_strength = gap_x_strength,
        within_ratio = within_ratio
    )

    q <- params$q[idx]
    gap <- params$h[idx] / sqrt(q * (1 - q))
    component_means <- data$Etau[idx] + c(-gap, 0, gap)
    component_sd <- sqrt(2) * params$within_sd[idx]

    lower <- min(component_means - 5 * component_sd)
    upper <- max(component_means + 5 * component_sd)

    n_grid_i <- min(
        max(
            n_grid,
            ceiling((upper - lower) / grid_resolution)
        ),
        max_grid_size
    )

    seq(lower, upper, length.out = n_grid_i)
}


true_outcome_density_S4 <- function(X, y_grid, arm = 0,
                                      dgp_path = default_dgp_path_S4(),
                                      q_min = 0.20,
                                      q_max = 0.80,
                                      gap_multiplier = 2.60,
                                      gap_x_strength = 0.35,
                                      within_ratio = 0.15,
                                      ...) {
    X <- as.matrix(X)

    if (nrow(X) != 1L) {
        stop("true_outcome_density_S4 expects one row of X.")
    }

    if (!(arm %in% c(0, 1))) {
        stop("arm must be 0 or 1.")
    }

    data <- load_dgp_S4(dgp_path)
    idx <- lookup_dgp_rows_S4(X, dgp_path)
    params <- bimodal_params_S4(
        X = data$X,
        IQR0 = data$IQR0,
        IQR1 = data$IQR1,
        reference_X = data$X,
        q_min = q_min,
        q_max = q_max,
        gap_multiplier = gap_multiplier,
        gap_x_strength = gap_x_strength,
        within_ratio = within_ratio
    )

    if (arm == 0) {
        mean_y <- data$EY0[idx]
    } else {
        mean_y <- data$EY0[idx] + data$Etau[idx]
    }

    component_means <- mean_y + c(
        params$shift_lower[idx],
        params$shift_upper[idx]
    )
    component_probs <- c(1 - params$q[idx], params$q[idx])
    component_sd <- params$within_sd[idx]

    density <- rep(0, length(y_grid))
    for (k in seq_along(component_means)) {
        density <- density + component_probs[k] * stats::dnorm(
            y_grid,
            mean = component_means[k],
            sd = component_sd
        )
    }
    density
}


oracle_y_grid_S4 <- function(X, arm = 0, n_grid = 500,
                               grid_resolution = 0.05,
                               max_grid_size = 2000,
                               dgp_path = default_dgp_path_S4(),
                               q_min = 0.20,
                               q_max = 0.80,
                               gap_multiplier = 2.60,
                               gap_x_strength = 0.35,
                               within_ratio = 0.15,
                               ...) {
    X <- as.matrix(X)

    if (nrow(X) != 1L) {
        stop("oracle_y_grid_S4 expects one row of X.")
    }

    if (!(arm %in% c(0, 1))) {
        stop("arm must be 0 or 1.")
    }

    data <- load_dgp_S4(dgp_path)
    idx <- lookup_dgp_rows_S4(X, dgp_path)
    params <- bimodal_params_S4(
        X = data$X,
        IQR0 = data$IQR0,
        IQR1 = data$IQR1,
        reference_X = data$X,
        q_min = q_min,
        q_max = q_max,
        gap_multiplier = gap_multiplier,
        gap_x_strength = gap_x_strength,
        within_ratio = within_ratio
    )

    if (arm == 0) {
        mean_y <- data$EY0[idx]
    } else {
        mean_y <- data$EY0[idx] + data$Etau[idx]
    }

    component_means <- mean_y + c(
        params$shift_lower[idx],
        params$shift_upper[idx]
    )
    component_sd <- params$within_sd[idx]

    lower <- min(component_means - 5 * component_sd)
    upper <- max(component_means + 5 * component_sd)

    n_grid_i <- min(
        max(
            n_grid,
            ceiling((upper - lower) / grid_resolution)
        ),
        max_grid_size
    )

    seq(lower, upper, length.out = n_grid_i)
}

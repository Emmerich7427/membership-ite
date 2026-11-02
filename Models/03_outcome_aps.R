################################################################################
# Module: 03_outcome_aps.R
# Purpose: Outcome-level APS prediction sets with fitted or oracle densities
################################################################################

source("01_interval_utils.R")
source("01_conformal_utils.R")
source("02_models.R")
source("02_dgp_registry.R")


logsumexp <- function(x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
}


scale_joint_data <- function(X, Y) {
    Z <- cbind(as.matrix(X), Y = as.numeric(Y))
    center <- colMeans(Z)
    scale <- apply(Z, 2, stats::sd)
    scale[scale == 0] <- 1
    list(Z = sweep(sweep(Z, 2, center, "-"), 2, scale, "/"), center = center, scale = scale)
}


niw_posterior <- function(Z, prior) {
    d <- ncol(prior$Lambda0)
    n <- if (is.null(Z)) 0 else nrow(Z)
    if (n == 0) {
        return(list(kappa = prior$kappa0, nu = prior$nu0, m = prior$m0, Lambda = prior$Lambda0))
    }

    zbar <- colMeans(Z)
    centered <- sweep(Z, 2, zbar, "-")
    S <- crossprod(centered)
    kappa_n <- prior$kappa0 + n
    nu_n <- prior$nu0 + n
    mean_diff <- matrix(zbar - prior$m0, ncol = 1)
    Lambda_n <- prior$Lambda0 + S +
        (prior$kappa0 * n / kappa_n) * tcrossprod(mean_diff)
    m_n <- (prior$kappa0 * prior$m0 + n * zbar) / kappa_n
    list(kappa = kappa_n, nu = nu_n, m = as.numeric(m_n), Lambda = Lambda_n)
}


log_mvt_density <- function(x, mean, scale, df, jitter = 1e-8) {
    x <- as.numeric(x)
    p <- length(x)
    scale <- as.matrix(scale) + diag(jitter, p)
    R <- chol(scale)
    logdet <- 2 * sum(log(diag(R)))
    delta <- sum(backsolve(R, x - mean, transpose = TRUE)^2)
    lgamma((df + p) / 2) - lgamma(df / 2) -
        0.5 * (p * log(df * pi) + logdet) -
        0.5 * (df + p) * log1p(delta / df)
}


log_niw_predictive <- function(z, Z_cluster, prior, jitter = 1e-8) {
    post <- niw_posterior(Z_cluster, prior)
    d <- length(z)
    df <- post$nu - d + 1
    scale <- ((post$kappa + 1) / (post$kappa * df)) * post$Lambda
    log_mvt_density(z, post$m, scale, df, jitter = jitter)
}


sample_inv_wishart <- function(nu, Lambda) {
    W <- stats::rWishart(1, df = nu, Sigma = solve(Lambda))[, , 1]
    solve(W)
}


sample_cluster_params <- function(Z_cluster, prior, jitter = 1e-8) {
    post <- niw_posterior(Z_cluster, prior)
    d <- length(post$m)
    Sigma <- sample_inv_wishart(post$nu, post$Lambda + diag(jitter, d))
    mu <- as.numeric(MASS::mvrnorm(1, mu = post$m, Sigma = Sigma / post$kappa))
    list(mu = mu, Sigma = Sigma)
}


compact_labels <- function(z) {
    active <- sort(unique(z[!is.na(z)]))
    z[!is.na(z)] <- match(z[!is.na(z)], active)
    z
}


initialize_labels <- function(Z, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    n <- nrow(Z)
    k0 <- min(5, max(1, floor(sqrt(n / 4))))
    if (k0 == 1) {
        return(rep(1L, n))
    }
    as.integer(stats::kmeans(Z, centers = k0, nstart = 10)$cluster)
}


fit_joint_dp_mixture_density <- function(X, Y,
                                         n_fit_max = 10000,
                                         n_iter = 250,
                                         burn = 125,
                                         thin = 5,
                                         alpha = 5,
                                         kappa0 = 0.05,
                                         nu0_extra = 3,
                                         prior_scale = 1,
                                         jitter = 1e-7,
                                         seed = NULL,
                                         verbose = TRUE) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package 'MASS' is required for the joint DP mixture.")
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }

    n_raw <- length(Y)
    if (!is.null(n_fit_max) && n_raw > n_fit_max) {
        fit_idx <- sample(seq_len(n_raw), n_fit_max)
        X <- X[fit_idx, , drop = FALSE]
        Y <- Y[fit_idx]
    } else {
        fit_idx <- seq_len(n_raw)
    }

    scaled <- scale_joint_data(X, Y)
    Z <- scaled$Z
    n <- nrow(Z)
    d <- ncol(Z)
    prior <- list(
        m0 = rep(0, d),
        kappa0 = kappa0,
        nu0 = d + nu0_extra,
        Lambda0 = diag(prior_scale^2, d)
    )

    z <- compact_labels(initialize_labels(Z, seed = if (!is.null(seed)) seed + 7 else NULL))
    keep_iters <- seq.int(burn + thin, n_iter, by = thin)
    draws <- vector("list", length(keep_iters))
    k_trace <- integer(n_iter)
    keep_pos <- 0

    if (verbose) {
        cat(sprintf("  [DP] n=%d/%d, d=%d, alpha=%.3f\n", n, n_raw, d, alpha))
    }

    for (iter in seq_len(n_iter)) {
        for (i in sample(seq_len(n))) {
            z[i] <- NA_integer_
            z <- compact_labels(z)
            K <- if (all(is.na(z))) 0 else max(z, na.rm = TRUE)
            log_probs <- numeric(K + 1)
            if (K > 0) {
                counts <- tabulate(z[!is.na(z)], nbins = K)
                for (k in seq_len(K)) {
                    log_probs[k] <- log(counts[k]) +
                        log_niw_predictive(Z[i, ], Z[z == k & !is.na(z), , drop = FALSE], prior, jitter)
                }
            }
            log_probs[K + 1] <- log(alpha) + log_niw_predictive(Z[i, ], NULL, prior, jitter)
            z[i] <- sample(seq_len(K + 1), 1, prob = exp(log_probs - logsumexp(log_probs)))
        }

        z <- compact_labels(z)
        k_trace[iter] <- max(z)
        if (iter %in% keep_iters) {
            keep_pos <- keep_pos + 1
            K <- max(z)
            counts <- as.integer(tabulate(z, nbins = K))
            params <- vector("list", K)
            for (k in seq_len(K)) {
                params[[k]] <- sample_cluster_params(Z[z == k, , drop = FALSE], prior, jitter)
            }
            draws[[keep_pos]] <- list(counts = counts, params = params)
        }

        if (verbose && (iter %% 10 == 0 || iter == 1)) {
            cat(sprintf("    DP iter %d/%d, K=%d, kept=%d\n", iter, n_iter, k_trace[iter], keep_pos))
        }
    }

    structure(
        list(
            draws = draws,
            k_trace = k_trace,
            center = scaled$center,
            scale = scaled$scale,
            n_total = n_raw,
            n_fit = n,
            fit_idx = fit_idx,
            params = list(
                n_fit_max = n_fit_max, n_iter = n_iter, burn = burn, thin = thin,
                alpha = alpha, kappa0 = kappa0, nu0_extra = nu0_extra,
                prior_scale = prior_scale, jitter = jitter
            )
        ),
        class = "joint_dp_mixture_density"
    )
}


dmvnorm_x <- function(x, mean, Sigma, jitter = 1e-8) {
    p <- length(x)
    Sigma <- as.matrix(Sigma) + diag(jitter, p)
    R <- chol(Sigma)
    logdet <- 2 * sum(log(diag(R)))
    delta <- sum(backsolve(R, x - mean, transpose = TRUE)^2)
    exp(-0.5 * (p * log(2 * pi) + logdet + delta))
}


conditional_normal_params <- function(mu, Sigma, x_scaled, jitter = 1e-8) {
    p <- length(x_scaled)
    S_xx <- Sigma[seq_len(p), seq_len(p), drop = FALSE] + diag(jitter, p)
    S_xy <- Sigma[seq_len(p), p + 1, drop = FALSE]
    S_yx <- Sigma[p + 1, seq_len(p), drop = FALSE]
    S_yy <- Sigma[p + 1, p + 1]
    solved <- solve(S_xx, x_scaled - mu[seq_len(p)])
    mean_y <- as.numeric(mu[p + 1] + S_yx %*% solved)
    var_y <- as.numeric(S_yy - S_yx %*% solve(S_xx, S_xy))
    list(mean = mean_y, sd = sqrt(max(var_y, jitter)))
}


predict_joint_dp_mixture_density <- function(model, X_new, y_grid, jitter = NULL) {
    if (is.null(jitter)) {
        jitter <- model$params$jitter
    }

    X_new <- as.matrix(X_new)
    p <- ncol(X_new)
    x_scaled <- as.numeric((X_new[1, ] - model$center[seq_len(p)]) / model$scale[seq_len(p)])
    y_grid_scaled <- as.numeric((y_grid - model$center[p + 1]) / model$scale[p + 1])
    y_scale <- model$scale[p + 1]

    density_draws <- matrix(NA_real_, nrow = length(model$draws), ncol = length(y_grid))
    for (s in seq_along(model$draws)) {
        draw <- model$draws[[s]]
        K <- length(draw$params)
        x_weights <- numeric(K)
        comp_density <- matrix(0, nrow = K, ncol = length(y_grid))
        for (k in seq_len(K)) {
            param <- draw$params[[k]]
            mu <- param$mu
            Sigma <- param$Sigma
            x_weights[k] <- draw$counts[k] *
                dmvnorm_x(x_scaled, mu[seq_len(p)], Sigma[seq_len(p), seq_len(p), drop = FALSE], jitter)
            cond <- conditional_normal_params(mu, Sigma, x_scaled, jitter)
            comp_density[k, ] <- stats::dnorm(y_grid_scaled, mean = cond$mean, sd = cond$sd) / y_scale
        }
        if (sum(x_weights) <= 0 || any(!is.finite(x_weights))) {
            x_weights <- rep(1 / K, K)
        } else {
            x_weights <- x_weights / sum(x_weights)
        }
        density_draws[s, ] <- as.numeric(x_weights %*% comp_density)
    }

    list(density_mean = colMeans(density_draws), density_draws = density_draws)
}


compute_aps_score <- function(density_grid, y_grid, Y_obs) {
    delta <- mean(diff(y_grid))
    f_y <- stats::approx(y_grid, density_grid, xout = Y_obs, rule = 2)$y
    min(sum(density_grid[density_grid > f_y]) * delta, 1)
}


compute_aps_scores_grid <- function(density_grid, y_grid) {
    delta <- mean(diff(y_grid))
    ord <- order(density_grid, decreasing = TRUE)
    scores <- numeric(length(density_grid))
    cumulative_mass <- 0
    last_density <- Inf
    last_score <- 0
    for (idx in ord) {
        if (density_grid[idx] < last_density) {
            last_score <- cumulative_mass
            last_density <- density_grid[idx]
        }
        scores[idx] <- last_score
        cumulative_mass <- cumulative_mass + density_grid[idx] * delta
    }
    pmin(scores, 1)
}


aps_density_threshold <- function(density_grid,
                                  scores_grid,
                                  cutoff,
                                  tol = 1e-12) {
    density_grid <- as.numeric(density_grid)
    scores_grid <- as.numeric(scores_grid)
    cutoff <- as.numeric(cutoff)[1]

    if (length(density_grid) != length(scores_grid)) {
        stop("density_grid and scores_grid must have the same length.")
    }
    if (length(density_grid) == 0) {
        stop("density_grid and scores_grid must be nonempty.")
    }
    if (any(is.na(density_grid)) || any(is.na(scores_grid))) {
        stop("density_grid and scores_grid must not contain NA or NaN.")
    }
    if (is.na(cutoff)) {
        stop("cutoff must not be NA or NaN.")
    }
    if (!is.finite(tol) || tol < 0) {
        stop("tol must be a nonnegative finite number.")
    }
    if (cutoff == Inf || cutoff >= 1) {
        return(-Inf)
    }

    accepted <- scores_grid <= cutoff + tol
    if (any(accepted)) {
        return(min(density_grid[accepted]))
    }

    warning("No APS grid nodes accepted at this cutoff; using maximum density as fallback.")
    max(density_grid)
}


normalize_outcome_grid_params <- function(n_grid = 500,
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


resolve_outcome_prediction_grid_params <- function(fit_obj,
                                                   n_grid = NULL,
                                                   grid_resolution = NULL,
                                                   max_grid_size = NULL) {
    default_params <- normalize_outcome_grid_params()
    fit_params <- if (is.null(fit_obj$grid_params)) {
        default_params
    } else {
        normalize_outcome_grid_params(
            fit_obj$grid_params$n_grid,
            fit_obj$grid_params$grid_resolution,
            fit_obj$grid_params$max_grid_size
        )
    }

    resolved <- fit_params
    supplied <- list(
        n_grid = n_grid,
        grid_resolution = grid_resolution,
        max_grid_size = max_grid_size
    )

    for (name in names(supplied)) {
        if (!is.null(supplied[[name]])) {
            candidate <- resolved
            candidate[[name]] <- supplied[[name]]
            candidate <- normalize_outcome_grid_params(
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

    resolved
}


dp_y_grid_for_x <- function(model, X_new, n_grid = 500,
                            grid_resolution = 0.05,
                            max_grid_size = 2000) {
    p <- ncol(as.matrix(X_new))
    x_scaled <- as.numeric((as.matrix(X_new)[1, ] - model$center[seq_len(p)]) /
        model$scale[seq_len(p)])
    means <- c()
    sds <- c()
    for (draw in model$draws) {
        for (param in draw$params) {
            cond <- conditional_normal_params(param$mu, param$Sigma, x_scaled, model$params$jitter)
            means <- c(means, cond$mean * model$scale[p + 1] + model$center[p + 1])
            sds <- c(sds, cond$sd * model$scale[p + 1])
        }
    }
    lower <- min(means - 4 * sds, na.rm = TRUE)
    upper <- max(means + 4 * sds, na.rm = TRUE)
    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
        lower <- model$center[p + 1] - 4 * model$scale[p + 1]
        upper <- model$center[p + 1] + 4 * model$scale[p + 1]
    }
    n_grid_i <- min(max(n_grid, ceiling((upper - lower) / grid_resolution)), max_grid_size)
    seq(lower, upper, length.out = n_grid_i)
}


fit_bart_predictive_density <- function(X, Y, seed = NULL, verbose = TRUE) {
    if (!requireNamespace("dbarts", quietly = TRUE)) {
        stop("Package 'dbarts' is required for BART Stage 1.")
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }

    X_train <- as.data.frame(X)
    X_names <- paste0("X", seq_len(ncol(X_train)))
    colnames(X_train) <- X_names

    if (verbose) {
        cat(sprintf("  [BART] n=%d, p=%d\n", length(Y), ncol(X_train)))
    }

    bart_fit <- dbarts::bart(
        x.train = X_train,
        y.train = Y,
        keeptrees = TRUE,
        verbose = FALSE
    )

    structure(
        list(
            bart_fit = bart_fit,
            X_names = X_names,
            y_center = mean(Y),
            y_scale = stats::sd(Y),
            n_fit = length(Y)
        ),
        class = "bart_predictive_density"
    )
}


predict_bart_ppd <- function(model, X_new) {
    X_new <- as.data.frame(X_new)
    colnames(X_new) <- model$X_names
    ppd <- predict(model$bart_fit, newdata = X_new, type = "ppd")
    if (is.null(dim(ppd))) {
        ppd <- matrix(ppd, ncol = nrow(X_new))
    }
    ppd
}


bart_density_bandwidth <- function(draws, fallback_scale) {
    draws <- draws[is.finite(draws)]
    bw <- if (length(unique(draws)) > 1) stats::bw.nrd0(draws) else NA_real_
    if (!is.finite(bw) || bw <= 0) {
        bw <- stats::sd(draws)
    }
    if (!is.finite(bw) || bw <= 0) {
        bw <- fallback_scale
    }
    if (!is.finite(bw) || bw <= 0) {
        bw <- 1
    }
    bw
}


bart_y_grid_for_x <- function(model, X_new, n_grid = 500,
                              grid_resolution = 0.05,
                              max_grid_size = 2000) {
    draws <- as.numeric(predict_bart_ppd(model, X_new)[, 1])
    draws <- draws[is.finite(draws)]
    if (length(draws) == 0) {
        draws <- model$y_center
    }
    bw <- bart_density_bandwidth(draws, model$y_scale)
    lower <- min(draws, na.rm = TRUE) - 4 * bw
    upper <- max(draws, na.rm = TRUE) + 4 * bw
    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
        fallback_scale <- max(c(model$y_scale, bw, 1), na.rm = TRUE)
        lower <- model$y_center - 4 * fallback_scale
        upper <- model$y_center + 4 * fallback_scale
    }
    n_grid_i <- min(max(n_grid, ceiling((upper - lower) / grid_resolution)), max_grid_size)
    y_grid <- seq(lower, upper, length.out = n_grid_i)
    attr(y_grid, "bart_ppd_draws") <- draws
    y_grid
}


predict_bart_density <- function(model, X_new, y_grid) {
    draws <- attr(y_grid, "bart_ppd_draws", exact = TRUE)
    if (is.null(draws)) {
        ppd <- predict_bart_ppd(model, X_new)
        draws <- as.numeric(ppd[, 1])
    }
    draws <- draws[is.finite(draws)]
    if (length(draws) == 0) {
        draws <- model$y_center
    }
    bw <- bart_density_bandwidth(draws, model$y_scale)
    density_grid <- colMeans(stats::dnorm(outer(draws, y_grid, "-") / bw)) / bw
    list(density_mean = density_grid, density_draws = NULL)
}


normalize_mdn_hidden <- function(hidden) {
    if (is.null(hidden)) {
        return(integer(0))
    }
    if (is.character(hidden)) {
        hidden <- unlist(strsplit(hidden, ",", fixed = TRUE))
        hidden <- trimws(hidden)
        hidden <- hidden[nzchar(hidden)]
        if (length(hidden) == 0) {
            return(integer(0))
        }
    }
    hidden <- as.integer(hidden)
    if (any(is.na(hidden)) || any(hidden < 1L)) {
        stop("mdn_hidden must contain positive hidden-layer sizes.")
    }
    hidden
}


validate_mdn_params <- function(n_components,
                                hidden,
                                epochs,
                                batch_size,
                                lr,
                                weight_decay,
                                patience,
                                val_prop,
                                min_sigma,
                                grid_tail_prob,
                                grid_sd_mult) {
    if (length(n_components) != 1L || is.na(n_components) || n_components < 2L) {
        stop("mdn_n_components must be at least 2.")
    }
    if (length(epochs) != 1L || is.na(epochs) || epochs < 1L) {
        stop("mdn_epochs must be a positive integer.")
    }
    if (length(batch_size) != 1L || is.na(batch_size) || batch_size < 1L) {
        stop("mdn_batch_size must be a positive integer.")
    }
    if (length(lr) != 1L || is.na(lr) || lr <= 0) {
        stop("mdn_lr must be positive.")
    }
    if (length(weight_decay) != 1L || is.na(weight_decay) || weight_decay < 0) {
        stop("mdn_weight_decay must be nonnegative.")
    }
    if (length(patience) != 1L || is.na(patience) || patience < 1L) {
        stop("mdn_patience must be a positive integer.")
    }
    if (length(val_prop) != 1L || is.na(val_prop) || val_prop < 0 || val_prop >= 0.5) {
        stop("mdn_val_prop must be in [0, 0.5).")
    }
    if (length(min_sigma) != 1L || is.na(min_sigma) || min_sigma <= 0) {
        stop("mdn_min_sigma must be positive.")
    }
    if (length(grid_tail_prob) != 1L || is.na(grid_tail_prob) ||
        grid_tail_prob <= 0 || grid_tail_prob >= 0.1) {
        stop("mdn_grid_tail_prob must be in (0, 0.1).")
    }
    if (length(grid_sd_mult) != 1L || is.na(grid_sd_mult) || grid_sd_mult <= 0) {
        stop("mdn_grid_sd_mult must be positive.")
    }

    list(
        n_components = as.integer(n_components),
        hidden = normalize_mdn_hidden(hidden),
        epochs = as.integer(epochs),
        batch_size = as.integer(batch_size),
        lr = as.numeric(lr),
        weight_decay = as.numeric(weight_decay),
        patience = as.integer(patience),
        val_prop = as.numeric(val_prop),
        min_sigma = as.numeric(min_sigma),
        grid_tail_prob = as.numeric(grid_tail_prob),
        grid_sd_mult = as.numeric(grid_sd_mult)
    )
}


build_mdn_torch_module <- function() {
    torch::nn_module(
        "MDNDensityNet",
        initialize = function(input_dim, hidden_units, n_components,
                              min_sigma_scaled) {
            self$n_components <- as.integer(n_components)
            self$min_sigma_scaled <- as.numeric(min_sigma_scaled)

            layers <- list()
            last_dim <- as.integer(input_dim)
            if (length(hidden_units) > 0L) {
                for (h in hidden_units) {
                    layers[[length(layers) + 1L]] <- torch::nn_linear(last_dim, h)
                    layers[[length(layers) + 1L]] <- torch::nn_relu()
                    last_dim <- h
                }
            }
            self$net <- if (length(layers) > 0L) {
                do.call(torch::nn_sequential, layers)
            } else {
                NULL
            }
            self$head <- torch::nn_linear(last_dim, 3L * self$n_components)
        },
        forward = function(x) {
            h <- if (is.null(self$net)) x else self$net(x)
            raw <- self$head(h)
            k <- self$n_components

            logits <- raw[, seq_len(k), drop = FALSE]
            mu <- raw[, k + seq_len(k), drop = FALSE]
            sigma_raw <- raw[, 2L * k + seq_len(k), drop = FALSE]
            sigma <- torch::nnf_softplus(sigma_raw) + self$min_sigma_scaled

            list(logits = logits, mu = mu, sigma = sigma)
        }
    )
}


mdn_nll <- function(logits, mu, sigma, y) {
    y_expanded <- y$view(c(-1L, 1L))$expand_as(mu)
    log_pi <- torch::nnf_log_softmax(logits, dim = 2L)
    log_pdf <- -0.5 * log(2 * pi) - torch::torch_log(sigma) -
        0.5 * ((y_expanded - mu) / sigma)^2
    -torch::torch_logsumexp(log_pi + log_pdf, dim = 2L)$mean()
}


clone_torch_state <- function(state) {
    lapply(state, function(x) x$clone())
}


tensor_to_matrix <- function(x) {
    z <- as.array(x$detach()$cpu())
    d <- dim(z)
    if (is.null(d)) {
        matrix(as.numeric(z), nrow = 1L)
    } else if (length(d) == 1L) {
        matrix(as.numeric(z), nrow = 1L)
    } else {
        matrix(as.numeric(z), nrow = d[1L], ncol = d[2L])
    }
}


resolve_mdn_device <- function(device = "cpu") {
    if (length(device) != 1L || is.na(device) ||
        !device %in% c("cpu", "cuda", "auto")) {
        stop("mdn_device must be cpu, cuda, or auto.")
    }
    if (device == "cpu") return("cpu")
    available <- torch::cuda_is_available()
    if (device == "cuda" && !available) {
        stop("mdn_device=cuda requires a working CUDA runtime and an allocated GPU.")
    }
    if (available) "cuda" else "cpu"
}


fit_mdn_density <- function(X, Y,
                            n_components = 5,
                            hidden = c(64L, 64L),
                            epochs = 500,
                            batch_size = 256,
                            lr = 1e-3,
                            weight_decay = 1e-4,
                            patience = 30,
                            val_prop = 0.2,
                            min_sigma = 0.02,
                            grid_tail_prob = 5e-4,
                            grid_sd_mult = 6,
                            seed = NULL,
                            verbose = TRUE,
                            device = "cpu") {
    if (!requireNamespace("torch", quietly = TRUE)) {
        stop("Package 'torch' is required for MDN Stage 1.")
    }
    if (!torch::torch_is_installed()) {
        stop("MDN requires the torch runtime; install it with torch::install_torch().")
    }
    device_requested <- device
    device <- resolve_mdn_device(device)

    params <- validate_mdn_params(
        n_components = n_components,
        hidden = hidden,
        epochs = epochs,
        batch_size = batch_size,
        lr = lr,
        weight_decay = weight_decay,
        patience = patience,
        val_prop = val_prop,
        min_sigma = min_sigma,
        grid_tail_prob = grid_tail_prob,
        grid_sd_mult = grid_sd_mult
    )

    if (!is.null(seed)) {
        set.seed(seed)
        torch::torch_manual_seed(seed)
    }

    X <- as.matrix(X)
    Y <- as.numeric(Y)
    complete <- stats::complete.cases(X) & is.finite(Y)
    X <- X[complete, , drop = FALSE]
    Y <- Y[complete]
    n <- length(Y)
    if (n < 10L) {
        stop("MDN Stage 1 needs at least 10 complete training observations.")
    }

    x_center <- colMeans(X)
    x_scale <- apply(X, 2, stats::sd)
    x_scale[!is.finite(x_scale) | x_scale == 0] <- 1
    X_scaled <- sweep(sweep(X, 2, x_center, "-"), 2, x_scale, "/")

    y_center <- mean(Y)
    y_scale <- stats::sd(Y)
    if (!is.finite(y_scale) || y_scale <= 0) {
        y_scale <- 1
    }
    Y_scaled <- as.numeric((Y - y_center) / y_scale)
    min_sigma_scaled <- max(params$min_sigma / y_scale, 1e-4)

    n_val <- floor(n * params$val_prop)
    if (n_val >= 1L && n - n_val >= 5L) {
        val_idx <- sort(sample.int(n, n_val))
        train_idx <- setdiff(seq_len(n), val_idx)
    } else {
        val_idx <- integer(0)
        train_idx <- seq_len(n)
    }

    # Keep the CUDA data resident across epochs instead of copying every batch.
    # The CPU branch retains the original tensor construction and RNG sequence.
    if (device == "cuda") {
        x_device <- torch::torch_tensor(X_scaled, dtype = torch::torch_float(), device = device)
        y_device <- torch::torch_tensor(matrix(Y_scaled, ncol = 1L),
                                       dtype = torch::torch_float(), device = device)
    }

    MDNDensityNet <- build_mdn_torch_module()
    net <- MDNDensityNet(
        input_dim = ncol(X_scaled),
        hidden_units = params$hidden,
        n_components = params$n_components,
        min_sigma_scaled = min_sigma_scaled
    )
    if (device == "cuda") net$to(device = device)
    optimizer <- torch::optim_adam(
        net$parameters,
        lr = params$lr,
        weight_decay = params$weight_decay
    )

    if (verbose) {
        cat(sprintf(
            "  [MDN] n=%d, p=%d, K=%d, hidden=%s, device=%s (requested=%s)\n",
            n,
            ncol(X_scaled),
            params$n_components,
            paste(params$hidden, collapse = ","),
            device,
            device_requested
        ))
    }

    eval_loss <- function(idx) {
        net$eval()
        if (device == "cuda") {
            x_t <- x_device[idx, , drop = FALSE]
            y_t <- y_device[idx, , drop = FALSE]
        } else {
            x_t <- torch::torch_tensor(X_scaled[idx, , drop = FALSE],
                                       dtype = torch::torch_float())
            y_t <- torch::torch_tensor(matrix(Y_scaled[idx], ncol = 1L),
                                       dtype = torch::torch_float())
        }
        torch::with_no_grad({
            out <- net(x_t)
            as.numeric(mdn_nll(out$logits, out$mu, out$sigma, y_t)$item())
        })
    }

    best_loss <- Inf
    best_state <- NULL
    stale_epochs <- 0L
    history <- data.frame(
        epoch = integer(0),
        train_loss = numeric(0),
        val_loss = numeric(0)
    )

    for (epoch in seq_len(params$epochs)) {
        net$train()
        idx_order <- sample(train_idx, length(train_idx))
        batch_starts <- seq.int(1L, length(idx_order), by = params$batch_size)
        batch_losses <- numeric(length(batch_starts))

        for (b in seq_along(batch_starts)) {
            start <- batch_starts[b]
            end <- min(start + params$batch_size - 1L, length(idx_order))
            batch <- idx_order[start:end]

            if (device == "cuda") {
                x_b <- x_device[batch, , drop = FALSE]
                y_b <- y_device[batch, , drop = FALSE]
            } else {
                x_b <- torch::torch_tensor(X_scaled[batch, , drop = FALSE],
                                           dtype = torch::torch_float())
                y_b <- torch::torch_tensor(matrix(Y_scaled[batch], ncol = 1L),
                                           dtype = torch::torch_float())
            }

            optimizer$zero_grad()
            out <- net(x_b)
            loss <- mdn_nll(out$logits, out$mu, out$sigma, y_b)
            loss$backward()
            optimizer$step()
            batch_losses[b] <- as.numeric(loss$item())
        }

        train_loss <- mean(batch_losses)
        val_loss <- if (length(val_idx) > 0L) eval_loss(val_idx) else train_loss
        history <- rbind(
            history,
            data.frame(epoch = epoch, train_loss = train_loss, val_loss = val_loss)
        )

        if (is.finite(val_loss) && val_loss < best_loss - 1e-5) {
            best_loss <- val_loss
            best_state <- clone_torch_state(net$state_dict())
            stale_epochs <- 0L
        } else {
            stale_epochs <- stale_epochs + 1L
        }

        if (verbose && (epoch == 1L || epoch %% 25L == 0L)) {
            cat(sprintf(
                "    MDN epoch %d/%d, train=%.4f, val=%.4f\n",
                epoch,
                params$epochs,
                train_loss,
                val_loss
            ))
        }

        if (stale_epochs >= params$patience) {
            if (verbose) {
                cat(sprintf("    MDN early stop at epoch %d\n", epoch))
            }
            break
        }
    }

    if (!is.null(best_state)) {
        net$load_state_dict(best_state)
    }

    structure(
        list(
            net = net,
            x_center = x_center,
            x_scale = x_scale,
            y_center = y_center,
            y_scale = y_scale,
            n_fit = n,
            params = params,
            min_sigma_scaled = min_sigma_scaled,
            device = device,
            device_requested = device_requested,
            history = history,
            best_val_loss = best_loss,
            epochs_trained = nrow(history)
        ),
        class = "mdn_density"
    )
}


predict_mdn_params <- function(model, X_new) {
    X_new <- as.matrix(X_new)
    X_scaled <- sweep(sweep(X_new, 2, model$x_center, "-"), 2, model$x_scale, "/")
    # Read the network's actual device, including models fitted before this option existed.
    device <- model$net$parameters[[1L]]$device
    x_t <- torch::torch_tensor(X_scaled, dtype = torch::torch_float(), device = device)

    model$net$eval()
    out <- torch::with_no_grad({
        model$net(x_t)
    })

    logits <- tensor_to_matrix(out$logits)
    mu_scaled <- tensor_to_matrix(out$mu)
    sigma_scaled <- tensor_to_matrix(out$sigma)

    log_norm <- apply(logits, 1L, logsumexp)
    pi <- exp(logits - log_norm)
    mu <- model$y_center + model$y_scale * mu_scaled
    sigma <- pmax(model$y_scale * sigma_scaled, model$params$min_sigma)

    list(pi = pi, mu = mu, sigma = sigma)
}


mdn_component_range <- function(pi, mu, sigma, sd_mult) {
    active <- is.finite(pi) & is.finite(mu) & is.finite(sigma) &
        pi > 1e-8 & sigma > 0
    if (!any(active)) {
        active <- is.finite(mu) & is.finite(sigma) & sigma > 0
    }
    if (!any(active)) {
        return(c(lower = NA_real_, upper = NA_real_))
    }
    c(
        lower = min(mu[active] - sd_mult * sigma[active]),
        upper = max(mu[active] + sd_mult * sigma[active])
    )
}


mdn_param_vectors <- function(params, row_index) {
    list(
        pi = as.numeric(params$pi[row_index, ]),
        mu = as.numeric(params$mu[row_index, ]),
        sigma = as.numeric(params$sigma[row_index, ])
    )
}


mdn_one_row_params <- function(vectors) {
    list(
        pi = matrix(vectors$pi, nrow = 1L),
        mu = matrix(vectors$mu, nrow = 1L),
        sigma = matrix(vectors$sigma, nrow = 1L)
    )
}


mdn_mixture_quantile <- function(prob, pi, mu, sigma, lower, upper) {
    cdf_minus_prob <- function(y) {
        sum(pi * stats::pnorm(y, mean = mu, sd = sigma)) - prob
    }

    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
        return(NA_real_)
    }

    width <- upper - lower
    for (attempt in seq_len(8L)) {
        f_lower <- cdf_minus_prob(lower)
        f_upper <- cdf_minus_prob(upper)
        if (is.finite(f_lower) && is.finite(f_upper) &&
            f_lower <= 0 && f_upper >= 0) {
            break
        }
        lower <- lower - width
        upper <- upper + width
        width <- width * 2
    }

    tryCatch(
        stats::uniroot(cdf_minus_prob, lower = lower, upper = upper)$root,
        error = function(e) NA_real_
    )
}


mdn_y_grid_from_params <- function(model, params, row_index,
                                   n_grid = 500,
                                   grid_resolution = 0.05,
                                   max_grid_size = 2000) {
    vectors <- mdn_param_vectors(params, row_index)
    pi <- vectors$pi
    mu <- vectors$mu
    sigma <- vectors$sigma

    range <- mdn_component_range(pi, mu, sigma, model$params$grid_sd_mult)
    lower <- mdn_mixture_quantile(
        model$params$grid_tail_prob,
        pi,
        mu,
        sigma,
        range["lower"],
        range["upper"]
    )
    upper <- mdn_mixture_quantile(
        1 - model$params$grid_tail_prob,
        pi,
        mu,
        sigma,
        range["lower"],
        range["upper"]
    )

    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
        lower <- range["lower"]
        upper <- range["upper"]
    }
    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
        fallback_scale <- max(c(model$y_scale, sigma, model$params$min_sigma, 1), na.rm = TRUE)
        lower <- model$y_center - model$params$grid_sd_mult * fallback_scale
        upper <- model$y_center + model$params$grid_sd_mult * fallback_scale
    }

    n_grid_i <- min(max(n_grid, ceiling((upper - lower) / grid_resolution)), max_grid_size)
    y_grid <- seq(lower, upper, length.out = n_grid_i)
    attr(y_grid, "mdn_params") <- mdn_one_row_params(vectors)
    y_grid
}


predict_mdn_density_from_params <- function(model, params, row_index, y_grid) {
    vectors <- mdn_param_vectors(params, row_index)
    pi <- vectors$pi
    mu <- vectors$mu
    sigma <- vectors$sigma
    density_grid <- rep(0, length(y_grid))
    for (k in seq_along(pi)) {
        density_grid <- density_grid +
            pi[k] * stats::dnorm(y_grid, mean = mu[k], sd = sigma[k])
    }

    if (any(!is.finite(density_grid)) || sum(density_grid) <= 0) {
        fallback_scale <- max(c(model$y_scale, model$params$min_sigma, 1), na.rm = TRUE)
        density_grid <- stats::dnorm(y_grid, mean = model$y_center, sd = fallback_scale)
    }

    list(density_mean = density_grid, density_draws = NULL)
}


mdn_y_grid_for_x <- function(model, X_new, n_grid = 500,
                             grid_resolution = 0.05,
                             max_grid_size = 2000) {
    params <- predict_mdn_params(model, X_new[1L, , drop = FALSE])
    mdn_y_grid_from_params(model, params, 1L, n_grid, grid_resolution, max_grid_size)
}


predict_mdn_density <- function(model, X_new, y_grid) {
    params <- attr(y_grid, "mdn_params", exact = TRUE)
    if (is.null(params)) {
        params <- predict_mdn_params(model, X_new[1L, , drop = FALSE])
    }
    predict_mdn_density_from_params(model, params, 1L, y_grid)
}


normalize_nsf_hidden <- function(hidden) {
    if (is.null(hidden)) {
        return(integer(0))
    }
    if (is.character(hidden)) {
        hidden <- unlist(strsplit(hidden, ",", fixed = TRUE))
        hidden <- trimws(hidden)
        hidden <- hidden[nzchar(hidden)]
        if (length(hidden) == 0) {
            return(integer(0))
        }
    }
    hidden <- as.integer(hidden)
    if (any(is.na(hidden)) || any(hidden < 1L)) {
        stop("nsf_hidden must contain positive hidden-layer sizes.")
    }
    hidden
}


validate_nsf_params <- function(n_bins,
                                hidden,
                                epochs,
                                batch_size,
                                lr,
                                weight_decay,
                                patience,
                                val_prop,
                                tail_bound,
                                min_bin_width,
                                min_bin_height,
                                min_derivative,
                                grid_tail_prob,
                                grid_sd_mult) {
    if (length(n_bins) != 1L || is.na(n_bins) || n_bins < 2L) {
        stop("nsf_n_bins must be at least 2.")
    }
    if (length(epochs) != 1L || is.na(epochs) || epochs < 1L) {
        stop("nsf_epochs must be a positive integer.")
    }
    if (length(batch_size) != 1L || is.na(batch_size) || batch_size < 1L) {
        stop("nsf_batch_size must be a positive integer.")
    }
    if (length(lr) != 1L || is.na(lr) || lr <= 0) {
        stop("nsf_lr must be positive.")
    }
    if (length(weight_decay) != 1L || is.na(weight_decay) || weight_decay < 0) {
        stop("nsf_weight_decay must be nonnegative.")
    }
    if (length(patience) != 1L || is.na(patience) || patience < 1L) {
        stop("nsf_patience must be a positive integer.")
    }
    if (length(val_prop) != 1L || is.na(val_prop) || val_prop < 0 || val_prop >= 0.5) {
        stop("nsf_val_prop must be in [0, 0.5).")
    }
    if (length(tail_bound) != 1L || is.na(tail_bound) || tail_bound <= 0) {
        stop("nsf_tail_bound must be positive.")
    }
    if (length(min_bin_width) != 1L || is.na(min_bin_width) || min_bin_width <= 0) {
        stop("nsf_min_bin_width must be positive.")
    }
    if (length(min_bin_height) != 1L || is.na(min_bin_height) || min_bin_height <= 0) {
        stop("nsf_min_bin_height must be positive.")
    }
    if (length(min_derivative) != 1L || is.na(min_derivative) || min_derivative <= 0) {
        stop("nsf_min_derivative must be positive.")
    }
    if (as.integer(n_bins) * min_bin_width >= 2 * tail_bound) {
        stop("nsf_n_bins * nsf_min_bin_width must be smaller than 2 * nsf_tail_bound.")
    }
    if (as.integer(n_bins) * min_bin_height >= 2 * tail_bound) {
        stop("nsf_n_bins * nsf_min_bin_height must be smaller than 2 * nsf_tail_bound.")
    }
    if (length(grid_tail_prob) != 1L || is.na(grid_tail_prob) ||
        grid_tail_prob <= 0 || grid_tail_prob >= 0.1) {
        stop("nsf_grid_tail_prob must be in (0, 0.1).")
    }
    if (length(grid_sd_mult) != 1L || is.na(grid_sd_mult) || grid_sd_mult <= 0) {
        stop("nsf_grid_sd_mult must be positive.")
    }

    list(
        n_bins = as.integer(n_bins),
        hidden = normalize_nsf_hidden(hidden),
        epochs = as.integer(epochs),
        batch_size = as.integer(batch_size),
        lr = as.numeric(lr),
        weight_decay = as.numeric(weight_decay),
        patience = as.integer(patience),
        val_prop = as.numeric(val_prop),
        tail_bound = as.numeric(tail_bound),
        min_bin_width = as.numeric(min_bin_width),
        min_bin_height = as.numeric(min_bin_height),
        min_derivative = as.numeric(min_derivative),
        grid_tail_prob = as.numeric(grid_tail_prob),
        grid_sd_mult = as.numeric(grid_sd_mult)
    )
}


build_nsf_torch_module <- function() {
    torch::nn_module(
        "NSFDensityNet",
        initialize = function(input_dim, hidden_units, n_bins) {
            self$n_bins <- as.integer(n_bins)

            layers <- list()
            last_dim <- as.integer(input_dim)
            if (length(hidden_units) > 0L) {
                for (h in hidden_units) {
                    layers[[length(layers) + 1L]] <- torch::nn_linear(last_dim, h)
                    layers[[length(layers) + 1L]] <- torch::nn_relu()
                    last_dim <- h
                }
            }
            self$net <- if (length(layers) > 0L) {
                do.call(torch::nn_sequential, layers)
            } else {
                NULL
            }
            self$head <- torch::nn_linear(last_dim, 3L * self$n_bins - 1L)
        },
        forward = function(x) {
            h <- if (is.null(self$net)) x else self$net(x)
            self$head(h)
        }
    )
}


nsf_normalize_spline_params_torch <- function(raw, params) {
    k <- params$n_bins
    tail_bound <- params$tail_bound
    raw_widths <- raw[, seq_len(k), drop = FALSE]
    raw_heights <- raw[, k + seq_len(k), drop = FALSE]
    raw_derivatives <- raw[, 2L * k + seq_len(k - 1L), drop = FALSE]

    widths <- torch::nnf_softmax(raw_widths, dim = 2L) *
        (2 * tail_bound - k * params$min_bin_width) + params$min_bin_width
    heights <- torch::nnf_softmax(raw_heights, dim = 2L) *
        (2 * tail_bound - k * params$min_bin_height) + params$min_bin_height
    internal_derivatives <- torch::nnf_softplus(raw_derivatives) + params$min_derivative
    ones <- raw[, 1L, drop = FALSE] * 0 + 1
    derivatives <- torch::torch_cat(list(ones, internal_derivatives, ones), dim = 2L)

    zero <- raw[, 1L, drop = FALSE] * 0
    x_knots <- torch::torch_cat(
        list(zero - tail_bound, zero - tail_bound + torch::torch_cumsum(widths, dim = 2L)),
        dim = 2L
    )
    y_knots <- torch::torch_cat(
        list(zero - tail_bound, zero - tail_bound + torch::torch_cumsum(heights, dim = 2L)),
        dim = 2L
    )

    list(
        widths = widths,
        heights = heights,
        derivatives = derivatives,
        x_knots = x_knots,
        y_knots = y_knots
    )
}


nsf_rqs_inverse_torch <- function(y, spline, params, eps = 1e-6) {
    k <- params$n_bins
    tail_bound <- params$tail_bound
    y <- y$view(c(-1L, 1L))
    y_inside <- torch::torch_clamp(y, min = -tail_bound + eps, max = tail_bound - eps)

    y_left <- spline$y_knots[, seq_len(k), drop = FALSE]
    y_right <- spline$y_knots[, seq_len(k) + 1L, drop = FALSE]
    x_left <- spline$x_knots[, seq_len(k), drop = FALSE]
    widths <- spline$widths
    heights <- spline$heights
    d_left <- spline$derivatives[, seq_len(k), drop = FALSE]
    d_right <- spline$derivatives[, seq_len(k) + 1L, drop = FALSE]
    slopes <- heights / widths

    bin_mask <- (y_inside >= y_left) & (y_inside < y_right)
    y_frac <- torch::torch_clamp((y_inside - y_left) / heights, min = eps, max = 1 - eps)
    slope_gap <- d_right + d_left - 2 * slopes
    a_coef <- -y_frac * slope_gap - slopes + d_left
    b_coef <- y_frac * slope_gap - d_left
    c_coef <- y_frac * slopes
    disc <- torch::torch_clamp(b_coef^2 - 4 * a_coef * c_coef, min = eps)
    sqrt_disc <- torch::torch_sqrt(disc)

    safe_a <- torch::torch_where(torch::torch_abs(a_coef) < 1e-8, a_coef * 0 + 1, a_coef)
    safe_b <- torch::torch_where(torch::torch_abs(b_coef) < 1e-8, b_coef * 0 + eps, b_coef)
    root1 <- (-b_coef + sqrt_disc) / (2 * safe_a)
    root2 <- (-b_coef - sqrt_disc) / (2 * safe_a)
    linear_root <- -c_coef / safe_b
    root <- torch::torch_where(
        (root1 >= -eps) & (root1 <= 1 + eps),
        root1,
        root2
    )
    root <- torch::torch_where(torch::torch_abs(a_coef) < 1e-8, linear_root, root)
    root <- torch::torch_clamp(root, min = 0, max = 1)

    z_all <- x_left + root * widths
    denominator <- slopes + slope_gap * root * (1 - root)
    derivative <- slopes^2 *
        (d_right * root^2 + 2 * slopes * root * (1 - root) + d_left * (1 - root)^2) /
        denominator^2
    log_abs_det_all <- -torch::torch_log(torch::torch_clamp(derivative, min = eps))

    zeros <- z_all * 0
    z_inside <- torch::torch_sum(torch::torch_where(bin_mask, z_all, zeros),
                                 dim = 2L, keepdim = TRUE)
    log_abs_det_inside <- torch::torch_sum(
        torch::torch_where(bin_mask, log_abs_det_all, zeros),
        dim = 2L,
        keepdim = TRUE
    )

    inside <- (y > -tail_bound) & (y < tail_bound)
    list(
        z = torch::torch_where(inside, z_inside, y),
        log_abs_det = torch::torch_where(inside, log_abs_det_inside, y * 0)
    )
}


nsf_nll <- function(raw, y, params) {
    spline <- nsf_normalize_spline_params_torch(raw, params)
    inv <- nsf_rqs_inverse_torch(y, spline, params)
    base_log <- -0.5 * log(2 * pi) - 0.5 * inv$z^2
    log_prob <- base_log + inv$log_abs_det
    -log_prob$mean()
}


fit_nsf_density <- function(X, Y,
                            n_bins = 12,
                            hidden = c(128L, 128L),
                            epochs = 500,
                            batch_size = 256,
                            lr = 1e-3,
                            weight_decay = 1e-4,
                            patience = 30,
                            val_prop = 0.2,
                            tail_bound = 5,
                            min_bin_width = 1e-3,
                            min_bin_height = 1e-3,
                            min_derivative = 1e-3,
                            grid_tail_prob = 5e-4,
                            grid_sd_mult = 6,
                            seed = NULL,
                            verbose = TRUE) {
    if (!requireNamespace("torch", quietly = TRUE)) {
        stop("Package 'torch' is required for NSF Stage 1.")
    }

    params <- validate_nsf_params(
        n_bins = n_bins,
        hidden = hidden,
        epochs = epochs,
        batch_size = batch_size,
        lr = lr,
        weight_decay = weight_decay,
        patience = patience,
        val_prop = val_prop,
        tail_bound = tail_bound,
        min_bin_width = min_bin_width,
        min_bin_height = min_bin_height,
        min_derivative = min_derivative,
        grid_tail_prob = grid_tail_prob,
        grid_sd_mult = grid_sd_mult
    )

    if (!is.null(seed)) {
        set.seed(seed)
        torch::torch_manual_seed(seed)
    }

    X <- as.matrix(X)
    Y <- as.numeric(Y)
    complete <- stats::complete.cases(X) & is.finite(Y)
    X <- X[complete, , drop = FALSE]
    Y <- Y[complete]
    n <- length(Y)
    if (n < 10L) {
        stop("NSF Stage 1 needs at least 10 complete training observations.")
    }

    x_center <- colMeans(X)
    x_scale <- apply(X, 2, stats::sd)
    x_scale[!is.finite(x_scale) | x_scale == 0] <- 1
    X_scaled <- sweep(sweep(X, 2, x_center, "-"), 2, x_scale, "/")

    y_center <- mean(Y)
    y_scale <- stats::sd(Y)
    if (!is.finite(y_scale) || y_scale <= 0) {
        y_scale <- 1
    }
    Y_scaled <- as.numeric((Y - y_center) / y_scale)
    observed_tail <- max(abs(Y_scaled), na.rm = TRUE) + 0.25
    if (is.finite(observed_tail)) {
        params$tail_bound <- max(params$tail_bound, observed_tail)
    }
    validate_nsf_params(
        n_bins = params$n_bins,
        hidden = params$hidden,
        epochs = params$epochs,
        batch_size = params$batch_size,
        lr = params$lr,
        weight_decay = params$weight_decay,
        patience = params$patience,
        val_prop = params$val_prop,
        tail_bound = params$tail_bound,
        min_bin_width = params$min_bin_width,
        min_bin_height = params$min_bin_height,
        min_derivative = params$min_derivative,
        grid_tail_prob = params$grid_tail_prob,
        grid_sd_mult = params$grid_sd_mult
    )

    n_val <- floor(n * params$val_prop)
    if (n_val >= 1L && n - n_val >= 5L) {
        val_idx <- sort(sample.int(n, n_val))
        train_idx <- setdiff(seq_len(n), val_idx)
    } else {
        val_idx <- integer(0)
        train_idx <- seq_len(n)
    }

    NSFDensityNet <- build_nsf_torch_module()
    net <- NSFDensityNet(
        input_dim = ncol(X_scaled),
        hidden_units = params$hidden,
        n_bins = params$n_bins
    )
    optimizer <- torch::optim_adam(
        net$parameters,
        lr = params$lr,
        weight_decay = params$weight_decay
    )

    if (verbose) {
        cat(sprintf(
            "  [NSF] n=%d, p=%d, bins=%d, hidden=%s, tail_bound=%.3f\n",
            n,
            ncol(X_scaled),
            params$n_bins,
            paste(params$hidden, collapse = ","),
            params$tail_bound
        ))
    }

    eval_loss <- function(idx) {
        net$eval()
        x_t <- torch::torch_tensor(X_scaled[idx, , drop = FALSE],
                                   dtype = torch::torch_float())
        y_t <- torch::torch_tensor(matrix(Y_scaled[idx], ncol = 1L),
                                   dtype = torch::torch_float())
        raw <- net(x_t)
        as.numeric(nsf_nll(raw, y_t, params)$item())
    }

    best_loss <- Inf
    best_state <- NULL
    stale_epochs <- 0L
    history <- data.frame(
        epoch = integer(0),
        train_loss = numeric(0),
        val_loss = numeric(0)
    )

    for (epoch in seq_len(params$epochs)) {
        net$train()
        idx_order <- sample(train_idx, length(train_idx))
        batch_starts <- seq.int(1L, length(idx_order), by = params$batch_size)
        batch_losses <- numeric(length(batch_starts))

        for (b in seq_along(batch_starts)) {
            start <- batch_starts[b]
            end <- min(start + params$batch_size - 1L, length(idx_order))
            batch <- idx_order[start:end]

            x_b <- torch::torch_tensor(X_scaled[batch, , drop = FALSE],
                                       dtype = torch::torch_float())
            y_b <- torch::torch_tensor(matrix(Y_scaled[batch], ncol = 1L),
                                       dtype = torch::torch_float())

            optimizer$zero_grad()
            raw <- net(x_b)
            loss <- nsf_nll(raw, y_b, params)
            loss$backward()
            optimizer$step()
            batch_losses[b] <- as.numeric(loss$item())
        }

        train_loss <- mean(batch_losses)
        val_loss <- if (length(val_idx) > 0L) eval_loss(val_idx) else train_loss
        history <- rbind(
            history,
            data.frame(epoch = epoch, train_loss = train_loss, val_loss = val_loss)
        )

        if (is.finite(val_loss) && val_loss < best_loss - 1e-5) {
            best_loss <- val_loss
            best_state <- clone_torch_state(net$state_dict())
            stale_epochs <- 0L
        } else {
            stale_epochs <- stale_epochs + 1L
        }

        if (verbose && (epoch == 1L || epoch %% 25L == 0L)) {
            cat(sprintf(
                "    NSF epoch %d/%d, train=%.4f, val=%.4f\n",
                epoch,
                params$epochs,
                train_loss,
                val_loss
            ))
        }

        if (stale_epochs >= params$patience) {
            if (verbose) {
                cat(sprintf("    NSF early stop at epoch %d\n", epoch))
            }
            break
        }
    }

    if (!is.null(best_state)) {
        net$load_state_dict(best_state)
    }

    structure(
        list(
            net = net,
            x_center = x_center,
            x_scale = x_scale,
            y_center = y_center,
            y_scale = y_scale,
            n_fit = n,
            params = params,
            history = history,
            best_val_loss = best_loss,
            epochs_trained = nrow(history)
        ),
        class = "nsf_density"
    )
}


row_softmax <- function(x) {
    x <- as.matrix(x)
    row_max <- apply(x, 1L, max)
    ex <- exp(x - row_max)
    ex / rowSums(ex)
}


softplus_numeric <- function(x) {
    log1p(exp(-abs(x))) + pmax(x, 0)
}


normalize_nsf_params_numeric <- function(raw, params) {
    raw <- as.matrix(raw)
    k <- params$n_bins
    raw_widths <- raw[, seq_len(k), drop = FALSE]
    raw_heights <- raw[, k + seq_len(k), drop = FALSE]
    raw_derivatives <- raw[, 2L * k + seq_len(k - 1L), drop = FALSE]

    widths <- row_softmax(raw_widths) *
        (2 * params$tail_bound - k * params$min_bin_width) + params$min_bin_width
    heights <- row_softmax(raw_heights) *
        (2 * params$tail_bound - k * params$min_bin_height) + params$min_bin_height
    internal_derivatives <- softplus_numeric(raw_derivatives) + params$min_derivative
    derivatives <- cbind(1, internal_derivatives, 1)

    x_knots <- cbind(-params$tail_bound,
                     -params$tail_bound + t(apply(widths, 1L, cumsum)))
    y_knots <- cbind(-params$tail_bound,
                     -params$tail_bound + t(apply(heights, 1L, cumsum)))

    list(
        widths = widths,
        heights = heights,
        derivatives = derivatives,
        x_knots = x_knots,
        y_knots = y_knots
    )
}


predict_nsf_params <- function(model, X_new) {
    X_new <- as.matrix(X_new)
    X_scaled <- sweep(sweep(X_new, 2, model$x_center, "-"), 2, model$x_scale, "/")
    x_t <- torch::torch_tensor(X_scaled, dtype = torch::torch_float())

    model$net$eval()
    raw <- torch::with_no_grad({
        model$net(x_t)
    })

    normalize_nsf_params_numeric(tensor_to_matrix(raw), model$params)
}


nsf_param_vectors <- function(params, row_index) {
    list(
        widths = as.numeric(params$widths[row_index, ]),
        heights = as.numeric(params$heights[row_index, ]),
        derivatives = as.numeric(params$derivatives[row_index, ]),
        x_knots = as.numeric(params$x_knots[row_index, ]),
        y_knots = as.numeric(params$y_knots[row_index, ])
    )
}


nsf_one_row_params <- function(vectors) {
    list(
        widths = matrix(vectors$widths, nrow = 1L),
        heights = matrix(vectors$heights, nrow = 1L),
        derivatives = matrix(vectors$derivatives, nrow = 1L),
        x_knots = matrix(vectors$x_knots, nrow = 1L),
        y_knots = matrix(vectors$y_knots, nrow = 1L)
    )
}


nsf_rqs_forward_numeric <- function(z, vectors, tail_bound, eps = 1e-8) {
    z <- as.numeric(z)
    out <- z
    inside <- z > -tail_bound & z < tail_bound
    if (!any(inside)) {
        return(out)
    }

    z_inside <- pmin(pmax(z[inside], -tail_bound + eps), tail_bound - eps)
    bin <- findInterval(z_inside, vectors$x_knots, rightmost.closed = FALSE)
    bin <- pmin(pmax(bin, 1L), length(vectors$widths))

    x_left <- vectors$x_knots[bin]
    y_left <- vectors$y_knots[bin]
    widths <- vectors$widths[bin]
    heights <- vectors$heights[bin]
    d_left <- vectors$derivatives[bin]
    d_right <- vectors$derivatives[bin + 1L]
    slopes <- heights / widths
    theta <- pmin(pmax((z_inside - x_left) / widths, 0), 1)
    theta_one_minus <- theta * (1 - theta)
    numerator <- slopes * theta^2 + d_left * theta_one_minus
    denominator <- slopes + (d_right + d_left - 2 * slopes) * theta_one_minus
    out[inside] <- y_left + heights * numerator / denominator
    out
}


nsf_rqs_inverse_numeric <- function(y, vectors, tail_bound, eps = 1e-8) {
    y <- as.numeric(y)
    z <- y
    log_abs_det <- rep(0, length(y))
    inside <- y > -tail_bound & y < tail_bound
    if (!any(inside)) {
        return(list(z = z, log_abs_det = log_abs_det))
    }

    y_inside <- pmin(pmax(y[inside], -tail_bound + eps), tail_bound - eps)
    bin <- findInterval(y_inside, vectors$y_knots, rightmost.closed = FALSE)
    bin <- pmin(pmax(bin, 1L), length(vectors$widths))

    y_left <- vectors$y_knots[bin]
    x_left <- vectors$x_knots[bin]
    widths <- vectors$widths[bin]
    heights <- vectors$heights[bin]
    d_left <- vectors$derivatives[bin]
    d_right <- vectors$derivatives[bin + 1L]
    slopes <- heights / widths
    y_frac <- pmin(pmax((y_inside - y_left) / heights, eps), 1 - eps)
    slope_gap <- d_right + d_left - 2 * slopes
    a_coef <- -y_frac * slope_gap - slopes + d_left
    b_coef <- y_frac * slope_gap - d_left
    c_coef <- y_frac * slopes
    disc <- pmax(b_coef^2 - 4 * a_coef * c_coef, eps)
    sqrt_disc <- sqrt(disc)
    safe_a <- ifelse(abs(a_coef) < 1e-8, 1, a_coef)
    safe_b <- ifelse(abs(b_coef) < 1e-8, ifelse(b_coef >= 0, eps, -eps), b_coef)
    root1 <- (-b_coef + sqrt_disc) / (2 * safe_a)
    root2 <- (-b_coef - sqrt_disc) / (2 * safe_a)
    root <- ifelse(root1 >= -eps & root1 <= 1 + eps, root1, root2)
    linear <- abs(a_coef) < 1e-8
    linear_root <- -c_coef / safe_b
    root[linear] <- linear_root[linear]
    root <- pmin(pmax(root, 0), 1)

    z_inside <- x_left + root * widths
    denominator <- slopes + slope_gap * root * (1 - root)
    derivative <- slopes^2 *
        (d_right * root^2 + 2 * slopes * root * (1 - root) + d_left * (1 - root)^2) /
        denominator^2
    z[inside] <- z_inside
    log_abs_det[inside] <- -log(pmax(derivative, eps))
    list(z = z, log_abs_det = log_abs_det)
}


nsf_y_grid_from_params <- function(model, params, row_index,
                                   n_grid = 500,
                                   grid_resolution = 0.05,
                                   max_grid_size = 2000) {
    vectors <- nsf_param_vectors(params, row_index)
    tail_z <- stats::qnorm(c(model$params$grid_tail_prob,
                             1 - model$params$grid_tail_prob))
    y_scaled <- nsf_rqs_forward_numeric(
        tail_z,
        vectors,
        tail_bound = model$params$tail_bound
    )
    lower <- model$y_center + model$y_scale * min(y_scaled, na.rm = TRUE)
    upper <- model$y_center + model$y_scale * max(y_scaled, na.rm = TRUE)

    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
        fallback_scale <- max(c(model$y_scale, 1), na.rm = TRUE)
        lower <- model$y_center - model$params$grid_sd_mult * fallback_scale
        upper <- model$y_center + model$params$grid_sd_mult * fallback_scale
    }

    n_grid_i <- min(max(n_grid, ceiling((upper - lower) / grid_resolution)), max_grid_size)
    y_grid <- seq(lower, upper, length.out = n_grid_i)
    attr(y_grid, "nsf_params") <- nsf_one_row_params(vectors)
    y_grid
}


predict_nsf_density_from_params <- function(model, params, row_index, y_grid) {
    vectors <- nsf_param_vectors(params, row_index)
    y_scaled <- as.numeric((y_grid - model$y_center) / model$y_scale)
    inv <- nsf_rqs_inverse_numeric(
        y_scaled,
        vectors,
        tail_bound = model$params$tail_bound
    )
    density_grid <- stats::dnorm(inv$z) * exp(inv$log_abs_det) / model$y_scale

    if (any(!is.finite(density_grid)) || sum(density_grid) <= 0) {
        fallback_scale <- max(c(model$y_scale, 1), na.rm = TRUE)
        density_grid <- stats::dnorm(y_grid, mean = model$y_center, sd = fallback_scale)
    }

    list(density_mean = density_grid, density_draws = NULL)
}


nsf_y_grid_for_x <- function(model, X_new, n_grid = 500,
                             grid_resolution = 0.05,
                             max_grid_size = 2000) {
    params <- predict_nsf_params(model, X_new[1L, , drop = FALSE])
    nsf_y_grid_from_params(model, params, 1L, n_grid, grid_resolution, max_grid_size)
}


predict_nsf_density <- function(model, X_new, y_grid) {
    params <- attr(y_grid, "nsf_params", exact = TRUE)
    if (is.null(params)) {
        params <- predict_nsf_params(model, X_new[1L, , drop = FALSE])
    }
    predict_nsf_density_from_params(model, params, 1L, y_grid)
}


outcome_y_grid_for_x <- function(model, X_new, n_grid = 500,
                                 grid_resolution = 0.05,
                                 max_grid_size = 2000) {
    if (inherits(model, "joint_dp_mixture_density")) {
        return(dp_y_grid_for_x(model, X_new, n_grid, grid_resolution, max_grid_size))
    }
    if (inherits(model, "bart_predictive_density")) {
        return(bart_y_grid_for_x(model, X_new, n_grid, grid_resolution, max_grid_size))
    }
    if (inherits(model, "mdn_density")) {
        return(mdn_y_grid_for_x(model, X_new, n_grid, grid_resolution, max_grid_size))
    }
    if (inherits(model, "nsf_density")) {
        return(nsf_y_grid_for_x(model, X_new, n_grid, grid_resolution, max_grid_size))
    }
    if (inherits(model, "drbart_density")) {
        return(drbart_y_grid_for_x(model, X_new, n_grid, grid_resolution, max_grid_size))
    }
    stop(sprintf("Unsupported outcome density model class: %s", paste(class(model), collapse = ", ")))
}


predict_outcome_density <- function(model, X_new, y_grid) {
    if (inherits(model, "joint_dp_mixture_density")) {
        return(predict_joint_dp_mixture_density(model, X_new, y_grid))
    }
    if (inherits(model, "bart_predictive_density")) {
        return(predict_bart_density(model, X_new, y_grid))
    }
    if (inherits(model, "mdn_density")) {
        return(predict_mdn_density(model, X_new, y_grid))
    }
    if (inherits(model, "nsf_density")) {
        return(predict_nsf_density(model, X_new, y_grid))
    }
    if (inherits(model, "drbart_density")) {
        return(predict_drbart_density(model, X_new, y_grid))
    }
    stop(sprintf("Unsupported outcome density model class: %s", paste(class(model), collapse = ", ")))
}


fit_outcome_aps <- function(X, Y,
                            weights = NULL,
                            X_cal = NULL,
                            Y_cal = NULL,
                            weights_cal = NULL,
                            train_prop = 0.75,
                            stage1_model = "dp",
                            n_grid = 500,
                            grid_resolution = 0.05,
                            max_grid_size = 2000,
                            dp_n_fit_max = 10000,
                            dp_n_iter = 250,
                            dp_burn = 125,
                            dp_thin = 5,
                            dp_alpha = 5,
                            dp_kappa0 = 0.05,
                            dp_nu0_extra = 3,
                            dp_prior_scale = 1,
                            dp_jitter = 1e-7,
                            mdn_n_components = 5,
                            mdn_hidden = c(64L, 64L),
                            mdn_epochs = 500,
                            mdn_batch_size = 256,
                            mdn_lr = 1e-3,
                            mdn_weight_decay = 1e-4,
                            mdn_patience = 30,
                            mdn_val_prop = 0.2,
                            mdn_min_sigma = 0.02,
                            mdn_grid_tail_prob = 5e-4,
                            mdn_grid_sd_mult = 6,
                            nsf_n_bins = 12,
                            nsf_hidden = c(128L, 128L),
                            nsf_epochs = 500,
                            nsf_batch_size = 256,
                            nsf_lr = 1e-3,
                            nsf_weight_decay = 1e-4,
                            nsf_patience = 30,
                            nsf_val_prop = 0.2,
                            nsf_tail_bound = 5,
                            nsf_min_bin_width = 1e-3,
                            nsf_min_bin_height = 1e-3,
                            nsf_min_derivative = 1e-3,
                            nsf_grid_tail_prob = 5e-4,
                            nsf_grid_sd_mult = 6,
                            seed = NULL,
                            verbose = TRUE,
                            drbart_nburn = 5000,
                            drbart_nsim = 5000,
                            drbart_nthin = 1,
                            drbart_m_mean = 200,
                            drbart_m_var = 100,
                            drbart_variance = "ux",
                            drbart_tree_alpha = 0.95,
                            drbart_tree_beta = 2,
                            drbart_lambda = 1,
                            drbart_nu = 2,
                            drbart_kfac = 2,
                            drbart_phi0 = 1,
                            drbart_grid_tail_prob = 5e-4,
                            drbart_predict_batch_size = 100,
                            mdn_device = "cpu") {
    if (!stage1_model %in% c("dp", "bart", "mdn", "nsf", "drbart")) {
        stop("stage1_model must be dp, bart, mdn, nsf, or drbart.")
    }
    grid_params <- normalize_outcome_grid_params(n_grid, grid_resolution, max_grid_size)
    n <- length(Y)
    has_external_cal <- !is.null(X_cal) || !is.null(Y_cal) || !is.null(weights_cal)
    if (has_external_cal) {
        if (is.null(X_cal) || is.null(Y_cal)) {
            stop("X_cal and Y_cal must both be supplied when using external calibration data.")
        }
        X_train <- X
        Y_train <- Y
        X_val <- X_cal
        Y_val <- Y_cal
        train_idx <- seq_len(n)
        val_idx <- seq_along(Y_cal)
        weights_val <- if (is.null(weights_cal)) rep(1, length(Y_cal)) else weights_cal
        if (length(weights_val) != length(Y_cal)) {
            stop("weights_cal must have the same length as Y_cal.")
        }
    } else {
        n_train <- floor(n * train_prop)
        if (!is.null(seed)) {
            set.seed(seed)
        }
        train_idx <- sample(n, n_train)
        val_idx <- setdiff(seq_len(n), train_idx)
        X_train <- X[train_idx, , drop = FALSE]
        Y_train <- Y[train_idx]
        X_val <- X[val_idx, , drop = FALSE]
        Y_val <- Y[val_idx]
        weights_val <- if (is.null(weights)) rep(1, length(val_idx)) else weights[val_idx]
    }

    if (stage1_model == "dp") {
        density_model <- fit_joint_dp_mixture_density(
            X_train, Y_train,
            n_fit_max = dp_n_fit_max, n_iter = dp_n_iter, burn = dp_burn, thin = dp_thin,
            alpha = dp_alpha, kappa0 = dp_kappa0, nu0_extra = dp_nu0_extra,
            prior_scale = dp_prior_scale, jitter = dp_jitter, seed = seed, verbose = verbose
        )
    } else if (stage1_model == "bart") {
        density_model <- fit_bart_predictive_density(
            X_train, Y_train,
            seed = seed, verbose = verbose
        )
    } else if (stage1_model == "mdn") {
        density_model <- fit_mdn_density(
            X_train, Y_train,
            n_components = mdn_n_components,
            hidden = mdn_hidden,
            epochs = mdn_epochs,
            batch_size = mdn_batch_size,
            lr = mdn_lr,
            weight_decay = mdn_weight_decay,
            patience = mdn_patience,
            val_prop = mdn_val_prop,
            min_sigma = mdn_min_sigma,
            grid_tail_prob = mdn_grid_tail_prob,
            grid_sd_mult = mdn_grid_sd_mult,
            device = mdn_device,
            seed = seed,
            verbose = verbose
        )
    } else if (stage1_model == "nsf") {
        density_model <- fit_nsf_density(
            X_train, Y_train,
            n_bins = nsf_n_bins,
            hidden = nsf_hidden,
            epochs = nsf_epochs,
            batch_size = nsf_batch_size,
            lr = nsf_lr,
            weight_decay = nsf_weight_decay,
            patience = nsf_patience,
            val_prop = nsf_val_prop,
            tail_bound = nsf_tail_bound,
            min_bin_width = nsf_min_bin_width,
            min_bin_height = nsf_min_bin_height,
            min_derivative = nsf_min_derivative,
            grid_tail_prob = nsf_grid_tail_prob,
            grid_sd_mult = nsf_grid_sd_mult,
            seed = seed,
            verbose = verbose
        )
    } else if (stage1_model == "drbart") {
        density_model <- fit_drbart_density(
            X_train, Y_train,
            nburn = drbart_nburn,
            nsim = drbart_nsim,
            nthin = drbart_nthin,
            m_mean = drbart_m_mean,
            m_var = drbart_m_var,
            variance = drbart_variance,
            tree_alpha = drbart_tree_alpha,
            tree_beta = drbart_tree_beta,
            lambda = drbart_lambda,
            nu = drbart_nu,
            kfac = drbart_kfac,
            phi0 = drbart_phi0,
            grid_tail_prob = drbart_grid_tail_prob,
            predict_batch_size = drbart_predict_batch_size,
            seed = seed, verbose = verbose
        )
    }

    mdn_params_val <- NULL
    nsf_params_val <- NULL
    if (inherits(density_model, "mdn_density")) {
        mdn_params_val <- predict_mdn_params(density_model, X_val)
    } else if (inherits(density_model, "nsf_density")) {
        nsf_params_val <- predict_nsf_params(density_model, X_val)
    }

    scores_val <- numeric(length(val_idx))
    for (i in seq_along(val_idx)) {
        if (!is.null(mdn_params_val)) {
            y_grid <- mdn_y_grid_from_params(
                density_model,
                mdn_params_val,
                i,
                n_grid = grid_params$n_grid,
                grid_resolution = grid_params$grid_resolution,
                max_grid_size = grid_params$max_grid_size
            )
            density <- predict_mdn_density_from_params(
                density_model,
                mdn_params_val,
                i,
                y_grid
            )$density_mean
        } else if (!is.null(nsf_params_val)) {
            y_grid <- nsf_y_grid_from_params(
                density_model,
                nsf_params_val,
                i,
                n_grid = grid_params$n_grid,
                grid_resolution = grid_params$grid_resolution,
                max_grid_size = grid_params$max_grid_size
            )
            density <- predict_nsf_density_from_params(
                density_model,
                nsf_params_val,
                i,
                y_grid
            )$density_mean
        } else {
            y_grid <- outcome_y_grid_for_x(
                density_model,
                X_val[i, , drop = FALSE],
                n_grid = grid_params$n_grid,
                grid_resolution = grid_params$grid_resolution,
                max_grid_size = grid_params$max_grid_size
            )
            density <- predict_outcome_density(
                density_model,
                X_val[i, , drop = FALSE],
                y_grid
            )$density_mean
        }
        scores_val[i] <- compute_aps_score(density, y_grid, Y_val[i])
    }

    if (stage1_model == "dp") {
        dp_diagnostics <- list(
            k_trace = density_model$k_trace,
            k_mean_last_half = mean(tail(density_model$k_trace, floor(length(density_model$k_trace) / 2))),
            k_min_last_half = min(tail(density_model$k_trace, floor(length(density_model$k_trace) / 2))),
            k_max_last_half = max(tail(density_model$k_trace, floor(length(density_model$k_trace) / 2))),
            posterior_draws = length(density_model$draws),
            n_fit_used = density_model$n_fit,
            n_fit_total = density_model$n_total
        )
        bart_diagnostics <- NULL
        mdn_diagnostics <- NULL
        nsf_diagnostics <- NULL
    } else if (stage1_model == "bart") {
        dp_diagnostics <- NULL
        bart_diagnostics <- list(
            n_fit_used = density_model$n_fit
        )
        mdn_diagnostics <- NULL
        nsf_diagnostics <- NULL
    } else if (stage1_model == "mdn") {
        dp_diagnostics <- NULL
        bart_diagnostics <- NULL
        mdn_diagnostics <- list(
            n_fit_used = density_model$n_fit,
            n_components = density_model$params$n_components,
            hidden = paste(density_model$params$hidden, collapse = ","),
            epochs_trained = density_model$epochs_trained,
            best_val_loss = density_model$best_val_loss,
            final_train_loss = if (nrow(density_model$history) > 0L) {
                tail(density_model$history$train_loss, 1L)
            } else NA_real_,
            final_val_loss = if (nrow(density_model$history) > 0L) {
                tail(density_model$history$val_loss, 1L)
            } else NA_real_
        )
        nsf_diagnostics <- NULL
    } else if (stage1_model == "nsf") {
        dp_diagnostics <- NULL
        bart_diagnostics <- NULL
        mdn_diagnostics <- NULL
        nsf_diagnostics <- list(
            n_fit_used = density_model$n_fit,
            n_bins = density_model$params$n_bins,
            hidden = paste(density_model$params$hidden, collapse = ","),
            tail_bound = density_model$params$tail_bound,
            epochs_trained = density_model$epochs_trained,
            best_val_loss = density_model$best_val_loss,
            final_train_loss = if (nrow(density_model$history) > 0L) {
                tail(density_model$history$train_loss, 1L)
            } else NA_real_,
            final_val_loss = if (nrow(density_model$history) > 0L) {
                tail(density_model$history$val_loss, 1L)
            } else NA_real_
        )
    } else if (stage1_model == "drbart") {
        dp_diagnostics <- NULL
        bart_diagnostics <- NULL
        mdn_diagnostics <- NULL
        nsf_diagnostics <- NULL
    }

    structure(
        list(
            density_model = density_model,
            stage1_model = stage1_model,
            scores_val = scores_val,
            weights_val = weights_val,
            train_idx = train_idx,
            val_idx = val_idx,
            split_source = if (has_external_cal) "external" else "internal",
            score_range = range(scores_val),
            score_mean = mean(scores_val),
            grid_params = grid_params,
            dp_diagnostics = dp_diagnostics,
            bart_diagnostics = bart_diagnostics,
            mdn_diagnostics = mdn_diagnostics,
            nsf_diagnostics = nsf_diagnostics
        ),
        class = "outcome_aps_fit"
    )
}


fit_outcome_aps_oracle <- function(X, Y, scenario, arm,
                                   weights = NULL,
                                   X_cal = NULL,
                                   Y_cal = NULL,
                                   weights_cal = NULL,
                                   train_prop = 0.75,
                                   n_grid = 500,
                                   grid_resolution = 0.05,
                                   max_grid_size = 2000,
                                   seed = NULL,
                                   ...) {
    grid_params <- normalize_outcome_grid_params(n_grid, grid_resolution, max_grid_size)
    n <- length(Y)
    has_external_cal <- !is.null(X_cal) || !is.null(Y_cal) || !is.null(weights_cal)
    if (has_external_cal) {
        if (is.null(X_cal) || is.null(Y_cal)) {
            stop("X_cal and Y_cal must both be supplied when using external calibration data.")
        }
        X_val <- X_cal
        Y_val <- Y_cal
        train_idx <- seq_len(n)
        val_idx <- seq_along(Y_cal)
        weights_val <- if (is.null(weights_cal)) rep(1, length(Y_cal)) else weights_cal
        if (length(weights_val) != length(Y_cal)) {
            stop("weights_cal must have the same length as Y_cal.")
        }
    } else {
        n_train <- floor(n * train_prop)
        if (!is.null(seed)) {
            set.seed(seed)
        }
        train_idx <- sample(n, n_train)
        val_idx <- setdiff(seq_len(n), train_idx)
        X_val <- X[val_idx, , drop = FALSE]
        Y_val <- Y[val_idx]
        weights_val <- if (is.null(weights)) rep(1, length(val_idx)) else weights[val_idx]
    }

    scores_val <- numeric(length(val_idx))
    for (i in seq_along(val_idx)) {
        y_grid <- get_oracle_y_grid(
            scenario,
            X_val[i, , drop = FALSE],
            arm,
            grid_params$n_grid,
            grid_params$grid_resolution,
            grid_params$max_grid_size,
            ...
        )
        density <- get_true_outcome_density(scenario, X_val[i, , drop = FALSE], y_grid, arm, ...)
        scores_val[i] <- compute_aps_score(density, y_grid, Y_val[i])
    }

    structure(
        list(
            scenario = scenario,
            arm = arm,
            scores_val = scores_val,
            weights_val = weights_val,
            train_idx = train_idx,
            val_idx = val_idx,
            split_source = if (has_external_cal) "external" else "internal",
            score_range = range(scores_val),
            score_mean = mean(scores_val),
            grid_params = grid_params
        ),
        class = "outcome_aps_oracle_fit"
    )
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
        cutoff_i <- weighted_quantile(fit_obj$scores_val, weights_val_final, qt, use_inf = TRUE)
        if (cutoff_i == Inf) {
            intervals_list[[i]] <- normalize_interval_union(matrix(c(-Inf, Inf), nrow = 1))
            next
        }
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


predict_outcome_aps_oracle <- function(fit_obj,
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

    intervals_list <- vector("list", n_test)
    for (i in seq_len(n_test)) {
        qt <- (1 + weights_test_trunc[i] / total_weight_val) * (1 - alpha)
        cutoff_i <- weighted_quantile(fit_obj$scores_val, weights_val_final, qt, use_inf = TRUE)
        if (cutoff_i == Inf) {
            intervals_list[[i]] <- normalize_interval_union(matrix(c(-Inf, Inf), nrow = 1))
            next
        }
        y_grid <- get_oracle_y_grid(fit_obj$scenario, X_test[i, , drop = FALSE], fit_obj$arm,
                                    grid_params$n_grid,
                                    grid_params$grid_resolution,
                                    grid_params$max_grid_size,
                                    ...)
        density_grid <- get_true_outcome_density(fit_obj$scenario, X_test[i, , drop = FALSE],
                                                 y_grid, fit_obj$arm, ...)
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

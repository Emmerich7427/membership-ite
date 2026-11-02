################################################################################
# Module: 02_models.R
# Purpose: Base nuisance-model interfaces for propensity and quantile models
################################################################################

fit_propensity_score <- function(X, T, n_trees = 100) {
    if (!requireNamespace("gbm", quietly = TRUE)) {
        stop("Package 'gbm' is required for estimated propensity scores.")
    }

    X <- as.data.frame(X)
    colnames(X) <- paste0("X", seq_len(ncol(X)))
    data <- data.frame(T = T, X)

    ps_model <- gbm::gbm(
        T ~ .,
        data = data,
        distribution = "bernoulli",
        n.trees = n_trees,
        interaction.depth = 3,
        verbose = FALSE
    )

    structure(
        list(gbm_model = ps_model, n_trees = n_trees, X_names = colnames(X)),
        class = "propensity_score_model"
    )
}


predict_propensity_score <- function(model, X_new) {
    X_new <- as.data.frame(X_new)
    colnames(X_new) <- model$X_names
    ps <- gbm::predict.gbm(
        model$gbm_model,
        newdata = X_new,
        n.trees = model$n_trees,
        type = "response"
    )
    as.numeric(pmin(pmax(ps, 0.01), 0.99))
}


fit_quantile_rf <- function(X_train, Y_train, quantile, ...) {
    if (!requireNamespace("grf", quietly = TRUE)) {
        stop("Package 'grf' is required for quantile random forests.")
    }
    X_train <- as.matrix(X_train)
    qrf_fit <- grf::quantile_forest(X = X_train, Y = Y_train, quantiles = quantile, ...)
    structure(list(qrf_fit = qrf_fit, quantile = quantile), class = "quantile_rf_model")
}


predict_quantile_rf <- function(model, X_new) {
    pred <- predict(model$qrf_fit, newdata = as.matrix(X_new), quantiles = model$quantile)
    if (is.list(pred) && !is.data.frame(pred)) {
        pred <- pred$predictions
    }
    as.numeric(pred)
}


fit_quantile_bart <- function(X_train, Y_train, quantile, n_post = 200, ...) {
    if (!requireNamespace("dbarts", quietly = TRUE)) {
        stop("Package 'dbarts' is required for quantile BART.")
    }
    if (quantile <= 0 || quantile >= 1) {
        stop("quantile must be in (0, 1).")
    }

    X_train <- as.data.frame(X_train)
    colnames(X_train) <- paste0("X", seq_len(ncol(X_train)))
    bart_fit <- dbarts::bart(
        x.train = X_train,
        y.train = Y_train,
        keeptrees = TRUE,
        verbose = FALSE,
        ...
    )

    structure(
        list(
            bart_fit = bart_fit,
            quantile = quantile,
            n_post = n_post,
            X_names = colnames(X_train)
        ),
        class = "quantile_bart_model"
    )
}


predict_quantile_bart <- function(model, X_new) {
    X_new <- as.data.frame(X_new)
    colnames(X_new) <- model$X_names
    post_samples <- predict(model$bart_fit, newdata = X_new, type = "ppd")

    if (nrow(post_samples) > model$n_post) {
        idx <- round(seq(1, nrow(post_samples), length.out = model$n_post))
        post_samples <- post_samples[idx, , drop = FALSE]
    }

    as.numeric(apply(
        post_samples,
        2,
        function(x) stats::quantile(x, probs = model$quantile, names = FALSE)
    ))
}

################################################################################
# Module: 02_dgp_S3_propensity_misspec.R
# Purpose: Misspecified propensity model used only for S3
################################################################################

fit_misspecified_ps_S3 <- function(X, T) {
    X <- as.matrix(X)
    data <- data.frame(T = T, X1 = X[, 1], X2 = X[, 2])
    model <- stats::glm(T ~ X1 + X2, data = data, family = stats::binomial())
    structure(list(glm_model = model), class = "misspecified_ps_S3")
}


predict_misspecified_ps_S3 <- function(model, X_new) {
    X_new <- as.matrix(X_new)
    data_new <- data.frame(X1 = X_new[, 1], X2 = X_new[, 2])
    ps <- stats::predict(model$glm_model, newdata = data_new, type = "response")
    as.numeric(pmin(pmax(ps, 0.01), 0.99))
}

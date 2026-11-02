#!/usr/bin/env Rscript
# Summaries for the manuscript, calculated from newly fitted replicate results.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), ".."))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[1] %in% c("simulation", "jtpa"))
    stop("Usage: Rscript Scripts/03_summarize_results.R simulation|jtpa")
dest <- file.path(root, "Results/summary")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
write <- function(x, name) write.csv(x, file.path(dest, paste0(name, ".csv")), row.names = FALSE)
stat <- function(x) c(mean = mean(x), sd = sd(x), se = sd(x) / sqrt(length(x)), min = min(x), max = max(x))

if (args[1] == "simulation") {
    cfg <- data.frame(scenario = rep(c("S1", "S2", "S3", "S4"), each = 2),
        method = rep(c("CQR_endpoint", "APS_classification"), 4),
        rho = 0, mode = "estimated_ps", gamma = 0)
    baseline <- cfg[6, ]
    rho <- baseline[rep(1, 4), ]; rho$rho <- c(-1, -0.5, 0.5, 1)
    ps <- baseline[rep(1, 3), ]; ps$mode <- c("true_ps", "misspecified_ps", "oracle")
    gamma <- baseline[rep(1, 5), ]; gamma$gamma <- (1:5) / 10
    cfg <- rbind(cfg, rho, ps, gamma)
    rownames(cfg) <- NULL
    cells <- lapply(seq_len(nrow(cfg)), function(i) {
        e <- cfg[i, ]
        path <- file.path(root, "Results", e$scenario, paste0("rho_", e$rho), e$method, e$mode)
        if (e$method == "APS_classification" && e$mode == "estimated_ps")
            path <- file.path(path, paste0("gamma_", e$gamma))
        rows <- lapply(0:99, function(seed) {
            env <- new.env()
            load(file.path(path, sprintf("result_%03d.RData", seed)), env)
            x <- env$results_df
            stopifnot(nrow(x) == 1L, x$seed == seed, x$scenario == e$scenario,
                      x$rho == e$rho, is.finite(x$coverage), is.finite(x$mean_length))
            x[c("seed", "coverage", "mean_length")]
        })
        do.call(rbind, rows)
    })
    summary <- cbind(cfg, do.call(rbind, lapply(cells, function(x) data.frame(
        coverage = mean(x$coverage), coverage_se = sd(x$coverage) / 10,
        mean_length = mean(x$mean_length), length_se = sd(x$mean_length) / 10))))
    write(summary[1:8, ], "table1_main")
    paired <- do.call(rbind, lapply(seq(1, 8, 2), function(i) {
        cqr <- cells[[i]]$mean_length; mls <- cells[[i + 1]]$mean_length
        stopifnot(all(cqr > 0))
        r <- 100 * (cqr - mls) / cqr
        ci <- t.test(r, conf.level = 0.95)$conf.int
        data.frame(scenario = cfg$scenario[i], mean = mean(r), lower = ci[1], upper = ci[2],
            p_value = t.test(r, alternative = "greater", mu = 0)$p.value)
    }))
    write(paired, "table1_paired_relative_reduction")
    write(summary[c(9, 10, 6, 11, 12), ], "table2_rho")
    write(summary[c(6, 16:20), ], "table3_gamma")
    write(summary[c(6, 13:15), ], "propensity_oracle")
} else {
    methods <- c("CQR_endpoint", "APS_classification")
    rows <- lapply(methods, function(method) do.call(rbind, lapply(1:10, function(seed) {
        x <- readRDS(file.path(root, "Results/JTPA", method, sprintf("result_%03d.rds", seed)))$row
        stopifnot(x$seed == seed, x$method == method)
        x
    })))
    b <- rows[[1]]; a <- rows[[2]]
    metrics <- c("stage1_coverage_control", "stage1_coverage_treated", "inclusion_rate",
                 "mean_length", "prop_disconnected")
    stopifnot(all(is.finite(as.matrix(b[metrics]))), all(is.finite(as.matrix(a[metrics]))))
    main <- do.call(rbind, lapply(metrics, function(metric) {
        scale <- if (metric == "mean_length") 1 else 100
        data.frame(metric = metric, cqr_mean = mean(b[[metric]]) * scale,
            cqr_sd = sd(b[[metric]]) * scale, membership_mean = mean(a[[metric]]) * scale,
            membership_sd = sd(a[[metric]]) * scale)
    }))
    write(main, "table4_JTPA")
    diagnostics <- data.frame(seed = 1:10,
        cqr_control = 100 * b$stage1_coverage_control, membership_control = 100 * a$stage1_coverage_control,
        cqr_treated = 100 * b$stage1_coverage_treated, membership_treated = 100 * a$stage1_coverage_treated,
        cqr_inclusion = 100 * b$inclusion_rate, membership_inclusion = 100 * a$inclusion_rate)
    lengths <- data.frame(seed = 1:10, cqr_length = b$mean_length, membership_length = a$mean_length,
        reduction = b$mean_length - a$mean_length, disconnected_percent = 100 * a$prop_disconnected)
    for (name in c("tableC2", "tableC3")) {
        x <- if (name == "tableC2") diagnostics else lengths
        write(x, paste0(name, "_repetitions"))
        write(data.frame(statistic = names(stat(x[[2]])), lapply(x[-1], stat)), paste0(name, "_summary"))
    }
    write(data.frame(mean_relative_reduction_percent = mean(100 * lengths$reduction / lengths$cqr_length),
        repetitions_membership_shorter = sum(lengths$reduction > 0)), "JTPA_relative_reduction")
}
cat("Saved manuscript summaries to ", dest, "\n", sep = "")

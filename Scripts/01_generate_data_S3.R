#!/usr/bin/env Rscript
################################################################################
# Script: 01_generate_data_S3.R
# Purpose: Generate S3 moderate-dimensional observational datasets
################################################################################

get_script_dir <- function() {
    file_arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
    normalizePath(getwd())
}

SCRIPT_DIR <- get_script_dir()
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))

parse_args <- function(args) {
    config <- list(
        rho = 0,
        ntrain = 20000,
        ntest = 5000,
        n_simulations = 100,
        p = 6,
        overlap_strength = 1.0
    )
    i <- 1
    while (i <= length(args)) {
        arg <- args[i]
        if (arg == "--rho" && i < length(args)) {
            config$rho <- as.numeric(args[i + 1]); i <- i + 2
        } else if (arg == "--ntrain" && i < length(args)) {
            config$ntrain <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--ntest" && i < length(args)) {
            config$ntest <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--n-sim" && i < length(args)) {
            config$n_simulations <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--p" && i < length(args)) {
            config$p <- as.integer(args[i + 1]); i <- i + 2
        } else if (arg == "--overlap-strength" && i < length(args)) {
            config$overlap_strength <- as.numeric(args[i + 1]); i <- i + 2
        } else {
            stop(sprintf("Unknown or incomplete argument: %s", arg))
        }
    }
    config
}

format_rho_label <- function(rho) paste0("rho_", format(rho, trim = TRUE, scientific = FALSE))

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(file.path(PROJECT_ROOT, "Models"))
source("02_dgp_registry.R")

config <- parse_args(commandArgs(trailingOnly = TRUE))
output_dir <- file.path(PROJECT_ROOT, "Data", "S3", format_rho_label(config$rho))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("=== Generate S3 data ===\n")
for (seed in 0:(config$n_simulations - 1)) {
    data_simul <- generate_dataset(
        "S3", config$ntrain, config$ntest, seed,
        rho = config$rho,
        p = config$p,
        overlap_strength = config$overlap_strength
    )
    save(data_simul, file = file.path(output_dir, sprintf("synthetic_%03d.RData", seed)))
    cat(sprintf("  saved seed %03d\n", seed))
}

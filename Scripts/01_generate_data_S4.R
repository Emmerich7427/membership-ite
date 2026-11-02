#!/usr/bin/env Rscript
################################################################################
# Script: 01_generate_data_S4.R
# Purpose: Generate S4 NLSM-based X-dependent bimodal datasets
################################################################################

get_script_dir <- function() {
    file_arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(file_arg) > 0) {
        return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
    }
    normalizePath(getwd())
}


SCRIPT_DIR <- get_script_dir()
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))


parse_args <- function(args) {
    config <- list(
        ntrain = 2000,
        ntest = 5000,
        n_simulations = 100,
        dgp_path = file.path(PROJECT_ROOT, "Data", "NLSM_simul.RData"),
        q_min = 0.20,
        q_max = 0.80,
        gap_multiplier = 2.60,
        gap_x_strength = 0.35,
        within_ratio = 0.15
    )

    i <- 1
    while (i <= length(args)) {
        arg <- args[i]

        if (arg == "--ntrain" && i < length(args)) {
            config$ntrain <- as.integer(args[i + 1])
            i <- i + 2
        } else if (arg == "--ntest" && i < length(args)) {
            config$ntest <- as.integer(args[i + 1])
            i <- i + 2
        } else if (arg == "--n-sim" && i < length(args)) {
            config$n_simulations <- as.integer(args[i + 1])
            i <- i + 2
        } else if (arg == "--dgp-path" && i < length(args)) {
            config$dgp_path <- normalizePath(
                args[i + 1],
                winslash = "/",
                mustWork = TRUE
            )
            i <- i + 2
        } else if (arg == "--q-min" && i < length(args)) {
            config$q_min <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--q-max" && i < length(args)) {
            config$q_max <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--gap-multiplier" && i < length(args)) {
            config$gap_multiplier <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--gap-x-strength" && i < length(args)) {
            config$gap_x_strength <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (arg == "--within-ratio" && i < length(args)) {
            config$within_ratio <- as.numeric(args[i + 1])
            i <- i + 2
        } else {
            stop(sprintf("Unknown or incomplete argument: %s", arg))
        }
    }

    config
}


format_rho_label <- function(rho) {
    paste0("rho_", format(rho, trim = TRUE, scientific = FALSE))
}


old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)

setwd(file.path(PROJECT_ROOT, "Models"))

# The registry supplies build_data_simul() and the common S1-S3 interface.
source("02_dgp_registry.R")
source("02_dgp_S4.R")


config <- parse_args(commandArgs(trailingOnly = TRUE))

if (!file.exists(config$dgp_path)) {
    stop(
        "Cannot find NLSM_simul.RData at:\n  ",
        config$dgp_path,
        "\nPlace the Lei-Candes DGP file there or use --dgp-path."
    )
}

if (config$ntrain < 1L || config$ntest < 1L) {
    stop("ntrain and ntest must be positive integers.")
}

if (config$n_simulations < 1L) {
    stop("n_simulations must be a positive integer.")
}

# S4 keeps the S4 convention and fixes rho = 0.
rho <- 0

output_dir <- file.path(
    PROJECT_ROOT,
    "Data",
    "S4",
    format_rho_label(rho)
)

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}


required_fields <- c(
    "Xtrain",
    "Ttrain",
    "Ytrain",
    "Xtest",
    "Ttest",
    "Ytest",
    "Y0test",
    "Y1test",
    "tautest",
    "CATEtest",
    "stdtest",
    "pstrain",
    "pstest"
)


cat("=== Generate S4 NLSM-based X-dependent bimodal data ===\n")
cat(sprintf("DGP file: %s\n", config$dgp_path))
cat(sprintf("ntrain: %d\n", config$ntrain))
cat(sprintf("ntest: %d\n", config$ntest))
cat(sprintf("replications: %d\n", config$n_simulations))
cat(sprintf("q_min: %.4f\n", config$q_min))
cat(sprintf("q_max: %.4f\n", config$q_max))
cat(sprintf("gap_multiplier: %.4f\n", config$gap_multiplier))
cat(sprintf("gap_x_strength: %.4f\n", config$gap_x_strength))
cat(sprintf("within_ratio: %.4f\n", config$within_ratio))
cat(sprintf("output: %s\n\n", output_dir))


for (seed in 0:(config$n_simulations - 1L)) {
    data_simul <- generate_data_S4(
        ntrain = config$ntrain,
        ntest = config$ntest,
        seed = seed,
        rho = rho,
        dgp_path = config$dgp_path,
        q_min = config$q_min,
        q_max = config$q_max,
        gap_multiplier = config$gap_multiplier,
        gap_x_strength = config$gap_x_strength,
        within_ratio = config$within_ratio
    )

    if (!identical(names(data_simul), required_fields)) {
        stop(
            sprintf(
                paste0(
                    "S4 output structure mismatch at seed %03d.\n",
                    "Expected: %s\n",
                    "Observed: %s"
                ),
                seed,
                paste(required_fields, collapse = ", "),
                paste(names(data_simul), collapse = ", ")
            )
        )
    }

    if (nrow(data_simul$Xtrain) != config$ntrain ||
        nrow(data_simul$Xtest) != config$ntest) {
        stop(sprintf("Unexpected S4 sample size at seed %03d.", seed))
    }

    save(
        data_simul,
        file = file.path(
            output_dir,
            sprintf("synthetic_%03d.RData", seed)
        )
    )

    cat(sprintf("  saved seed %03d\n", seed))
}


cat("\n=== S4 data generation complete ===\n")

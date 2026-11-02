#!/usr/bin/env Rscript
# Reconstruct the ten final-analysis datasets from the original Stata files.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "Models/02_data_JTPA.R"))
dest <- file.path(root, "Data/JTPA")
pool <- load_jtpa_pool(file.path(dest, "raw/expbif.dta"),
                       file.path(dest, "raw/earns.dta"),
                       sample_group = "full", outcome = "30m", y_scale = 1000)
stopifnot(nrow(pool$X) == 12136L)
for (seed in 1:10) {
    set.seed(seed)
    ids <- sample.int(nrow(pool$X), 12000, replace = FALSE)
    test_ids <- pool$ids[ids[7000 + seq_len(5000)]]
    data_simul <- build_data_real(pool$X, pool$T, pool$Y, 7000, 5000, seed)
    save(data_simul, test_ids, file = file.path(dest, sprintf("real_%03d.RData", seed)))
}
cat("Saved JTPA datasets for seeds 1..10: 7000 development and 5000 test subjects.\n")

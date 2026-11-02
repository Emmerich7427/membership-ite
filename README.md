# Conformal prediction for individual treatment effects via membership level sets

Minimal code and inputs for reproducing the experiments in the paper. Simulations and the JTPA analysis use R. Figure 1 uses the original numerical calculation and plotting scripts in Python.
The main experiments use `rho = 0`, the runners default to `estimated_ps`, and the fourth simulation scenario is named `S4` throughout.

## Files and dependencies

- `Models/`: data generation, propensity estimation, the two-stage methods, and JTPA diagnostics.
- `Scripts/01_generate_data_*.R`: simulation data generation and the ten JTPA data splits.
- `Scripts/02_run_*.R`: Membership level set and CQR experiment runners.
- `Scripts/03_summarize_results.R`: manuscript summaries calculated from replicate results.
- `Data/NLSM_simul.RData`: the fixed empirical-pool object used in S4.
- `Data/JTPA/raw/`: the original JTPA baseline and earnings files, `expbif.dta` and `earns.dta`.
- `oracle_sweep/`: numerical calculations and plotting for Figure 1.

Simulation data, JTPA splits, `Results/`, and `Figures/` are created when the scripts are run. Saved paper results, manuscript copies, and migration utilities are not included.

The experiments use **R 4.5.2** and the following package versions. Use separate R package environments for the simulations and JTPA analysis.

| Package | Simulations | Final JTPA fits |
| --- | --- | --- |
| dbarts | 0.9-32 | 0.9-32 |
| gbm | 2.2.2 | 2.2.3 |
| grf | 2.5.0 | 2.6.1 |
| MASS | 7.3-65 | 7.3-65 |

The final JTPA results also record Matrix 1.7-4, survival 3.8-6, Rcpp 1.1.1, and lattice 0.22-9.
Reading the Stata files additionally requires haven; data reconstruction has been checked with version 2.5.5. The scripts do not install or replace packages automatically.

## Simulation experiments

In an **R console**, set the working directory to this repository and run the blocks below.
The `run()` helper is shorthand for calling `Rscript`; no additional script file is needed.

```r
setwd("/path/to/Biometrika_Code")  # Replace with the local repository path.
rscript <- file.path(R.home("bin"), "Rscript")
run <- function(script, args = character()) {
    status <- system2(rscript, c(shQuote(script), args))
    if (status != 0) stop("Script failed: ", script)
}
```

Generate all datasets for `seed = 0,...,99`. The training sample sizes for S1, S2, S3, and S4 are
5000, 20000, 20000, and 2000, respectively; the test sample size is 5000 throughout.
S4 uses `q_min=0.2`, `q_max=0.8`, `gap_multiplier=2.6`, `gap_x_strength=0.35`, and `within_ratio=0.15`.

```r
for (s in c("S1", "S2", "S3"))
    run(paste0("Scripts/01_generate_data_", s, ".R"), c("--rho", "0"))
run("Scripts/01_generate_data_S4.R")
for (rho in c(-1, -0.5, 0.5, 1))
    run("Scripts/01_generate_data_S3.R", c("--rho", as.character(rho)))
```

The main comparison has eight experiment cells: four scenarios and two methods. The first stage of the Membership method uses BART for S1/S4 and a DP mixture density for S2/S3.

```r
for (s in c("S1", "S2", "S3", "S4")) {
    run("Scripts/02_run_cqr_endpoint.R",
        c("--scenario", s, "--rho", "0", "--ps-mode", "estimated_ps", "--n-sim", "100"))
    run("Scripts/02_run_aps_classification.R",
        c("--scenario", s, "--rho", "0", "--stage1-mode", "estimated_ps",
          "--stage1-model", if (s %in% c("S1", "S4")) "bart" else "dp", "--n-sim", "100"))
}
```

The remaining twelve cells use S3. They share the `rho=0`, `gamma=0`, estimated-propensity Membership baseline from the main comparison.

```r
# Sensitivity to rho: Table 2.
for (rho in c(-1, -0.5, 0.5, 1))
    run("Scripts/02_run_aps_classification.R",
        c("--scenario", "S3", "--rho", as.character(rho), "--n-sim", "100"))

# Propensity-score and oracle comparisons: Section 6.2.
for (mode in c("true_ps", "misspecified_ps", "oracle"))
    run("Scripts/02_run_aps_classification.R",
        c("--scenario", "S3", "--rho", "0", "--stage1-mode", mode, "--n-sim", "100"))

# Relaxed scores indexed by gamma: Table 3.
for (gamma in (1:5) / 10)
    run("Scripts/02_run_aps_classification.R",
        c("--scenario", "S3", "--rho", "0", "--gamma_score", as.character(gamma), "--n-sim", "100"))

run("Scripts/03_summarize_results.R", "simulation")
```

These commands use the paper settings in the original runners: `alpha=0.1`, allocated as 0.05 to each stage;
the outer training split has proportion 0.5, and the fitting proportion within each stage is 0.75.
Membership uses `stage2_method=mixture`, `strategy=region_uniform`, and `n_samples=100`.
The first-stage grid settings are 500/0.05/2000 and the second-stage settings are 1000/0.02/5000
(number of grid points, resolution, and maximum number of grid points, respectively).
The DP sampler uses 250 iterations, a burn-in of 125, and thinning interval 5.

For a single replicate, replace `--n-sim 100` with `--seed 0`, choosing any seed from 0 to 99.
Rerunning the same configuration and seed overwrites the corresponding output in this repository.
The full simulation study requires 2000 fits and can take substantial time.
Summaries require 100 valid results per cell and produce CSV files for Tables 1-3 and the propensity/oracle paragraph.
Parenthetical simulation statistics are Monte Carlo standard errors. For Table 1, relative reductions are calculated within matched replicates before computing the mean, t confidence interval, and one-sided p-value.

## JTPA

Switch to the JTPA package environment listed above. In an **R console** at the repository root, define the same `run()` helper and execute:

```r
run("Scripts/01_generate_data_JTPA.R")
for (seed in 1:10) {
    run("Scripts/02_run_JTPA.R", c("membership", as.character(seed)))
    run("Scripts/02_run_JTPA.R", c("cqr", as.character(seed)))
}
run("Scripts/03_summarize_results.R", "jtpa")
```

The analysis uses seven original covariates. Baseline and earnings records are matched by `recid`, leaving 12136 complete observations.
The outcome is total earnings over 30 months divided by 1000. Each repetition uses 7000 development subjects and 5000 test subjects, with shared splits for the two methods and seeds 1-10.
Both methods use estimated propensity scores. The Membership first stage uses a DP mixture, with `gamma=0` and the grid and DP settings given above.

`Models/03_outcome_JTPA.R` contains the three quantile/prediction functions used for the final JTPA results.
The JTPA runner loads them after the shared model modules. `Models/06_JTPA_diagnostics.R` contains the corresponding test-set calculations.
This reuses the common model components while preserving the rules used for the final JTPA analysis.

Results are saved in `Results/JTPA/`; statistics for Table 4 and Supplementary Tables C2-C3 are written to `Results/summary/`.
Parenthetical values in Table 4 are standard deviations across the ten splits. Supplementary summaries include individual repetitions and the mean, SD, SE, minimum, and maximum.
Coverage, inclusion, and disconnected-set proportions are reported as percentages. Supplementary Table C1 describes the variables and requires no numerical calculation.

## Figure 1

Only this part requires Python. The scripts have been checked with Python 3.12, NumPy 2.4.4, pandas 3.0.2,
SciPy 1.17.1, and matplotlib 3.10.8. These versions describe the environment used to validate the repository.

Run the following command in a **terminal** at the repository root:

```sh
python oracle_sweep/oracle_sweep.py --alpha-cf 0.025 --alpha-inc 0.025 --mu-step 0.01 --output-dir Figures
```

The command evaluates 601 points from 2 to 8 and saves Figure 1 as PDF and PNG, together with numerical results and transition locations.
The two Python files retain the original numerical calculation and plotting implementations.

## Input data

S4 uses the fixed `NLSM_simul.RData` object containing X, EY0, Etau, IQR0, IQR1, and ps.
The code generates bimodal residuals as described in the supplement; it does not refit this input object.
The upstream NLSM materials and preparation script are described in the [Lei and Candes paper repository](https://github.com/lihualei71/cfcausalPaper).
Use the fixed input accompanying this repository when reproducing the reported experiments.

The two Stata files are from the [National JTPA Study public-use data](https://www.upjohn.org/data-tools/employment-research-data-center/national-jtpa-study).
The required files are `expbif.dta` and `earns.dta`, stored in `Data/JTPA/raw/`.
See the paper's references and Supplementary Sections B1 and C1 for the data sources and analysis definitions.

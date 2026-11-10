# Simulation study: does zinbreg() recover known parameters? (handout
# Section 5.1)
#
# What it does:
# * simulates R data sets of size n from a ZINB model with known parameters,
#   in a fixed design (the covariates are drawn once; only the response is
#   redrawn): count part log mu = 1 + 0.5 x1 - 0.5 x2, zero part
#   logit pi = -0.8 + 0.6 x1 - 0.5 z, theta = 1.5;
# * fits each with zinbreg(y ~ x1 + x2 | x1 + z) at its default settings;
# * records every estimate, standard error and 95% Wald interval, the fit
#   diagnostics, and a class for each replication (usable, boundary, invalid
#   covariance, invalid values, not converged, error), with its reason;
# * computes the profile-likelihood interval for theta on the first
#   `profile-R` replications;
# * summarises bias, Monte Carlo SD, mean estimated SE, their ratio, RMSE
#   and coverage, each with a Monte Carlo standard error;
# * records the provenance (the code fingerprint and the git commit);
# * saves everything, also when the run fails.
#
# Seeds: the covariates are drawn with set.seed(seed) and replication r with
# set.seed(seed + r), so any single replication can be reproduced alone.
#
# How to run it, from the repository root, with a clean working tree:
#   devtools::install()        # install the code being studied
#   Rscript inst/benchmark/simulation.R
# Options: --R=N (data sets, default 1000), --n=N (size, default 500),
# --profile-R=N (profile intervals, default 200), --B=N (bootstrap
# resamples for the Monte Carlo SEs, default 2000), --seed=N (default 651),
# --out=DIR (default inst/benchmark).
#
# Output (in --out): simulation_results.rds (everything), simulation_summary.csv,
# simulation_replicates.csv (one row per replication) and
# simulation_status.txt, whose first line is COMPLETE or FAILED. The run is
# complete when its provenance is verified and every replication was
# attempted; fits that fail are results, counted and reported. A failed run
# still saves and then exits with status 1.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1L]) else default
}
arg_int <- function(name, default, min) {
  v <- suppressWarnings(as.integer(arg_value(name, as.character(default))))
  if (is.na(v) || v < min) stop("--", name, " must be a whole number of at least ", min)
  v
}
R <- arg_int("R", 1000L, 1L)
n <- arg_int("n", 500L, 20L)
profile_R <- arg_int("profile-R", 200L, 0L)
B <- arg_int("B", 2000L, 0L)
seed <- arg_int("seed", 651L, 0L)
out_dir <- arg_value("out", file.path("inst", "benchmark"))

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1, 1]),
               "nbzinbreg")) {
  stop("Run this script from the root of the nbzinbreg source tree.")
}
if (!requireNamespace("nbzinbreg", quietly = TRUE)) {
  stop("Install nbzinbreg first (devtools::install()).")
}
suppressPackageStartupMessages(library(nbzinbreg))
source(system.file("benchmark", "simulation_functions.R", package = "nbzinbreg"))

results <- run_simulation(R = R, n = n, seed = seed,
                          profile_R = min(profile_R, R), B = B,
                          source_path = ".")
results <- save_simulation(results, out_dir)

cat("\n", readLines(file.path(out_dir, "simulation_status.txt")), sep = "\n")
cat("\n== Summary ==\n")
print(results$summary[, c("parameter", "truth", "mean", "bias", "bias_mcse",
                          "mc_sd", "mean_se", "se_ratio", "rmse", "coverage",
                          "coverage_mcse")],
      digits = 4, row.names = FALSE)
if (nrow(results$classes$non_usable)) {
  cat("\n== Non-usable replications ==\n")
  print(results$classes$non_usable, row.names = FALSE)
}

if (results$status != "complete") {
  cat("\nSIMULATION FAILED: the results were saved in ", out_dir,
      ", marked FAILED.\n", sep = "")
  quit(status = 1)
}
cat("\nSimulation complete; results saved in ", out_dir, ".\n", sep = "")

# Benchmark: nbzinbreg against the reference implementations (handout
# Section 5).
#
# What it does:
# * fits NB with nbreg() and MASS::glm.nb(), and ZINB with zinbreg() (both
#   Hessian methods) and pscl::zeroinfl(dist = "negbin"), on
#   pscl::bioChemists and AER::NMES1988 (our functions at their default
#   settings, the references at tight settings);
# * checks that every pair is valid (both fits converged, finite estimates
#   and standard errors, the same observations and parameterization, a
#   valid covariance matrix);
# * compares the estimates side by side: the four required quantities of
#   handout Table 4 against their tolerances, and the supplementary
#   quantities (SE(theta), SE(log theta), the numeric Hessian) for
#   information;
# * times every fit with bench::mark() (10 runs each by default);
# * records the provenance: a fingerprint of the installed package's code,
#   checked against the source tree, and the git commit;
# * saves everything, also when a check fails.
#
# How to run it, from the repository root, with a clean working tree:
#   devtools::install()        # install the code being benchmarked
#   Rscript inst/benchmark/benchmark.R
# Options: --out=DIR (where to save; default inst/benchmark) and
# --iterations=N (timing runs per fit; default 10, 0 skips the timing).
#
# Output (in --out): benchmark_results.rds (everything), the readable copies
# benchmark_required.csv, benchmark_supplementary.csv and
# benchmark_timing.csv, and benchmark_status.txt, whose first line is PASSED
# or FAILED. The run passes only if every pair is valid, every required
# quantity is within its tolerance, and the provenance is verified. A failed
# run still saves its results (marked FAILED) and then exits with status 1,
# so it cannot be mistaken for a passed one. Only a passed run is committed
# as the benchmark of record.
#
# The shared functions are in inst/benchmark/benchmark_functions.R; the
# report sources the same file from the installed package.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1L]) else default
}
out_dir <- arg_value("out", file.path("inst", "benchmark"))
iterations <- as.integer(arg_value("iterations", "10"))
if (is.na(iterations) || iterations < 0L) {
  stop("--iterations must be a whole number of at least 0")
}

# Run from the repository root of nbzinbreg.
if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1, 1]),
               "nbzinbreg")) {
  stop("Run this script from the root of the nbzinbreg source tree.")
}
for (p in c("nbzinbreg", "MASS", "pscl", "AER", "bench")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Package ", p, " is needed for the benchmark; install it first.")
  }
}
suppressPackageStartupMessages(library(nbzinbreg))
source(system.file("benchmark", "benchmark_functions.R",
                   package = "nbzinbreg"))

results <- run_benchmark(source_path = ".", iterations = iterations)
save_benchmark(results, out_dir)

# Report.
show <- function(title, x) {
  cat("\n== ", title, " ==\n", sep = "")
  print(x, row.names = FALSE, digits = 4)
}
for (id in names(results$validity)) {
  show(paste("Validity:", id), results$validity[[id]])
  if (!is.null(results$tables[[id]])) {
    show(paste("Side by side:", id), results$tables[[id]])
  }
}
show("Required quantities (handout Table 4)", results$required)
show("Supplementary quantities (not graded)", results$supplementary)
if (!is.null(results$timing)) show("Timing (seconds)", results$timing)
cat("\n", readLines(file.path(out_dir, "benchmark_status.txt")), sep = "\n")

if (results$status != "passed") {
  cat("\nBENCHMARK FAILED: the results were saved for diagnosis in ",
      out_dir, ", marked FAILED.\n", sep = "")
  quit(status = 1)
}
cat("\nBenchmark passed; results saved in ", out_dir, ".\n", sep = "")

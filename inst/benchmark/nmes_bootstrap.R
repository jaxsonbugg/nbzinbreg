# Parametric bootstrap likelihood-ratio test of NB against ZINB on NMES1988
# (the report's data analysis), run once with saved results because it is
# slow (about 24 s per replicate, so about 7 hours for B = 999).
#
# What it does:
# * fits the NB model and the ZINB model of the report (count part hospital +
#   health + chronic + gender + school + insurance; zero part hospital +
#   chronic + gender + school + insurance);
# * runs vuong_test(B = ...) in chunks of --chunk replicates, saving a
#   checkpoint after every chunk: the replicates so far, the failures, the
#   warnings and the state of the random number stream;
# * pools the chunks into one bootstrap and computes the p-value and its
#   bounds with the same formulas as vuong_test();
# * records the provenance (the code fingerprint and the git commit) and
#   saves everything.
#
# Resuming. If the run stops (the computer sleeps, a crash), run the same
# command again: it reads the checkpoint, checks that the settings, the
# installed code and the fitted models are unchanged (it refuses otherwise),
# restores the random number stream and continues with the next chunk. The
# finished replicates are kept unchanged, and the result is identical to an
# uninterrupted run with the same seed.
#
# How to run it, from the repository root, with a clean working tree:
#   devtools::install()
#   Rscript inst/benchmark/nmes_bootstrap.R
# Options: --B=N (default 999), --chunk=N (default 25), --seed=N (default
# 651), --alpha=X (default 0.05), --maxit=N (the ZINB EM limit, default
# 5000), --out=DIR (default inst/benchmark), --no-resume (refuse an
# existing checkpoint), --max-chunks=N (stop after N new chunks, for pilots),
# --pilot-default=FILE and --pilot-high=FILE (the results files of the
# pilots: maxit is then chosen by the fixed rule in nmes_maxit_decision(),
# which uses only the pilots' failure shares).
#
# Output (in --out): nmes_bootstrap_checkpoint.rds (kept for resuming),
# nmes_bootstrap_results.rds, nmes_bootstrap_replicates.csv and
# nmes_bootstrap_status.txt. Its first line is COMPLETE or FAILED (did the
# run finish, with verified provenance?); its second line is CONCLUSIVE,
# INCONCLUSIVE (with the reasons) or NOT EVALUATED (can the test decide?).
# An inconclusive test is a valid result of a complete run and exits with
# status 0; a failed run exits with status 1.

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
B <- arg_int("B", 999L, 1L)
chunk <- arg_int("chunk", 25L, 1L)
seed <- arg_int("seed", 651L, 0L)
alpha <- as.numeric(arg_value("alpha", "0.05"))
out_dir <- arg_value("out", file.path("inst", "benchmark"))
max_chunks <- if (any(grepl("^--max-chunks=", args))) arg_int("max-chunks", NA, 1L) else Inf
resume <- !"--no-resume" %in% args
pilot_default <- arg_value("pilot-default", NA)
pilot_high <- arg_value("pilot-high", NA)
maxit_given <- any(grepl("^--maxit=", args))

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1, 1]),
               "nbzinbreg")) {
  stop("Run this script from the root of the nbzinbreg source tree.")
}
if (!requireNamespace("nbzinbreg", quietly = TRUE)) {
  stop("Install nbzinbreg first (devtools::install()).")
}
suppressPackageStartupMessages(library(nbzinbreg))
source(system.file("benchmark", "nmes_functions.R", package = "nbzinbreg"))

# The EM iteration limit: from the pilots by the fixed rule, or as given.
decision <- NULL
if (!is.na(pilot_default)) {
  decision <- nmes_maxit_decision(
    readRDS(pilot_default),
    if (!is.na(pilot_high)) readRDS(pilot_high))
  if (decision$needs_high_pilot) {
    stop("The default-maxit pilot's failure share is ",
         format(round(decision$share_default, 3)), " (above ",
         decision$threshold, "); run a pilot with --maxit=20000 and pass it ",
         "with --pilot-high.")
  }
  if (maxit_given && arg_int("maxit", 5000L, 1L) != decision$maxit) {
    stop("--maxit disagrees with the rule's choice (", decision$maxit, ").")
  }
  maxit <- decision$maxit
  message(decision$summary)
} else {
  maxit <- arg_int("maxit", 5000L, 1L)
}

results <- run_nmes_bootstrap(B = B, chunk = chunk, seed = seed,
                              alpha = alpha, maxit = maxit, out_dir = out_dir,
                              source_path = ".", resume = resume,
                              max_chunks = max_chunks,
                              maxit_decision = decision)
results <- save_nmes_bootstrap(results, out_dir)

cat("\n", readLines(file.path(out_dir, "nmes_bootstrap_status.txt")), sep = "\n")
if (results$run_status != "COMPLETE") {
  cat("\nBOOTSTRAP RUN NOT COMPLETE: the results were saved in ", out_dir,
      ", marked FAILED. If chunks remain, re-run the same command to resume.\n",
      sep = "")
  quit(status = 1)
}
cat("\nBootstrap complete (", results$inference, "); results saved in ",
    out_dir, ".\n", sep = "")

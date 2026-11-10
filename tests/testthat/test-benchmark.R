# The benchmark's shared functions (inst/benchmark/benchmark_functions.R):
# provenance fingerprints, validity checks, comparisons, status, saving and
# timing. The functions are sourced from the package's installed files
# (under devtools::test(), pkgload maps system.file() to inst/).

bench_env <- function() {
  path <- system.file("benchmark", "benchmark_functions.R",
                      package = "nbzinbreg")
  if (!nzchar(path)) testthat::skip("benchmark functions not available")
  env <- new.env()
  sys.source(path, envir = env, keep.source = FALSE)
  env
}

# The package source tree, available under devtools::test() but not in an
# R CMD check of the installed package.
source_root <- function() {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  if (!file.exists(file.path(root, "DESCRIPTION")) ||
      !dir.exists(file.path(root, "R"))) {
    testthat::skip("package source tree not available")
  }
  root
}

copy_source_r <- function(root) {
  dest <- file.path(tempfile("srccopy"), "R")
  dir.create(dest, recursive = TRUE)
  file.copy(list.files(file.path(root, "R"), full.names = TRUE), dest)
  dirname(dest)
}

# Standardized NB pair on the 8-observation worked example.
nb_std_pair <- function(b) {
  skip_if_not_installed("MASS")
  d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
  list(ours = b$standardize_fit(nbreg(y ~ x, data = d), "NB", "ours"),
       ref = b$standardize_fit(MASS::glm.nb(y ~ x, data = d,
                                            control = b$bench_reference_control("NB")),
                               "NB", "reference"))
}

zinb_sim <- function() {
  set.seed(5)
  x <- stats::rnorm(300)
  y <- ifelse(stats::runif(300) < stats::plogis(-0.5 + 0.8 * x), 0,
              stats::rnbinom(300, size = 2, mu = exp(1 + 0.5 * x)))
  data.frame(y = y, x = x)
}

# ---- provenance ----------------------------------------------------------------

test_that("the source fingerprint equals the loaded package's", {
  b <- bench_env()
  root <- source_root()
  src <- b$code_fingerprint_source(root)
  ins <- b$code_fingerprint_installed()
  cmp <- b$compare_fingerprints(src, ins)
  expect_true(cmp$match)
  expect_identical(src$md5, ins$md5)
  expect_identical(cmp$differ, character())
  expect_gt(length(src$objects), 150)
  expect_true(all(c("nbreg", "zinb_em", "lr_default_maxit") %in%
                    names(src$objects)))
})

test_that("a changed body, or a missing file, changes the fingerprint and is named", {
  b <- bench_env()
  root <- source_root()
  ins <- b$code_fingerprint_installed()
  # a constant changed in a copy of the source
  copy <- copy_source_r(root)
  f <- file.path(copy, "R", "vuong.R")
  txt <- readLines(f)
  txt <- sub("^lr_default_maxit <- 500L$", "lr_default_maxit <- 501L", txt)
  writeLines(txt, f)
  cmp <- b$compare_fingerprints(b$code_fingerprint_source(copy), ins)
  expect_false(cmp$match)
  expect_identical(cmp$differ, "lr_default_maxit")
  # one function body changed
  copy2 <- copy_source_r(root)
  f2 <- file.path(copy2, "R", "nbreg-fit.R")
  txt2 <- readLines(f2)
  i <- grep("^rel_change <- function", txt2)
  txt2[i] <- "rel_change <- function(new, old) abs(new - old) / (abs(old) + 0.2)"
  writeLines(txt2, f2)
  cmp2 <- b$compare_fingerprints(b$code_fingerprint_source(copy2), ins)
  expect_identical(cmp2$differ, "rel_change")
  # a file missing from the source: its objects are named
  copy3 <- copy_source_r(root)
  file.remove(file.path(copy3, "R", "diagnostics.R"))
  cmp3 <- b$compare_fingerprints(b$code_fingerprint_source(copy3), ins)
  expect_false(cmp3$match)
  expect_true(all(c("expected_counts", "rootogram") %in% cmp3$extra_installed))
})

test_that("provenance records the code, the commit and the machine", {
  b <- bench_env()
  root <- source_root()
  p <- b$provenance(root)
  expect_true(p$status %in% c("verified", "unverified"))
  expect_identical(p$status == "verified", length(p$reasons) == 0L)
  expect_identical(p$fingerprint_installed, p$fingerprint_source)
  for (nm in c("git_commit", "r_version", "platform", "os", "date",
               "package_version", "benchmark_functions_md5")) {
    expect_true(is.character(p[[nm]]) && length(p[[nm]]) == 1L, info = nm)
  }
  expect_named(p$versions, c("nbzinbreg", "MASS", "pscl", "AER", "bench"))
  # an unreadable source tree is a reason, not an error
  q <- b$provenance(file.path(tempdir(), "no-such-tree"))
  expect_identical(q$status, "unverified")
  expect_match(q$reasons, "benchmark functions differ", all = FALSE)
  # several helper files: each is compared, and its MD5 recorded
  p2 <- b$provenance(root, files = c("benchmark_functions.R",
                                     "simulation_functions.R"))
  expect_named(p2$helper_md5, c("benchmark_functions.R",
                                "simulation_functions.R"))
  expect_true(all(!is.na(p2$helper_md5)))
  q2 <- b$provenance(file.path(tempdir(), "no-such-tree"),
                     files = c("benchmark_functions.R", "simulation_functions.R"))
  expect_match(q2$reasons, "simulation_functions.R differs", all = FALSE)
  q3 <- b$provenance(root, files = "no_such_helper.R")
  expect_match(q3$reasons, "no_such_helper.R differs", all = FALSE)
})

test_that("saved provenance is checked against the installed code", {
  b <- bench_env()
  now <- b$code_fingerprint_installed()
  saved <- list(r_version = R.version.string, fingerprint_installed = now$md5,
                fingerprint_objects = now$objects, git_commit = "abc123")
  expect_identical(b$check_saved_provenance(saved)$status, "matches")
  changed <- saved
  changed$fingerprint_installed <- "0000"
  changed$fingerprint_objects[["nbreg"]] <- "0000"
  r <- b$check_saved_provenance(changed)
  expect_identical(r$status, "differs")
  expect_match(r$message, "^WARNING.*differing objects: nbreg")
  other <- saved
  other$r_version <- "R version 0.0.0"
  expect_identical(b$check_saved_provenance(other)$status, "cannot verify")
})

# ---- validity checks -------------------------------------------------------------

test_that("a valid pair passes every check", {
  b <- bench_env()
  p <- nb_std_pair(b)
  chk <- b$validity_checks(p$ours, p$ref)
  expect_true(all(chk$passed))
  expect_identical(nrow(chk), 10L)
  expect_true(is.na(b$pair_failure_reason(chk)))
})

test_that("each validity failure is detected with an explicit reason", {
  b <- bench_env()
  p <- nb_std_pair(b)
  reason <- function(o, r) b$pair_failure_reason(b$validity_checks(o, r))
  o <- p$ours
  o$converged <- FALSE
  expect_match(reason(o, p$ref), "^our fit converged")
  r <- p$ref
  r$converged <- FALSE
  r$converged_detail <- "glm.nb converged = FALSE"
  expect_match(reason(p$ours, r), "^reference fit converged: glm.nb converged = FALSE")
  r <- p$ref
  r$vcov[1, 1] <- NaN
  expect_match(reason(p$ours, r), "^finite values \\(reference\\): not finite: standard errors")
  o <- p$ours
  o$loglik <- -Inf
  expect_match(reason(o, p$ref), "not finite: log-likelihood")
  o <- p$ours
  o$nobs <- 9
  expect_match(reason(o, p$ref), "^same observations: different numbers")
  o <- p$ours
  o$y[3] <- o$y[3] + 1
  expect_match(reason(o, p$ref), "^same observations: the response vectors differ")
  o <- p$ours
  o$coef <- rev(o$coef)
  expect_match(reason(o, p$ref), "^same parameterization: parameter names differ")
  o <- p$ours
  names(o$coef)[2] <- "z"
  expect_match(reason(o, p$ref), "^same parameterization")
  o <- p$ours
  o$vcov_check[1, 2] <- o$vcov_check[1, 2] * 1.5
  expect_match(reason(o, p$ref), "^valid covariance \\(ours\\): not symmetric")
  o <- p$ours
  o$vcov_check <- diag(c(1, -1))
  expect_match(reason(o, p$ref), "^valid covariance \\(ours\\): not positive definite")
  # a fit that failed outright
  chk <- b$validity_checks(NULL, p$ref, ours_error = "boom")
  expect_identical(nrow(chk), 2L)
  expect_identical(b$pair_failure_reason(chk), "our fit completed: error: boom")
})

test_that("a ZINB covariance with an unestimated parameter is invalid", {
  skip_if_not_installed("pscl")
  b <- bench_env()
  d <- zinb_sim()
  o <- b$standardize_fit(zinbreg(y ~ x | x, d), "ZINB", "ours")
  r <- b$standardize_fit(pscl::zeroinfl(y ~ x | x, d, dist = "negbin",
                                        control = b$bench_reference_control("ZINB")),
                         "ZINB", "reference")
  expect_true(all(b$validity_checks(o, r)$passed))
  # the full covariance (with log(theta)) is the one checked
  expect_identical(rownames(o$vcov_check)[ncol(o$vcov_check)], "Log(theta)")
  o$vcov_check["Log(theta)", ] <- NA
  expect_match(b$pair_failure_reason(b$validity_checks(o, r)),
               "^valid covariance \\(ours\\): non-finite entries in the rows of Log\\(theta\\)")
  # an unconverged gradient check fails our convergence check
  o2 <- b$standardize_fit(
    suppressWarnings(zinbreg(y ~ x | x, d, control = zinbreg_control(maxit = 2))),
    "ZINB", "ours")
  expect_match(b$pair_failure_reason(b$validity_checks(o2, r)),
               "^our fit converged: converged = FALSE")
})

# ---- comparisons -----------------------------------------------------------------

test_that("the side-by-side table and the required summary are exact", {
  b <- bench_env()
  p <- nb_std_pair(b)
  tab <- b$side_by_side(p$ours, p$ref, "NB")
  expect_identical(tab$parameter, c("(Intercept)", "x", "theta", "log-likelihood"))
  expect_equal(tab$estimate_abs_diff,
               abs(c(p$ours$coef, p$ours$theta, p$ours$loglik) -
                     c(p$ref$coef, p$ref$theta, p$ref$loglik)),
               ignore_attr = TRUE)
  expect_equal(tab$signed_diff[4], p$ours$loglik - p$ref$loglik)
  expect_equal(tab$se_ours[3], p$ours$se_theta)
  expect_true(is.na(tab$se_ours[4]))
  req <- b$required_summary(tab, "NB", "worked")
  expect_identical(req$quantity, c("coefficients", "standard errors", "theta",
                                   "log-likelihood"))
  expect_identical(req$tolerance, c(1e-5, 1e-5, 1e-4, 1e-4))
  expect_equal(req$max_abs_diff[1], max(tab$estimate_abs_diff[1:2]))
  expect_equal(req$max_abs_diff[2], max(tab$se_abs_diff[1:2]))
  expect_true(all(req$passed))
  expect_true(req$note[4] %in% c("log-likelihood: ours is higher",
                                 "log-likelihood: reference is higher",
                                 "log-likelihoods equal"))
  # supplementary: absolute and relative differences, no pass/fail
  sup <- b$supplementary_rows(tab, "worked", "NB", NA)
  expect_identical(sup$quantity, "SE(theta)")
  expect_false("passed" %in% names(sup))
  expect_equal(sup$max_rel_diff, abs(p$ours$se_theta - p$ref$se_theta) /
                 abs(p$ref$se_theta))
  expect_match(sup$note, "not a graded requirement")
})

test_that("a tolerance failure gives status failed and keeps the true numbers", {
  b <- bench_env()
  p <- nb_std_pair(b)
  tab <- b$side_by_side(p$ours, p$ref, "NB")
  tab$estimate_abs_diff[1] <- 0.5     # doctored: far beyond 1e-5
  req <- b$required_summary(tab, "NB", "worked")
  expect_false(req$passed[1])
  expect_identical(req$max_abs_diff[1], 0.5)
  cmp <- list(validity = list(worked = b$validity_checks(p$ours, p$ref)),
              required = req)
  prov <- list(status = "verified", reasons = character())
  st <- b$benchmark_status(cmp, prov)
  expect_identical(st$status, "failed")
  expect_identical(st$failures$item, "coefficients")
  expect_identical(st$failures$value, 0.5)
  expect_match(st$failures$reason, "above the tolerance")
  # the same comparison with verified provenance and no doctoring passes
  cmp$required <- b$required_summary(b$side_by_side(p$ours, p$ref, "NB"),
                                     "NB", "worked")
  expect_identical(b$benchmark_status(cmp, prov)$status, "passed")
  # unverified provenance alone fails the run
  bad <- list(status = "unverified", reasons = "uncommitted changes")
  st2 <- b$benchmark_status(cmp, bad)
  expect_identical(st2$status, "failed")
  expect_identical(st2$failures$reason, "uncommitted changes")
})

test_that("an invalid pair is not assessed and fails the run", {
  b <- bench_env()
  skip_if_not_installed("MASS")
  d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
  pair <- list(id = "broken", spec_id = "broken", data_name = "worked",
               model = "NB", hessian = NA_character_,
               ours = list(fit = NULL, error = "boom", warnings = character()),
               ref = b$capture_fit(MASS::glm.nb(y ~ x, data = d)))
  cmp <- b$compare_pairs(list(pair))
  expect_null(cmp$tables$broken)
  expect_true(all(is.na(cmp$required$max_abs_diff)))
  expect_false(any(cmp$required$passed))
  expect_match(cmp$required$note, "not assessed \\(invalid: our fit completed: error: boom",
               all = TRUE)
  st <- b$benchmark_status(cmp, list(status = "verified", reasons = character()))
  expect_identical(st$status, "failed")
})

test_that("a numeric-Hessian pair is supplementary, with analytic against numeric", {
  skip_if_not_installed("pscl")
  b <- bench_env()
  d <- zinb_sim()
  spec <- b$bench_spec("sim_ZINB", "sim", "ZINB", y ~ x | x, d)
  pairs <- list(b$fit_pair(spec, "analytic"), b$fit_pair(spec, "numeric"))
  cmp <- b$compare_pairs(pairs)
  # only the analytic pair is graded
  expect_identical(unique(cmp$required$pair), "sim_ZINB_analytic")
  expect_true(all(cmp$required$passed))
  sup <- cmp$supplementary
  expect_setequal(unique(sup$pair), c("sim_ZINB_analytic", "sim_ZINB_numeric",
                                      "sim_ZINB_analytic_vs_numeric"))
  num <- sup[sup$pair == "sim_ZINB_numeric", ]
  expect_identical(num$quantity, c("coefficients", "standard errors",
                                   "log(theta)", "log-likelihood",
                                   "SE(log theta)"))
  expect_true(all(is.finite(num$max_rel_diff)))
  avn <- sup[sup$pair == "sim_ZINB_analytic_vs_numeric", ]
  expect_lt(avn$max_abs_diff, 1e-4)
  expect_false("passed" %in% names(sup))
})

# ---- saving and timing -------------------------------------------------------------

test_that("a failed run is saved in full and marked FAILED", {
  b <- bench_env()
  p <- nb_std_pair(b)
  tab <- b$side_by_side(p$ours, p$ref, "NB")
  tab$estimate_abs_diff[1] <- 0.5
  req <- b$required_summary(tab, "NB", "worked")
  cmp <- list(validity = list(worked = b$validity_checks(p$ours, p$ref)),
              required = req)
  prov <- list(status = "verified", reasons = character(),
               fingerprint_installed = "f", fingerprint_source = "f",
               git_commit = "abc", git_clean = TRUE, package_version = "0",
               package_built = "x", r_version = R.version.string,
               platform = "p", os = "o", cpu = "c", cores = 1L,
               versions = c(nbzinbreg = "0"), date = "d")
  st <- b$benchmark_status(cmp, prov)
  res <- list(status = st$status, failures = st$failures, provenance = prov,
              tables = list(worked = tab), validity = cmp$validity,
              required = req,
              supplementary = b$supplementary_rows(tab, "worked", "NB", NA),
              timing = NULL, iterations = 0L)
  out <- tempfile("bench_out")
  b$save_benchmark(res, out)
  expect_setequal(list.files(out), c("benchmark_results.rds",
                                     "benchmark_required.csv",
                                     "benchmark_supplementary.csv",
                                     "benchmark_status.txt"))
  status <- readLines(file.path(out, "benchmark_status.txt"))
  expect_identical(status[1], "FAILED")
  expect_match(status, "worked / coefficients: maximum absolute difference above the tolerance \\(0.5 > 1e-05\\)",
               all = FALSE)
  back <- readRDS(file.path(out, "benchmark_results.rds"))
  expect_identical(back$status, "failed")
  expect_identical(back$tables$worked$estimate_abs_diff[1], 0.5)
  csv <- utils::read.csv(file.path(out, "benchmark_required.csv"))
  expect_false(csv$passed[1])
})

test_that("time_fits() times every fit of a specification", {
  skip_on_cran()
  skip_if_not_installed("bench")
  skip_if_not_installed("MASS")
  b <- bench_env()
  d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
  tm <- b$time_fits(b$bench_spec("worked", "worked", "NB", y ~ x, d),
                    iterations = 2L)
  expect_identical(tm$fit, c("nbreg() (default settings)",
                             "glm.nb() (tight settings)",
                             "glm.nb() (default settings)"))
  expect_identical(tm$runs, rep(2L, 3))
  expect_true(all(tm$min_s > 0 & tm$min_s <= tm$median_s &
                    tm$median_s <= tm$max_s))
})

test_that("the benchmark specifications are the handout's", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("AER")
  b <- bench_env()
  s <- b$bench_specs()
  expect_identical(vapply(s, `[[`, "", "id"),
                   c("bioChemists_NB", "bioChemists_ZINB", "NMES1988_NB",
                     "NMES1988_ZINB"))
  expect_identical(vapply(s, function(x) nrow(x$data), 1L),
                   c(915L, 915L, 4406L, 4406L))
  expect_identical(deparse(s[[4]]$formula[[3]][[3]]),
                   "hospital + chronic + gender + school + insurance")
  tol <- b$bench_tolerances()
  expect_identical(tol$tolerance,
                   c(1e-5, 1e-5, 1e-4, 1e-4, 1e-3, 1e-3, 1e-3, 1e-3))
  expect_identical(tol$quantity[tol$model == "ZINB"],
                   c("coefficients", "standard errors", "log(theta)",
                     "log-likelihood"))
})

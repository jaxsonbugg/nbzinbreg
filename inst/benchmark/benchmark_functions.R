# Functions shared by the benchmark script (inst/benchmark/benchmark.R), the
# report (vignettes/articles/report.qmd), the supplementary benchmark article
# and the tests (tests/testthat/test-benchmark.R).
#
# What the benchmark does (handout Section 5): it fits each model with our
# package and with the reference implementation (MASS::glm.nb() for the
# negative binomial (NB) model, pscl::zeroinfl(dist = "negbin") for the
# zero-inflated NB (ZINB) model) on two data sets, shows the two sets of
# estimates side by side with their maximum absolute differences, checks
# those against the required tolerances (handout Table 4), and times the
# fits.
#
# This file is not part of the package code in R/. The reference functions
# are called here on purpose: they are what we compare against. They are
# never called by the fitting functions, which the prohibited-call test
# checks. The file is installed with the package (inst/), so the report can
# source the same code from the installed package with
#   source(system.file("benchmark", "benchmark_functions.R",
#                      package = "nbzinbreg"))
#
# Terms used below:
# * pair: our fit and the reference fit of one model to one data set.
# * maximum absolute difference: the largest |ours - reference| over a group
#   of numbers, for example over all coefficients of one fit.
# * required quantity: one of the four rows of handout Table 4
#   (coefficients, their standard errors, theta or log theta, the
#   log-likelihood). Only these have pass/fail tolerances.
# * supplementary quantity: an additional comparison reported for
#   information only (SE(theta), SE(log theta), the numeric ZINB Hessian),
#   with absolute and relative differences and no threshold.
# * code fingerprint: an MD5 hash of the text of every object the package
#   defines, used to show that the installed code is the source code of a
#   given commit (see code_fingerprint_text()).

# ---------------------------------------------------------------------------
# Specifications and tolerances
# ---------------------------------------------------------------------------

# Required agreement (handout Table 4): maximum absolute differences.
bench_tolerances <- function() {
  data.frame(
    model = rep(c("NB", "ZINB"), each = 4),
    quantity = c("coefficients", "standard errors", "theta", "log-likelihood",
                 "coefficients", "standard errors", "log(theta)",
                 "log-likelihood"),
    tolerance = c(1e-5, 1e-5, 1e-4, 1e-4, 1e-3, 1e-3, 1e-3, 1e-3),
    stringsAsFactors = FALSE)
}

# The reference settings (handout tips): tight tolerances, so the
# references are as close to the exact maximum as they can get.
bench_reference_control <- function(model) {
  if (model == "NB") {
    stats::glm.control(epsilon = 1e-10, maxit = 100)
  } else {
    pscl::zeroinfl.control(reltol = 1e-12)
  }
}

# Load a benchmark data set from the package that provides it.
load_bench_data <- function(name) {
  env <- new.env()
  pkg <- switch(name, bioChemists = "pscl", NMES1988 = "AER",
                stop("unknown data set ", name))
  utils::data(list = name, package = pkg, envir = env)
  env[[name]]
}

# One benchmark specification: an identifier, the data, the model ("NB" or
# "ZINB") and the formula.
bench_spec <- function(id, data_name, model, formula, data) {
  list(id = id, data_name = data_name, model = model, formula = formula,
       data = data)
}

# The four required specifications (handout Section 5 and Appendix A). The
# ZINB zero part for bioChemists uses every covariate, as in Appendix A; for
# NMES1988 it is the handout's minimum set.
bench_specs <- function() {
  bc <- load_bench_data("bioChemists")
  nm <- load_bench_data("NMES1988")
  list(
    bench_spec("bioChemists_NB", "bioChemists", "NB",
               art ~ fem + mar + kid5 + phd + ment, bc),
    bench_spec("bioChemists_ZINB", "bioChemists", "ZINB",
               art ~ fem + mar + kid5 + phd + ment |
                 fem + mar + kid5 + phd + ment, bc),
    bench_spec("NMES1988_NB", "NMES1988", "NB",
               visits ~ hospital + health + chronic + gender + school +
                 insurance, nm),
    bench_spec("NMES1988_ZINB", "NMES1988", "ZINB",
               visits ~ hospital + health + chronic + gender + school +
                 insurance | hospital + chronic + gender + school +
                 insurance, nm))
}

# ---------------------------------------------------------------------------
# Fitting a pair
# ---------------------------------------------------------------------------

# Run `expr`, keeping its value, or its error message, and every warning.
capture_fit <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (inherits(value, "error")) {
    list(fit = NULL, error = conditionMessage(value), warnings = warnings)
  } else {
    list(fit = value, error = NA_character_, warnings = warnings)
  }
}

# Our fit (default settings; for ZINB with the given Hessian) and the
# reference fit (tight settings) for one specification. A fit that fails is
# recorded with its error message instead of stopping the benchmark.
fit_pair <- function(spec, hessian = c("analytic", "numeric")) {
  hessian <- match.arg(hessian)
  if (spec$model == "NB") {
    ours <- capture_fit(nbzinbreg::nbreg(spec$formula, data = spec$data))
    ref <- capture_fit(MASS::glm.nb(spec$formula, data = spec$data,
                                    control = bench_reference_control("NB")))
  } else {
    ours <- capture_fit(nbzinbreg::zinbreg(
      spec$formula, data = spec$data,
      control = nbzinbreg::zinbreg_control(hessian = hessian)))
    ref <- capture_fit(pscl::zeroinfl(spec$formula, data = spec$data,
                                      dist = "negbin",
                                      control = bench_reference_control("ZINB")))
  }
  list(id = if (spec$model == "ZINB") paste0(spec$id, "_", hessian) else spec$id,
       spec_id = spec$id, data_name = spec$data_name, model = spec$model,
       hessian = if (spec$model == "ZINB") hessian else NA_character_,
       ours = ours, ref = ref)
}

# Every quantity the comparison needs, in one standard form for either
# package. `theta` is theta for NB and log(theta) for ZINB (the scale of
# handout Table 4), `se_theta` its standard error. The ZINB coefficient
# names are zeroinfl's (count_*, zero_*). `vcov` is the covariance matrix of
# the coefficients; `vcov_check` is the matrix whose validity is checked (for
# our ZINB fit the full one, including log(theta)).
standardize_fit <- function(fit, model, source = c("ours", "reference")) {
  source <- match.arg(source)
  if (source == "ours" && model == "NB") {
    list(coef = fit$coefficients, vcov = fit$vcov, vcov_check = fit$vcov,
         theta = fit$theta, se_theta = fit$SE.theta, loglik = fit$loglik,
         nobs = fit$nobs, y = unname(as.numeric(fit$y)),
         converged = isTRUE(fit$converged),
         converged_detail = sprintf("converged = %s, %d iterations",
                                    fit$converged, fit$iter))
  } else if (source == "ours") {
    ok <- isTRUE(fit$converged) && isTRUE(fit$gradient_check$ok) &&
      length(fit$diverging) == 0L
    list(coef = stats::coef(fit), vcov = fit$vcov, vcov_check = fit$vcov_full,
         theta = log(fit$theta), se_theta = fit$SE.logtheta,
         loglik = fit$loglik, nobs = fit$nobs, y = unname(as.numeric(fit$y)),
         converged = ok,
         converged_detail = sprintf(paste(
           "converged = %s, gradient check passed = %s (Newton decrement",
           "%.3g), diverging coefficients: %s, %d EM iterations"),
           fit$converged, isTRUE(fit$gradient_check$ok),
           fit$gradient_check$decrement,
           if (length(fit$diverging)) paste(fit$diverging, collapse = ", ")
           else "none", fit$iter))
  } else if (model == "NB") {
    ok <- isTRUE(fit$converged) && is.null(fit$th.warn)
    list(coef = stats::coef(fit), vcov = stats::vcov(fit),
         vcov_check = stats::vcov(fit), theta = fit$theta,
         se_theta = fit$SE.theta, loglik = as.numeric(stats::logLik(fit)),
         nobs = stats::nobs(fit), y = unname(as.numeric(fit$y)),
         converged = ok,
         converged_detail = sprintf(
           "glm.nb converged = %s, theta warning: %s", fit$converged,
           if (is.null(fit$th.warn)) "none" else fit$th.warn))
  } else {
    ok <- identical(as.integer(fit$optim$convergence), 0L) &&
      isTRUE(fit$converged)
    list(coef = stats::coef(fit), vcov = fit$vcov, vcov_check = fit$vcov,
         theta = log(fit$theta), se_theta = fit$SE.logtheta,
         loglik = fit$loglik, nobs = fit$n, y = unname(as.numeric(fit$y)),
         converged = ok,
         converged_detail = sprintf(
           "zeroinfl optim convergence code %s, converged = %s",
           fit$optim$convergence, fit$converged))
  }
}

# ---------------------------------------------------------------------------
# Validity checks
# ---------------------------------------------------------------------------

# Is a covariance matrix usable? Finite, symmetric to 1e-8 relative to its
# largest entry, and positive definite (its Cholesky factorization exists).
# Returns list(ok, detail).
covariance_valid <- function(V) {
  if (is.null(V) || !is.matrix(V) || !length(V)) {
    return(list(ok = FALSE, detail = "no covariance matrix"))
  }
  if (any(!is.finite(V))) {
    bad <- rownames(V)[rowSums(!is.finite(V)) > 0]
    return(list(ok = FALSE, detail = paste(
      "non-finite entries in the rows of", paste(bad, collapse = ", "))))
  }
  scale <- max(abs(V))
  asym <- max(abs(V - t(V)))
  if (asym > 1e-8 * scale) {
    return(list(ok = FALSE, detail = sprintf(
      "not symmetric (largest asymmetry %.3g, relative %.3g)", asym,
      asym / scale)))
  }
  R <- tryCatch(chol((V + t(V)) / 2), error = function(e) NULL)
  if (is.null(R)) {
    return(list(ok = FALSE, detail = "not positive definite (chol() fails)"))
  }
  list(ok = TRUE, detail = sprintf("finite, symmetric, positive definite; condition number %.3g",
                                   kappa((V + t(V)) / 2, exact = TRUE)))
}

# The validity checks of one pair, given the standardized fits (or NULL for
# a fit that failed, with its error message). Returns a data frame with one
# row per check: check, passed, detail. A pair is compared only if every
# check passes; otherwise the first failed check is its failure reason.
validity_checks <- function(ours, ref, ours_error = NA_character_,
                            ref_error = NA_character_) {
  rows <- list()
  add <- function(check, passed, detail) {
    rows[[length(rows) + 1L]] <<- data.frame(
      check = check, passed = isTRUE(passed), detail = detail,
      stringsAsFactors = FALSE)
  }
  add("our fit completed", is.null(ours) == FALSE,
      if (is.null(ours)) paste("error:", ours_error) else "yes")
  add("reference fit completed", is.null(ref) == FALSE,
      if (is.null(ref)) paste("error:", ref_error) else "yes")
  if (is.null(ours) || is.null(ref)) return(do.call(rbind, rows))

  add("our fit converged", ours$converged, ours$converged_detail)
  add("reference fit converged", ref$converged, ref$converged_detail)

  finite_detail <- function(s) {
    se <- sqrt(diag(s$vcov))
    parts <- c(coefficients = all(is.finite(s$coef)),
               `standard errors` = all(is.finite(se)),
               theta = is.finite(s$theta),
               `SE(theta)` = is.finite(s$se_theta),
               `log-likelihood` = is.finite(s$loglik))
    list(ok = all(parts),
         detail = if (all(parts)) "all finite" else paste(
           "not finite:", paste(names(parts)[!parts], collapse = ", ")))
  }
  fo <- finite_detail(ours)
  fr <- finite_detail(ref)
  add("finite values (ours)", fo$ok, fo$detail)
  add("finite values (reference)", fr$ok, fr$detail)

  same_n <- identical(as.numeric(ours$nobs), as.numeric(ref$nobs))
  same_y <- length(ours$y) == length(ref$y) && all(ours$y == ref$y)
  add("same observations", same_n && same_y,
      if (!same_n) sprintf("different numbers of observations (%s and %s)",
                           ours$nobs, ref$nobs)
      else if (!same_y) "the response vectors differ"
      else sprintf("%s observations, identical responses", ours$nobs))

  same_names <- identical(names(ours$coef), names(ref$coef))
  vcov_names <- identical(rownames(ours$vcov), names(ours$coef)) &&
    identical(rownames(ref$vcov), names(ref$coef))
  add("same parameterization", same_names && vcov_names,
      if (!same_names) {
        sprintf("parameter names differ: ours %s; reference %s",
                paste(names(ours$coef), collapse = ", "),
                paste(names(ref$coef), collapse = ", "))
      } else if (!vcov_names) {
        "covariance matrix rows do not match the coefficient names"
      } else {
        sprintf("%d coefficients with identical names in the same order",
                length(ours$coef))
      })

  vo <- covariance_valid(ours$vcov_check)
  vr <- covariance_valid(ref$vcov_check)
  add("valid covariance (ours)", vo$ok, vo$detail)
  add("valid covariance (reference)", vr$ok, vr$detail)
  do.call(rbind, rows)
}

pair_valid <- function(checks) all(checks$passed)

pair_failure_reason <- function(checks) {
  bad <- checks[!checks$passed, , drop = FALSE]
  if (!nrow(bad)) return(NA_character_)
  paste0(bad$check[1L], ": ", bad$detail[1L])
}

# ---------------------------------------------------------------------------
# Comparisons
# ---------------------------------------------------------------------------

# Side-by-side table of one valid pair: one row per coefficient, then the
# theta row ("theta" for NB, "log(theta)" for ZINB, with SE(theta) or
# SE(log theta) in the standard-error columns), then the log-likelihood.
# `signed_diff` is ours - reference for the estimate, so a higher
# log-likelihood of ours shows as a positive number.
side_by_side <- function(ours, ref, model) {
  theta_label <- if (model == "NB") "theta" else "log(theta)"
  se_o <- sqrt(diag(ours$vcov))
  se_r <- sqrt(diag(ref$vcov))
  est_o <- c(ours$coef, ours$theta, ours$loglik)
  est_r <- c(ref$coef, ref$theta, ref$loglik)
  out <- data.frame(
    parameter = c(names(ours$coef), theta_label, "log-likelihood"),
    kind = c(rep("coefficient", length(ours$coef)), "theta",
             "log-likelihood"),
    estimate_ours = unname(est_o), estimate_reference = unname(est_r),
    estimate_abs_diff = unname(abs(est_o - est_r)),
    signed_diff = unname(est_o - est_r),
    se_ours = unname(c(se_o, ours$se_theta, NA)),
    se_reference = unname(c(se_r, ref$se_theta, NA)),
    stringsAsFactors = FALSE)
  out$se_abs_diff <- abs(out$se_ours - out$se_reference)
  out
}

# The four required quantities (handout Table 4) for one pair: maximum
# absolute difference, tolerance, and whether it is met. For an invalid pair
# the rows are "not assessed" and count as failed.
required_summary <- function(tab, model, pair_id, invalid_reason = NA_character_) {
  tol <- bench_tolerances()
  tol <- tol[tol$model == model, ]
  if (!is.na(invalid_reason)) {
    return(data.frame(
      pair = pair_id, model = model, quantity = tol$quantity,
      max_abs_diff = NA_real_, tolerance = tol$tolerance, passed = FALSE,
      note = paste("not assessed (invalid:", invalid_reason, ")"),
      stringsAsFactors = FALSE))
  }
  coef_rows <- tab$kind == "coefficient"
  d <- c(max(tab$estimate_abs_diff[coef_rows]),
         max(tab$se_abs_diff[coef_rows]),
         tab$estimate_abs_diff[tab$kind == "theta"],
         tab$estimate_abs_diff[tab$kind == "log-likelihood"])
  ll <- tab$signed_diff[tab$kind == "log-likelihood"]
  note <- c("", "", "",
            if (ll > 0) "log-likelihood: ours is higher"
            else if (ll < 0) "log-likelihood: reference is higher"
            else "log-likelihoods equal")
  data.frame(pair = pair_id, model = model, quantity = tol$quantity,
             max_abs_diff = d, tolerance = tol$tolerance,
             passed = is.finite(d) & d <= tol$tolerance, note = note,
             stringsAsFactors = FALSE)
}

# Relative difference |a - b| / max(|b|, 1e-12), elementwise.
rel_diff <- function(a, b) abs(a - b) / pmax(abs(b), 1e-12)

# Supplementary quantities (information only, no thresholds): SE(theta) for
# NB, SE(log theta) for ZINB, and for a numeric-Hessian ZINB pair the four
# Table 4 quantities as well. Absolute and relative differences.
supplementary_rows <- function(tab, pair_id, model, hessian) {
  th <- tab[tab$kind == "theta", ]
  out <- data.frame(
    pair = pair_id,
    quantity = if (model == "NB") "SE(theta)" else "SE(log theta)",
    max_abs_diff = th$se_abs_diff,
    max_rel_diff = rel_diff(th$se_ours, th$se_reference),
    stringsAsFactors = FALSE)
  if (identical(hessian, "numeric")) {
    cr <- tab$kind == "coefficient"
    ll <- tab$kind == "log-likelihood"
    extra <- data.frame(
      pair = pair_id,
      quantity = c("coefficients", "standard errors", "log(theta)",
                   "log-likelihood"),
      max_abs_diff = c(max(tab$estimate_abs_diff[cr]),
                       max(tab$se_abs_diff[cr]), th$estimate_abs_diff,
                       tab$estimate_abs_diff[ll]),
      max_rel_diff = c(max(rel_diff(tab$estimate_ours[cr],
                                    tab$estimate_reference[cr])),
                       max(rel_diff(tab$se_ours[cr], tab$se_reference[cr])),
                       rel_diff(th$estimate_ours, th$estimate_reference),
                       rel_diff(tab$estimate_ours[ll],
                                tab$estimate_reference[ll])),
      stringsAsFactors = FALSE)
    out <- rbind(extra, out)
  }
  out$note <- "supplementary (not a graded requirement)"
  out
}

# Our analytic against our numeric Hessian: the standard errors of every
# coefficient and of log(theta).
analytic_vs_numeric <- function(analytic, numeric, spec_id) {
  se_a <- c(sqrt(diag(analytic$vcov)), analytic$se_theta)
  se_n <- c(sqrt(diag(numeric$vcov)), numeric$se_theta)
  data.frame(pair = paste0(spec_id, "_analytic_vs_numeric"),
             quantity = "standard errors (ours: analytic against numeric Hessian)",
             max_abs_diff = max(abs(se_a - se_n)),
             max_rel_diff = max(rel_diff(se_a, se_n)),
             note = "supplementary (not a graded requirement)",
             stringsAsFactors = FALSE)
}

# Compare every pair. Returns the side-by-side tables, the validity tables,
# the required and the supplementary summaries.
compare_pairs <- function(pairs) {
  tables <- list()
  validity <- list()
  required <- list()
  supplementary <- list()
  std <- list()
  for (p in pairs) {
    o <- if (!is.null(p$ours$fit)) standardize_fit(p$ours$fit, p$model, "ours")
    r <- if (!is.null(p$ref$fit)) standardize_fit(p$ref$fit, p$model, "reference")
    std[[p$id]] <- list(ours = o, ref = r)
    chk <- validity_checks(o, r, p$ours$error, p$ref$error)
    validity[[p$id]] <- chk
    reason <- pair_failure_reason(chk)
    tab <- if (is.na(reason)) side_by_side(o, r, p$model)
    tables[[p$id]] <- tab
    graded <- p$model == "NB" || identical(p$hessian, "analytic")
    if (graded) {
      required[[p$id]] <- required_summary(tab, p$model, p$id, reason)
    }
    if (!is.null(tab)) {
      supplementary[[p$id]] <- supplementary_rows(tab, p$id, p$model, p$hessian)
    } else if (!graded) {
      # an invalid numeric-Hessian pair is reported, not graded
      supplementary[[p$id]] <- data.frame(
        pair = p$id, quantity = "all", max_abs_diff = NA_real_,
        max_rel_diff = NA_real_,
        note = paste("not assessed (invalid:", reason, ")"),
        stringsAsFactors = FALSE)
    }
  }
  # analytic against numeric Hessian, where both of our ZINB fits exist
  for (id in unique(vapply(pairs, function(p) p$spec_id, ""))) {
    a <- std[[paste0(id, "_analytic")]]$ours
    n <- std[[paste0(id, "_numeric")]]$ours
    if (!is.null(a) && !is.null(n)) {
      supplementary[[paste0(id, "_avn")]] <- analytic_vs_numeric(a, n, id)
    }
  }
  list(tables = tables, validity = validity,
       required = do.call(rbind, unname(required)),
       supplementary = do.call(rbind, unname(supplementary)))
}

# ---------------------------------------------------------------------------
# Timing
# ---------------------------------------------------------------------------

# Time the fits of one specification with bench::mark(): `iterations` runs
# of each fit, without memory profiling (it adds overhead) and keeping runs
# with garbage collections (dropping them would bias toward fast runs).
# Fits: ours (ZINB with each Hessian), the reference at its tight settings
# (the settings of the agreement tables), and the reference at its defaults
# (tight settings make it slower, so both are shown).
time_fits <- function(spec, iterations = 10L,
                      hessians = c("analytic", "numeric")) {
  f <- spec$formula
  d <- spec$data
  fits <- if (spec$model == "NB") {
    list(
      `nbreg() (default settings)` = function() nbzinbreg::nbreg(f, data = d),
      `glm.nb() (tight settings)` = function() {
        MASS::glm.nb(f, data = d, control = bench_reference_control("NB"))
      },
      `glm.nb() (default settings)` = function() MASS::glm.nb(f, data = d))
  } else {
    ours <- lapply(hessians, function(h) {
      force(h)
      function() nbzinbreg::zinbreg(
        f, data = d, control = nbzinbreg::zinbreg_control(hessian = h))
    })
    names(ours) <- sprintf("zinbreg() (default settings, %s Hessian)",
                           hessians)
    c(ours, list(
      `zeroinfl() (tight settings)` = function() {
        pscl::zeroinfl(f, data = d, dist = "negbin",
                       control = bench_reference_control("ZINB"))
      },
      `zeroinfl() (default settings)` = function() {
        pscl::zeroinfl(f, data = d, dist = "negbin")
      }))
  }
  rows <- lapply(names(fits), function(nm) {
    fn <- fits[[nm]]
    bm <- suppressWarnings(bench::mark(
      fn(), iterations = iterations, check = FALSE, memory = FALSE,
      filter_gc = FALSE, time_unit = "s"))
    t <- as.numeric(bm$time[[1L]])
    data.frame(spec = spec$id, data_name = spec$data_name, model = spec$model,
               fit = nm, median_s = stats::median(t), min_s = min(t),
               max_s = max(t), runs = length(t), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Provenance: which code produced the results
# ---------------------------------------------------------------------------

# R's own bookkeeping objects in a namespace, which are not package code.
fingerprint_ignored <- function(nms) {
  grepl("^\\.__", nms) | nms %in% c(".packageName")
}

# The text that is hashed: for each object, in sorted order of names, a
# header line with its name and then deparse() of the object with any
# source references removed (a function's arguments and body; a constant's
# value). Byte-compiled functions in an installed package deparse to the
# same text as the source, and deparse() never shows a function's
# environment, so the same code gives the same text whether it is read from
# R/*.R or from the installed package. (deparse() output may change between
# R versions, so fingerprints are only compared within one R version.)
code_fingerprint_text <- function(objects) {
  # radix sorting is the same in every locale, so the order (and the hash)
  # does not depend on the machine's collation settings
  nms <- sort(names(objects), method = "radix")
  unlist(lapply(nms, function(nm) {
    obj <- objects[[nm]]
    if (is.function(obj)) obj <- utils::removeSource(obj)
    c(paste0("### ", nm), deparse(obj))
  }), use.names = FALSE)
}

# MD5 of a character vector, via a temporary file (tools::md5sum() hashes
# files). Lines are written as bytes, so the encoding cannot change them.
md5_of_text <- function(text) {
  f <- tempfile(fileext = ".txt")
  on.exit(unlink(f))
  # the raw bytes of the lines, each ended by "\n"; writeBin() opens and
  # closes the file itself, so no connection can be left open
  writeBin(charToRaw(paste0(paste(text, collapse = "\n"), "\n")), f)
  unname(tools::md5sum(f))
}

# Per-object hashes and the overall hash of a named list of objects.
fingerprint_of <- function(objects) {
  per <- vapply(sort(names(objects), method = "radix"), function(nm) {
    md5_of_text(code_fingerprint_text(objects[nm]))
  }, character(1))
  list(md5 = md5_of_text(code_fingerprint_text(objects)), objects = per)
}

# Fingerprint of a source tree: evaluate every R/*.R file (in the collation
# order R uses without a Collate field: sorted names in the C locale) into a
# fresh environment, without keeping source references.
code_fingerprint_source <- function(path = ".") {
  files <- list.files(file.path(path, "R"), pattern = "\\.[RrSsq]$",
                      full.names = TRUE)
  files <- sort(files, method = "radix")   # C-locale order
  env <- new.env(parent = globalenv())
  for (f in files) sys.source(f, envir = env, keep.source = FALSE)
  nms <- ls(env, all.names = TRUE)
  nms <- nms[!fingerprint_ignored(nms)]
  fingerprint_of(mget(nms, envir = env))
}

# Fingerprint of the installed (or loaded) package's namespace.
code_fingerprint_installed <- function(package = "nbzinbreg") {
  ns <- asNamespace(package)
  nms <- ls(ns, all.names = TRUE)
  nms <- nms[!fingerprint_ignored(nms)]
  fingerprint_of(mget(nms, envir = ns))
}

# Compare two fingerprints: equal overall hash, and the objects that differ,
# are missing from the installed code, or are extra in it.
compare_fingerprints <- function(source, installed) {
  s <- source$objects
  i <- installed$objects
  common <- intersect(names(s), names(i))
  list(match = identical(source$md5, installed$md5),
       differ = common[s[common] != i[common]],
       missing_installed = setdiff(names(s), names(i)),
       extra_installed = setdiff(names(i), names(s)))
}

# Run git in `path`; NA if git or the repository is not available.
git_output <- function(path, args) {
  out <- tryCatch(suppressWarnings(system2(
    "git", c("-C", shQuote(path), args), stdout = TRUE, stderr = FALSE)),
    error = function(e) NA_character_)
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) {
    return(NA_character_)
  }
  out
}

# A short CPU description.
cpu_description <- function() {
  id <- Sys.getenv("PROCESSOR_IDENTIFIER")
  if (nzchar(id)) return(id)
  if (file.exists("/proc/cpuinfo")) {
    m <- grep("^model name", readLines("/proc/cpuinfo", warn = FALSE),
              value = TRUE)
    if (length(m)) return(trimws(sub(".*:", "", m[1L])))
  }
  out <- tryCatch(system2("sysctl", c("-n", "machdep.cpu.brand_string"),
                          stdout = TRUE, stderr = FALSE),
                  error = function(e) character())
  if (length(out)) out[1L] else NA_character_
}

# Where the results come from. `source_path` is the package source tree
# (the repository root). The provenance is "verified" only if
#   * the installed code's fingerprint equals the source tree's,
#   * the installed benchmark functions equal the source tree's copy,
#   * R/, inst/benchmark/benchmark_functions.R, DESCRIPTION and NAMESPACE
#     are unchanged from the git commit (so the fingerprint belongs to that
#     commit, and anyone can recompute it from it),
# and otherwise "unverified" with the reasons. Both fingerprints are
# computed here, in one R session, so the R version is the same for both.
provenance <- function(source_path = ".", package = "nbzinbreg") {
  reasons <- character()
  src <- tryCatch(code_fingerprint_source(source_path),
                  error = function(e) {
                    reasons <<- c(reasons, paste(
                      "the source tree could not be read:", conditionMessage(e)))
                    NULL
                  })
  inst <- code_fingerprint_installed(package)
  cmp <- if (!is.null(src)) compare_fingerprints(src, inst)
  if (!is.null(cmp) && !cmp$match) {
    reasons <- c(reasons, paste0(
      "the installed code differs from the source tree",
      if (length(cmp$differ)) paste0(" (differs: ",
                                     paste(cmp$differ, collapse = ", "), ")"),
      if (length(cmp$missing_installed)) paste0(
        " (missing from the installed package: ",
        paste(cmp$missing_installed, collapse = ", "), ")"),
      if (length(cmp$extra_installed)) paste0(
        " (only in the installed package: ",
        paste(cmp$extra_installed, collapse = ", "), ")")))
  }
  bf_src <- file.path(source_path, "inst", "benchmark", "benchmark_functions.R")
  bf_inst <- system.file("benchmark", "benchmark_functions.R", package = package)
  bf_match <- file.exists(bf_src) && nzchar(bf_inst) &&
    identical(unname(tools::md5sum(bf_src)), unname(tools::md5sum(bf_inst)))
  if (!bf_match) {
    reasons <- c(reasons, paste(
      "the installed benchmark functions differ from the source tree's",
      "inst/benchmark/benchmark_functions.R"))
  }
  commit <- git_output(source_path, c("rev-parse", "HEAD"))
  dirty <- git_output(source_path, c(
    "status", "--porcelain", "--", "R", "inst/benchmark/benchmark_functions.R",
    "DESCRIPTION", "NAMESPACE"))
  if (anyNA(commit)) {
    reasons <- c(reasons, "the git commit of the source tree is unknown")
  } else if (anyNA(dirty) || length(dirty)) {
    reasons <- c(reasons, paste(
      "R/, the benchmark functions, DESCRIPTION or NAMESPACE have",
      "uncommitted changes"))
  }
  desc <- utils::packageDescription(package)
  si <- Sys.info()
  versions <- vapply(c(package, "MASS", "pscl", "AER", "bench"), function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p))
    else NA_character_
  }, character(1))
  list(
    status = if (length(reasons)) "unverified" else "verified",
    reasons = reasons,
    fingerprint_installed = inst$md5,
    fingerprint_source = if (!is.null(src)) src$md5 else NA_character_,
    fingerprint_objects = inst$objects,
    benchmark_functions_md5 = if (nzchar(bf_inst))
      unname(tools::md5sum(bf_inst)) else NA_character_,
    git_commit = if (anyNA(commit)) NA_character_ else commit[1L],
    git_clean = !anyNA(dirty) && !length(dirty),
    package_version = desc$Version,
    package_packaged = if (is.null(desc$Packaged)) NA_character_ else desc$Packaged,
    package_built = if (is.null(desc$Built)) NA_character_ else desc$Built,
    package_remote_sha = if (is.null(desc$RemoteSha)) NA_character_ else desc$RemoteSha,
    r_version = R.version.string,
    platform = R.version$platform,
    os = paste(si[["sysname"]], si[["release"]]),
    cpu = cpu_description(),
    cores = parallel::detectCores(),
    versions = versions,
    date = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"))
}

# Check saved results against the code installed now (used by the report and
# the supplementary article before showing the saved timings). Returns
# list(status, message): "matches" (same code, same R version), "differs"
# (the installed code is not the code that produced the timings), or
# "cannot verify" (another R version, so deparse() may differ).
check_saved_provenance <- function(saved, package = "nbzinbreg") {
  now <- code_fingerprint_installed(package)
  if (!identical(saved$r_version, R.version.string)) {
    return(list(status = "cannot verify", message = sprintf(paste(
      "The timings were produced under %s, and this document is rendered",
      "under %s; the code fingerprints are only comparable within one R",
      "version, so it cannot be verified here that the installed code is",
      "the code that produced the timings."), saved$r_version,
      R.version.string)))
  }
  if (identical(now$md5, saved$fingerprint_installed)) {
    return(list(status = "matches", message = sprintf(paste(
      "The installed package's code fingerprint (%s) matches the code that",
      "produced the timings (commit %s)."), now$md5, saved$git_commit)))
  }
  common <- intersect(names(now$objects), names(saved$fingerprint_objects))
  differ <- common[now$objects[common] != saved$fingerprint_objects[common]]
  list(status = "differs", message = sprintf(paste(
    "WARNING: the installed package's code differs from the code that",
    "produced the timings (commit %s)%s. The timings below were measured on",
    "different code; re-run inst/benchmark/benchmark.R."), saved$git_commit,
    if (length(differ)) paste0("; differing objects: ",
                               paste(differ, collapse = ", ")) else ""))
}

# ---------------------------------------------------------------------------
# Status, running and saving
# ---------------------------------------------------------------------------

# Overall status: "passed" only if every pair is valid, every required
# quantity meets its tolerance, and the provenance is verified. The failures
# table lists every reason. Estimates are never adjusted to make a check
# pass: the saved tables hold the actual numbers either way.
benchmark_status <- function(comparison, prov) {
  fails <- list()
  for (id in names(comparison$validity)) {
    chk <- comparison$validity[[id]]
    bad <- chk[!chk$passed, , drop = FALSE]
    if (nrow(bad)) {
      fails[[length(fails) + 1L]] <- data.frame(
        pair = id, item = bad$check, value = NA_real_, tolerance = NA_real_,
        reason = bad$detail, stringsAsFactors = FALSE)
    }
  }
  req <- comparison$required
  bad <- req[!req$passed, , drop = FALSE]
  if (nrow(bad)) {
    fails[[length(fails) + 1L]] <- data.frame(
      pair = bad$pair, item = bad$quantity, value = bad$max_abs_diff,
      tolerance = bad$tolerance,
      reason = ifelse(is.na(bad$max_abs_diff), bad$note,
                      "maximum absolute difference above the tolerance"),
      stringsAsFactors = FALSE)
  }
  if (prov$status != "verified") {
    fails[[length(fails) + 1L]] <- data.frame(
      pair = "(provenance)", item = "code provenance", value = NA_real_,
      tolerance = NA_real_, reason = prov$reasons, stringsAsFactors = FALSE)
  }
  failures <- if (length(fails)) do.call(rbind, fails) else data.frame(
    pair = character(), item = character(), value = numeric(),
    tolerance = numeric(), reason = character(), stringsAsFactors = FALSE)
  list(status = if (nrow(failures)) "failed" else "passed",
       failures = failures)
}

# Fit, compare and (optionally) time everything. `iterations = 0` skips the
# timing (the report computes only the agreement live).
run_benchmark <- function(source_path = ".", iterations = 10L,
                          hessians = c("analytic", "numeric"),
                          progress = TRUE) {
  say <- function(...) if (progress) message(...)
  say("Provenance ...")
  prov <- provenance(source_path)
  specs <- bench_specs()
  pairs <- list()
  for (s in specs) {
    hs <- if (s$model == "NB") "analytic" else hessians
    for (h in hs) {
      say("Fitting ", s$id, if (s$model == "ZINB") paste0(" (", h, " Hessian)"),
          " ...")
      pairs[[length(pairs) + 1L]] <- fit_pair(s, h)
    }
  }
  comparison <- compare_pairs(pairs)
  timing <- NULL
  if (iterations > 0L) {
    valid_specs <- Filter(function(s) {
      ids <- names(comparison$validity)[startsWith(names(comparison$validity),
                                                   s$id)]
      all(vapply(comparison$validity[ids], pair_valid, logical(1)))
    }, specs)
    skipped <- setdiff(vapply(specs, `[[`, "", "id"),
                       vapply(valid_specs, `[[`, "", "id"))
    timing <- do.call(rbind, lapply(valid_specs, function(s) {
      say("Timing ", s$id, " (", iterations, " runs per fit) ...")
      time_fits(s, iterations, if (s$model == "ZINB") hessians else "analytic")
    }))
    attr(timing, "skipped") <- skipped
  }
  st <- benchmark_status(comparison, prov)
  list(status = st$status, failures = st$failures, provenance = prov,
       tables = comparison$tables, validity = comparison$validity,
       required = comparison$required,
       supplementary = comparison$supplementary, timing = timing,
       iterations = iterations,
       warnings = lapply(pairs, function(p) {
         list(ours = p$ours$warnings, reference = p$ref$warnings)
       }))
}

# Save the results: the full object as .rds, readable CSV copies, and a
# status file whose first line is PASSED or FAILED. Always called, also for
# a failed run, so the actual comparisons are kept for diagnosis.
save_benchmark <- function(results, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(results, file.path(out_dir, "benchmark_results.rds"))
  utils::write.csv(results$required,
                   file.path(out_dir, "benchmark_required.csv"),
                   row.names = FALSE)
  utils::write.csv(results$supplementary,
                   file.path(out_dir, "benchmark_supplementary.csv"),
                   row.names = FALSE)
  if (!is.null(results$timing)) {
    utils::write.csv(results$timing, file.path(out_dir, "benchmark_timing.csv"),
                     row.names = FALSE)
  }
  p <- results$provenance
  lines <- c(
    toupper(results$status),
    "",
    if (nrow(results$failures)) {
      c("Failures:", sprintf("- %s / %s: %s%s", results$failures$pair,
                             results$failures$item, results$failures$reason,
                             ifelse(is.na(results$failures$value), "",
                                    sprintf(" (%.3g > %.3g)",
                                            results$failures$value,
                                            results$failures$tolerance))), "")
    } else c("All validity checks passed and every required quantity is within its tolerance.", ""),
    "Provenance:",
    sprintf("- status: %s%s", p$status,
            if (length(p$reasons)) paste0(" (", paste(p$reasons, collapse = "; "), ")") else ""),
    sprintf("- git commit: %s (relevant paths clean: %s)", p$git_commit, p$git_clean),
    sprintf("- code fingerprint (installed / source): %s / %s",
            p$fingerprint_installed, p$fingerprint_source),
    sprintf("- nbzinbreg %s, built %s", p$package_version, p$package_built),
    sprintf("- %s on %s (%s; %s, %s cores)", p$r_version, p$platform, p$os,
            p$cpu, p$cores),
    sprintf("- packages: %s", paste(names(p$versions), p$versions,
                                    collapse = ", ")),
    sprintf("- run at %s, %s timing runs per fit", p$date, results$iterations))
  writeLines(lines, file.path(out_dir, "benchmark_status.txt"), useBytes = TRUE)
  invisible(out_dir)
}

# Functions for the simulation study (handout Section 5.1), shared by the
# script (inst/benchmark/simulation.R), the report and the tests
# (tests/testthat/test-simulation.R).
#
# What the study does: it simulates R data sets of size n from a
# zero-inflated negative binomial (ZINB) model with known parameters, fits
# each with zinbreg(), and checks how well the known values are recovered:
#   * bias: the average estimate minus the true value (0 for an unbiased
#     estimator);
#   * Monte Carlo SD: the actual spread of the estimates over the data sets;
#   * mean estimated SE: the average of the standard errors zinbreg()
#     reports. If those are right, mean SE is about equal to the Monte Carlo
#     SD (their ratio is about 1);
#   * coverage: the share of 95% Wald intervals (estimate +- 1.96 SE) that
#     contain the true value; about 95% if the intervals are right.
# Every summary comes with a Monte Carlo standard error (MC SE): the
# uncertainty that remains because only R data sets were simulated.
#
# The data-generating model (fixed by the owner's plan):
#   covariates, drawn once (a fixed design): x1 ~ N(0, 1), x2 ~ Bernoulli(0.5)
#   and z ~ N(0, 1);
#   count part: log mu = 1 + 0.5 x1 - 0.5 x2;
#   zero part:  logit pi = -0.8 + 0.6 x1 - 0.5 z;
#   theta = 1.5.
# Each observation is a structural zero with probability pi, and otherwise
# a negative binomial count with mean mu and dispersion theta. The fitted
# model, zinbreg(y ~ x1 + x2 | x1 + z), is the model that generated the
# data.
#
# rbinom() and rnbinom() are used here to draw the simulated data; this file
# is not part of the package code in R/.

bench_file <- system.file("benchmark", "benchmark_functions.R",
                          package = "nbzinbreg")
if (!nzchar(bench_file)) {
  stop("inst/benchmark/benchmark_functions.R was not found in the installed ",
       "nbzinbreg package.")
}
# provenance(), check_saved_provenance(), covariance_valid(), git helpers
sys.source(bench_file, envir = environment(), keep.source = FALSE)

# ---------------------------------------------------------------------------
# The model and the data
# ---------------------------------------------------------------------------

# The true parameter values, named as coef() of a zinbreg fit names them,
# plus log(theta) and theta.
sim_truth <- function() {
  beta <- c(`(Intercept)` = 1, x1 = 0.5, x2 = -0.5)
  gamma <- c(`(Intercept)` = -0.8, x1 = 0.6, z = -0.5)
  theta <- 1.5
  list(beta = beta, gamma = gamma, theta = theta,
       params = c(stats::setNames(beta, paste0("count_", names(beta))),
                  stats::setNames(gamma, paste0("zero_", names(gamma))),
                  `log(theta)` = log(theta), theta = theta))
}

sim_param_names <- function() names(sim_truth()$params)

# Run `expr` with the random-number seed `seed`, and restore the caller's
# random-number state afterwards (so a function call never changes the
# stream of the code around it).
with_rng_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv())
    else rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed)
  expr
}

# The fixed design: the covariates, drawn once with `seed`.
sim_design <- function(n = 500L, seed = 651L) {
  with_rng_seed(seed, {
    x1 <- stats::rnorm(n)
    x2 <- stats::rbinom(n, 1L, 0.5)
    z <- stats::rnorm(n)
    data.frame(x1 = x1, x2 = x2, z = z)
  })
}

# The means mu, the structural-zero probabilities pi and the NB probability
# of a zero f0 for every row of the design.
sim_linear <- function(design, truth = sim_truth()) {
  b <- truth$beta
  g <- truth$gamma
  mu <- exp(b[[1]] + b[["x1"]] * design$x1 + b[["x2"]] * design$x2)
  pi <- stats::plogis(g[[1]] + g[["x1"]] * design$x1 + g[["z"]] * design$z)
  f0 <- (truth$theta / (truth$theta + mu))^truth$theta
  list(mu = mu, pi = pi, f0 = f0)
}

# One simulated response: a structural zero with probability pi, otherwise a
# negative binomial count. Both vectors are always drawn in full, so the
# random-number stream does not depend on the outcomes.
sim_response <- function(design, truth = sim_truth(), seed) {
  p <- sim_linear(design, truth)
  with_rng_seed(seed, {
    structural <- stats::rbinom(nrow(design), 1L, p$pi)
    counts <- stats::rnbinom(nrow(design), size = truth$theta, mu = p$mu)
    list(y = ifelse(structural == 1L, 0, counts), structural = structural)
  })
}

# The design's exact expected values (averaged over its rows): the share of
# zeros, P(Y = 0) = pi + (1 - pi) f0; the share of structural zeros among
# the zeros; and the mean of Y, (1 - pi) mu.
design_moments <- function(design, truth = sim_truth()) {
  p <- sim_linear(design, truth)
  p0 <- p$pi + (1 - p$pi) * p$f0
  c(zero_share = mean(p0),
    structural_share_of_zeros = mean(p$pi) / mean(p0),
    structural_share = mean(p$pi),
    mean_y = mean((1 - p$pi) * p$mu))
}

# ---------------------------------------------------------------------------
# One replication: the record of a fit
# ---------------------------------------------------------------------------

# 95% Wald intervals use the normal quantile 1.959964.
wald_z <- function() stats::qnorm(0.975)

# The record of one fitted replication (or of a fit that failed with
# `error`). Every parameter (the six coefficients, log(theta) and theta) gets
# its estimate, SE, 95% Wald limits, whether the interval covers the truth,
# and a validity flag (estimate, SE and limits finite, SE > 0, lower <
# upper). theta's interval is the theta-scale one, theta +- 1.96 SE(theta);
# the log-scale interval for theta is exp() of log(theta)'s interval, so its
# coverage is that of log(theta).
#
# The replication class (the first that applies, with a reason):
#   error              zinbreg() failed (the message is kept)
#   not converged      converged = FALSE or the gradient check did not pass
#   boundary           a diverging coefficient, a collapsed zero part (every
#                      pi < 1e-8), or theta on a bound
#   invalid covariance vcov_full not finite, not symmetric or not positive
#                      definite
#   invalid values     an estimate, SE or interval limit of a parameter is
#                      not valid
#   usable             otherwise
replication_record <- function(r, fit, error = NA_character_,
                               warnings = character(), seconds = NA_real_,
                               truth = sim_truth()) {
  pn <- sim_param_names()
  na <- stats::setNames(rep(NA_real_, length(pn)), pn)
  rec <- list(replication = r, class = "error", reason = error,
              estimate = na, se = na, lower = na, upper = na,
              covers = stats::setNames(rep(NA, length(pn)), pn),
              valid = stats::setNames(rep(FALSE, length(pn)), pn),
              converged = NA, gradient_ok = NA, newton_decrement = NA_real_,
              iterations = NA_integer_, diverging = NA_character_,
              theta_bound = NA_character_, zero_collapsed = NA,
              seconds = seconds, n_warnings = length(warnings),
              first_warning = if (length(warnings)) warnings[1L]
                              else NA_character_)
  if (is.null(fit)) return(rec)

  cf <- c(stats::setNames(fit$coefficients$count,
                          paste0("count_", names(fit$coefficients$count))),
          stats::setNames(fit$coefficients$zero,
                          paste0("zero_", names(fit$coefficients$zero))))
  se_coef <- sqrt(diag(fit$vcov))[names(cf)]
  est <- c(cf, `log(theta)` = log(fit$theta), theta = fit$theta)
  se <- c(se_coef, `log(theta)` = fit$SE.logtheta, theta = fit$SE.theta)
  est <- est[pn]
  se <- se[pn]
  z <- wald_z()
  lower <- est - z * se
  upper <- est + z * se
  tv <- truth$params[pn]
  rec$estimate <- est
  rec$se <- se
  rec$lower <- lower
  rec$upper <- upper
  rec$covers <- lower <= tv & tv <= upper
  rec$valid <- is.finite(est) & is.finite(se) & se > 0 & is.finite(lower) &
    is.finite(upper) & lower < upper
  rec$valid[is.na(rec$valid)] <- FALSE

  rec$converged <- isTRUE(fit$converged)
  rec$gradient_ok <- isTRUE(fit$gradient_check$ok)
  rec$newton_decrement <- if (is.null(fit$gradient_check$decrement)) NA_real_
                          else fit$gradient_check$decrement
  rec$iterations <- if (is.null(fit$iter)) NA_integer_ else as.integer(fit$iter)
  rec$diverging <- paste(fit$diverging, collapse = ", ")
  rec$theta_bound <- if (is.null(fit$theta_bound)) NA_character_
                     else fit$theta_bound
  rec$zero_collapsed <- !is.null(fit$pi) && all(fit$pi < 1e-8)

  cv <- covariance_valid(fit$vcov_full)
  bad_values <- pn[!rec$valid]
  cls <- if (!rec$converged || !rec$gradient_ok) {
    c("not converged", sprintf(
      "converged = %s, gradient check passed = %s", rec$converged,
      rec$gradient_ok))
  } else if (length(fit$diverging) || rec$zero_collapsed ||
             !identical(rec$theta_bound, "none")) {
    c("boundary", paste(c(
      if (length(fit$diverging)) paste("diverging:",
                                       paste(fit$diverging, collapse = ", ")),
      if (rec$zero_collapsed) "zero part collapsed (every pi < 1e-8)",
      if (!identical(rec$theta_bound, "none")) paste("theta on its",
                                                     rec$theta_bound, "bound")),
      collapse = "; "))
  } else if (!cv$ok) {
    c("invalid covariance", cv$detail)
  } else if (length(bad_values)) {
    c("invalid values", paste("not valid:", paste(bad_values, collapse = ", ")))
  } else {
    c("usable", NA_character_)
  }
  rec$class <- cls[1L]
  rec$reason <- cls[2L]
  rec
}

# The profile-likelihood interval for theta of one replication, classified:
#   finite          both endpoints finite: covers if lower <= truth <= upper
#   unbounded       upper endpoint Inf with status "unbounded" (certified by
#                   confint()) and a finite lower endpoint: a valid interval
#                   [lower, Inf), which covers if lower <= truth
#   unresolved      an endpoint NA (beyond theta_max, below theta_min, inner
#                   fit failed): cannot be evaluated
#   error           confint() failed: cannot be evaluated
#   fit not usable  the fit itself was not usable: not attempted
# `ci` is the result of confint(fit, "theta") (or NULL), `error` its message.
profile_record <- function(ci, error = NA_character_, fit_usable = TRUE,
                           truth = sim_truth(), warnings = character()) {
  out <- list(class = "fit not usable", lower = NA_real_, upper = NA_real_,
              status_lower = NA_character_, status_upper = NA_character_,
              covers = NA, detail = NA_character_,
              n_warnings = length(warnings))
  if (!fit_usable) return(out)
  if (is.null(ci)) {
    out$class <- "error"
    out$detail <- error
    return(out)
  }
  lim <- unclass(ci)["theta", ]
  st <- attr(ci, "status")["theta", ]
  out$lower <- unname(lim[1L])
  out$upper <- unname(lim[2L])
  out$status_lower <- unname(st[1L])
  out$status_upper <- unname(st[2L])
  th <- truth$theta
  if (identical(unname(st), c("finite", "finite")) &&
      all(is.finite(lim))) {
    out$class <- "finite"
    out$covers <- out$lower <= th && th <= out$upper
  } else if (identical(unname(st), c("finite", "unbounded")) &&
             is.finite(out$lower) && identical(out$upper, Inf)) {
    out$class <- "unbounded"
    out$covers <- out$lower <= th
  } else {
    out$class <- "unresolved"
    out$detail <- paste("endpoint status:", paste(st, collapse = " / "))
  }
  out
}

# Fit one replication: zinbreg() at its default settings, with warnings and
# errors captured, then the profile interval if `profile` is TRUE.
fit_replication <- function(r, design, y, truth = sim_truth(),
                            profile = FALSE) {
  d <- data.frame(y = y, design)
  warns <- character()
  t0 <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
    tryCatch(nbzinbreg::zinbreg(y ~ x1 + x2 | x1 + z, data = d),
             error = function(e) e),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  seconds <- proc.time()[["elapsed"]] - t0
  err <- NA_character_
  if (inherits(fit, "error")) {
    err <- conditionMessage(fit)
    fit <- NULL
  }
  rec <- replication_record(r, fit, err, warns, seconds, truth)
  rec$zeros <- sum(y == 0)
  rec$profile <- NULL
  if (profile) {
    pw <- character()
    pt0 <- proc.time()[["elapsed"]]
    ci <- if (rec$class == "usable") {
      withCallingHandlers(
        tryCatch(stats::confint(fit, "theta"), error = function(e) e),
        warning = function(w) {
          pw <<- c(pw, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
    }
    perr <- NA_character_
    if (inherits(ci, "error")) {
      perr <- conditionMessage(ci)
      ci <- NULL
    }
    rec$profile <- profile_record(ci, perr, rec$class == "usable", truth, pw)
    rec$profile$seconds <- proc.time()[["elapsed"]] - pt0
  }
  rec
}

# ---------------------------------------------------------------------------
# Records as a table
# ---------------------------------------------------------------------------

# One row per replication: diagnostics, then per parameter p the columns
# est_p, se_p, lower_p, upper_p, covers_p and valid_p, then the profile
# interval (NA where not computed).
records_table <- function(records) {
  pn <- sim_param_names()
  rows <- lapply(records, function(rec) {
    base <- data.frame(
      replication = rec$replication, class = rec$class,
      reason = rec$reason, zeros = rec$zeros, converged = rec$converged,
      gradient_ok = rec$gradient_ok, newton_decrement = rec$newton_decrement,
      iterations = rec$iterations, diverging = rec$diverging,
      theta_bound = rec$theta_bound, zero_collapsed = rec$zero_collapsed,
      seconds = rec$seconds, n_warnings = rec$n_warnings,
      first_warning = rec$first_warning, stringsAsFactors = FALSE)
    for (what in c("estimate", "se", "lower", "upper", "covers", "valid")) {
      vals <- rec[[what]][pn]
      names(vals) <- paste0(what, "_", pn)
      base <- cbind(base, as.data.frame(as.list(vals), check.names = FALSE))
    }
    p <- rec$profile
    base$profile_class <- if (is.null(p)) NA_character_ else p$class
    base$profile_lower <- if (is.null(p)) NA_real_ else p$lower
    base$profile_upper <- if (is.null(p)) NA_real_ else p$upper
    base$profile_status <- if (is.null(p)) NA_character_ else
      paste(p$status_lower, p$status_upper, sep = " / ")
    base$profile_covers <- if (is.null(p)) NA else p$covers
    base$profile_detail <- if (is.null(p)) NA_character_ else p$detail
    base$profile_seconds <- if (is.null(p)) NA_real_ else p$seconds
    base
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Summaries with Monte Carlo standard errors
# ---------------------------------------------------------------------------

# The summaries of one parameter from its usable estimates `e`, standard
# errors `s` and coverage flags `cv`, with true value `t`. RMSE is computed
# directly, sqrt(mean((e - t)^2)).
param_summaries <- function(e, s, cv, t) {
  sd_e <- stats::sd(e)
  c(mean = mean(e), bias = mean(e) - t,
    rel_bias = if (t != 0) (mean(e) - t) / t else NA_real_,
    mc_sd = sd_e, mean_se = mean(s), se_ratio = mean(s) / sd_e,
    rmse = sqrt(mean((e - t)^2)), coverage = mean(cv))
}

# Summarise the records table. Over the R_u usable replications, per
# parameter:
#   analytic MC SEs: mean and bias SD/sqrt(R_u); mean SE sd(SE)/sqrt(R_u);
#     coverage sqrt(c (1 - c) / R_u) (binomial);
#   bootstrap MC SEs (B resamples of the usable rows, the same rows for every
#     parameter, seed `seed`): relative bias, MC SD, SE ratio and RMSE;
#   the normal-theory MC SE of the SD, SD / sqrt(2 (R_u - 1)), reported for
#     comparison (it assumes normally distributed estimates);
#   coverage bounds over all R replications: every non-usable replication
#     counted as not covering (lower) or covering (upper).
summarise_simulation <- function(tab, truth = sim_truth(), B = 2000L,
                                 seed = 651L) {
  pn <- sim_param_names()
  R <- nrow(tab)
  u <- tab$class == "usable"
  Ru <- sum(u)
  nu <- R - Ru
  E <- as.matrix(tab[u, paste0("estimate_", pn), drop = FALSE])
  S <- as.matrix(tab[u, paste0("se_", pn), drop = FALSE])
  CV <- as.matrix(tab[u, paste0("covers_", pn), drop = FALSE])
  colnames(E) <- colnames(S) <- colnames(CV) <- pn
  tv <- truth$params[pn]

  point <- vapply(pn, function(p) param_summaries(E[, p], S[, p], CV[, p],
                                                  tv[[p]]),
                  numeric(8))
  boot <- if (Ru >= 2L) {
    with_rng_seed(seed, {
      sims <- replicate(B, {
        i <- sample.int(Ru, Ru, replace = TRUE)
        vapply(pn, function(p) {
          param_summaries(E[i, p], S[i, p], CV[i, p], tv[[p]])[
            c("rel_bias", "mc_sd", "se_ratio", "rmse")]
        }, numeric(4))
      }, simplify = "array")
      apply(sims, c(1, 2), stats::sd)
    })
  }
  covered <- colSums(CV)
  conv_ok <- !(tab$class %in% c("error", "not converged"))
  valid_counts <- vapply(pn, function(p) sum(tab[conv_ok, paste0("valid_", p)]),
                         numeric(1))
  out <- data.frame(
    parameter = pn, truth = unname(tv), n_usable = Ru,
    mean = point["mean", ], mean_mcse = point["mc_sd", ] / sqrt(Ru),
    bias = point["bias", ], bias_mcse = point["mc_sd", ] / sqrt(Ru),
    rel_bias = point["rel_bias", ],
    rel_bias_mcse = if (!is.null(boot)) boot["rel_bias", ] else NA_real_,
    mc_sd = point["mc_sd", ],
    mc_sd_mcse = if (!is.null(boot)) boot["mc_sd", ] else NA_real_,
    mc_sd_mcse_normal = point["mc_sd", ] / sqrt(2 * (Ru - 1)),
    mean_se = point["mean_se", ],
    mean_se_mcse = apply(S, 2, stats::sd) / sqrt(Ru),
    se_ratio = point["se_ratio", ],
    se_ratio_mcse = if (!is.null(boot)) boot["se_ratio", ] else NA_real_,
    rmse = point["rmse", ],
    rmse_mcse = if (!is.null(boot)) boot["rmse", ] else NA_real_,
    coverage = point["coverage", ],
    coverage_mcse = sqrt(point["coverage", ] * (1 - point["coverage", ]) / Ru),
    coverage_lower_bound = covered / R,
    coverage_upper_bound = (covered + nu) / R,
    n_valid_converged = valid_counts,
    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

# The replication classes: counts, and the replication numbers and reasons
# of every non-usable one.
class_summary <- function(tab) {
  levels <- c("usable", "boundary", "invalid covariance", "invalid values",
              "not converged", "error")
  counts <- table(factor(tab$class, levels = levels))
  bad <- tab[tab$class != "usable", c("replication", "class", "reason"),
             drop = FALSE]
  list(counts = counts, non_usable = bad)
}

# The profile-interval classes over the replications that had one, coverage
# over the evaluable (finite or unbounded) intervals with its binomial MC
# SE, and sensitivity bounds over all of them: every non-evaluable interval
# counted as not covering (lower) or covering (upper).
profile_summary <- function(tab) {
  p <- tab[!is.na(tab$profile_class), , drop = FALSE]
  levels <- c("finite", "unbounded", "unresolved", "error", "fit not usable")
  counts <- table(factor(p$profile_class, levels = levels))
  n_all <- nrow(p)
  ev <- p$profile_class %in% c("finite", "unbounded")
  P <- sum(ev)
  covered <- sum(p$profile_covers[ev])
  c_hat <- if (P > 0L) covered / P else NA_real_
  list(counts = counts, n = n_all, evaluable = P,
       coverage = c_hat,
       coverage_mcse = if (P > 0L) sqrt(c_hat * (1 - c_hat) / P) else NA_real_,
       lower_bound = if (n_all > 0L) covered / n_all else NA_real_,
       upper_bound = if (n_all > 0L) (covered + n_all - P) / n_all else NA_real_,
       median_seconds = if (n_all > 0L) stats::median(p$profile_seconds,
                                                      na.rm = TRUE)
                        else NA_real_)
}

# ---------------------------------------------------------------------------
# Running, status and saving
# ---------------------------------------------------------------------------

# Run the whole study. `source_path` is the package source tree (for the
# provenance).
run_simulation <- function(R = 1000L, n = 500L, seed = 651L,
                           profile_R = 200L, B = 2000L, source_path = ".",
                           progress = TRUE) {
  say <- function(...) if (progress) message(...)
  prov <- provenance(source_path,
                     files = c("benchmark_functions.R", "simulation_functions.R"))
  truth <- sim_truth()
  design <- sim_design(n, seed)
  moments <- design_moments(design, truth)
  records <- vector("list", R)
  t0 <- proc.time()[["elapsed"]]
  for (r in seq_len(R)) {
    sim <- sim_response(design, truth, seed + r)
    records[[r]] <- fit_replication(r, design, sim$y, truth,
                                    profile = r <= profile_R)
    records[[r]]$structural <- sum(sim$structural)
    if (r %% 50L == 0L || r == R) {
      say(sprintf("replication %d of %d (%.0f s)", r, R,
                  proc.time()[["elapsed"]] - t0))
    }
  }
  tab <- records_table(records)
  tab$structural <- vapply(records, function(x) as.numeric(x$structural), 1)
  summary <- summarise_simulation(tab, truth, B, seed)
  list(
    settings = list(R = R, n = n, seed = seed, profile_R = profile_R, B = B,
                    formula = "y ~ x1 + x2 | x1 + z",
                    rng_kind = RNGkind(), r_version = R.version.string),
    truth = truth, design = design, design_moments = moments,
    observed = c(zero_share = sum(tab$zeros) / (R * n),
                 structural_share = sum(tab$structural) / (R * n)),
    replicates = tab, summary = summary,
    classes = class_summary(tab), profile = profile_summary(tab),
    provenance = prov,
    elapsed_seconds = proc.time()[["elapsed"]] - t0)
}

# The run is "complete" when its provenance is verified and every
# replication was attempted. Fits that fail are results of the study (they
# are counted and reported), not a failure of the run.
simulation_status <- function(results) {
  reasons <- c(
    if (results$provenance$status != "verified") {
      paste("provenance unverified:",
            paste(results$provenance$reasons, collapse = "; "))
    },
    if (nrow(results$replicates) != results$settings$R) {
      sprintf("only %d of %d replications were recorded",
              nrow(results$replicates), results$settings$R)
    })
  list(status = if (length(reasons)) "failed" else "complete",
       reasons = reasons)
}

# Save the results: the full object as .rds, the summary and the
# per-replication table as CSV, and a status file whose first line is
# COMPLETE or FAILED. Always called, also for a failed run.
save_simulation <- function(results, out_dir) {
  st <- simulation_status(results)
  results$status <- st$status
  results$status_reasons <- st$reasons
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(results, file.path(out_dir, "simulation_results.rds"))
  utils::write.csv(results$summary, file.path(out_dir, "simulation_summary.csv"),
                   row.names = FALSE)
  utils::write.csv(results$replicates,
                   file.path(out_dir, "simulation_replicates.csv"),
                   row.names = FALSE)
  s <- results$settings
  cl <- results$classes$counts
  pr <- results$profile
  p <- results$provenance
  lines <- c(
    toupper(st$status),
    if (length(st$reasons)) paste("-", st$reasons),
    "",
    sprintf("R = %d data sets of n = %d, seed %d, model %s; profile intervals on the first %d; %d bootstrap resamples for the MC SEs.",
            s$R, s$n, s$seed, s$formula, s$profile_R, s$B),
    sprintf("Replication classes: %s.",
            paste(names(cl), as.integer(cl), sep = " ", collapse = ", ")),
    sprintf("Profile intervals: %s; coverage %.3f over %d evaluable.",
            paste(names(pr$counts), as.integer(pr$counts), sep = " ",
                  collapse = ", "), pr$coverage, pr$evaluable),
    sprintf("Zero share: expected %.4f, observed %.4f; structural share: expected %.4f, observed %.4f.",
            results$design_moments[["zero_share"]],
            results$observed[["zero_share"]],
            results$design_moments[["structural_share"]],
            results$observed[["structural_share"]]),
    sprintf("RNG: %s; elapsed %.0f s.", paste(s$rng_kind, collapse = ", "),
            results$elapsed_seconds),
    "",
    "Provenance:",
    sprintf("- status: %s%s", p$status,
            if (length(p$reasons)) paste0(" (", paste(p$reasons, collapse = "; "), ")") else ""),
    sprintf("- git commit: %s (relevant paths clean: %s)", p$git_commit, p$git_clean),
    sprintf("- code fingerprint (installed / source): %s / %s",
            p$fingerprint_installed, p$fingerprint_source),
    sprintf("- %s on %s (%s; %s, %s cores)", p$r_version, p$platform, p$os,
            p$cpu, p$cores),
    sprintf("- run at %s", p$date))
  writeLines(lines, file.path(out_dir, "simulation_status.txt"),
             useBytes = TRUE)
  invisible(results)
}

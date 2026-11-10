# Functions for the NMES1988 data analysis (handout Section 5.2), shared by
# the bootstrap script (inst/benchmark/nmes_bootstrap.R), the report and the
# tests (tests/testthat/test-nmes.R).
#
# The data: AER::NMES1988, 4406 people aged 66 and over from the US National
# Medical Expenditure Survey of 1987-88 (Deb and Trivedi 1997). The response
# is `visits`, the number of physician office visits.
#
# The models (fixed by the owner's plan):
#   * Poisson regression, fitted with glm() as the reference;
#   * negative binomial (NB) regression, nbreg();
#   * zero-inflated NB (ZINB), zinbreg(), with the handout's minimum zero
#     part (hospital, chronic, gender, school, insurance);
#   * a sensitivity fit: the same ZINB with `health` added to the zero part.
# The count part is always hospital + health + chronic + gender + school +
# insurance.
#
# Two different "Poisson zeros" appear in the report, and they answer
# different questions:
#   * the UNCONDITIONAL reference: one Poisson distribution with the sample
#     mean for everybody, ignoring the covariates. It predicts
#     n * exp(-mean) zeros. It is a quick descriptive first look only;
#   * the FITTED Poisson regression: each person i has their own mean mu_i
#     from the covariates, and the expected number of zeros is the sum over
#     people of P(Y_i = 0) = exp(-mu_i). This is the model-based number used
#     in the model comparison.
# Covariates spread the means out, and a spread of means produces more zeros
# than a single mean does, so the fitted regression predicts more zeros
# than the unconditional reference (but, for these data, still far too few).
#
# The parametric bootstrap likelihood-ratio (LR) test of NB against ZINB
# (vuong_test(B = ...)) is slow on these data (about 24 s per replicate), so
# it is run once by the script, in chunks, with a checkpoint after each
# chunk; see run_nmes_bootstrap() below for how an interrupted run resumes
# exactly.

bench_file <- system.file("benchmark", "benchmark_functions.R",
                          package = "nbzinbreg")
if (!nzchar(bench_file)) {
  stop("inst/benchmark/benchmark_functions.R was not found in the installed ",
       "nbzinbreg package.")
}
# load_bench_data(), provenance(), check_saved_provenance(), git helpers
sys.source(bench_file, envir = environment(), keep.source = FALSE)

# ---------------------------------------------------------------------------
# Data, formulas and fits
# ---------------------------------------------------------------------------

nmes_data <- function() load_bench_data("NMES1988")

nmes_formulas <- function() {
  list(
    count = visits ~ hospital + health + chronic + gender + school + insurance,
    zinb = visits ~ hospital + health + chronic + gender + school + insurance |
      hospital + chronic + gender + school + insurance,
    zinb_health = visits ~ hospital + health + chronic + gender + school +
      insurance | hospital + health + chronic + gender + school + insurance)
}

# Fit the four models. `maxit` is the EM iteration limit of the ZINB fits
# (the bootstrap refits inherit it). Each fit is checked: the Poisson and NB
# fits must have converged, and the ZINB fits must have converged with a
# passed gradient check. The checks are returned in attr(, "checks"); a fit
# that fails them stops with an error, because every later number depends
# on the fits being maximum likelihood estimates.
nmes_fits <- function(data = nmes_data(), maxit = 5000L) {
  f <- nmes_formulas()
  ctrl <- zinbreg_control(maxit = maxit)
  fits <- list(
    poisson = stats::glm(f$count, family = stats::poisson, data = data),
    nb = nbreg(f$count, data = data),
    zinb = zinbreg(f$zinb, data = data, control = ctrl),
    zinb_health = zinbreg(f$zinb_health, data = data, control = ctrl))
  checks <- data.frame(
    model = names(fits),
    converged = vapply(fits, function(m) isTRUE(m$converged), TRUE),
    gradient_ok = vapply(fits, function(m) {
      if (inherits(m, "zinbreg")) isTRUE(m$gradient_check$ok) else NA
    }, TRUE),
    row.names = NULL)
  bad <- !checks$converged | (!is.na(checks$gradient_ok) & !checks$gradient_ok)
  if (any(bad)) {
    stop("NMES1988 fit(s) did not converge: ",
         paste(checks$model[bad], collapse = ", "))
  }
  attr(fits, "checks") <- checks
  fits
}

# ---------------------------------------------------------------------------
# Descriptive summary and model comparison
# ---------------------------------------------------------------------------

# Descriptive evidence: mean, variance, zeros, and the frequency table of
# visits 0, 1, ..., top - 1 and "top or more", with the frequencies an
# UNCONDITIONAL Poisson distribution (one mean for everybody) would give.
nmes_descriptives <- function(data = nmes_data(), top = 30L) {
  y <- data$visits
  n <- length(y)
  ybar <- mean(y)
  k <- 0:(top - 1L)
  freq <- data.frame(
    visits = c(as.character(k), paste0(top, "+")),
    observed = c(tabulate(y + 1L, nbins = top), sum(y >= top)),
    poisson_unconditional = n * c(stats::dpois(k, ybar),
                                  stats::ppois(top - 1L, ybar,
                                               lower.tail = FALSE)))
  list(
    n = n, mean = ybar, variance = stats::var(y),
    variance_to_mean = stats::var(y) / ybar,
    max = max(y),
    observed_zeros = sum(y == 0), observed_zero_share = mean(y == 0),
    poisson_unconditional_zero_share = exp(-ybar),
    poisson_unconditional_zeros = n * exp(-ybar),
    poisson_unconditional_label =
      "single Poisson with the sample mean, ignoring covariates",
    freq = freq)
}

# Model comparison: df, log-likelihood, AIC, BIC and the expected number of
# zeros under each FITTED model, the sum over people of the model's
# P(Y_i = 0 | x_i):
#   Poisson regression: exp(-mu_i);
#   NB:   (theta / (theta + mu_i))^theta, i.e. dnbinom(0, theta, mu_i);
#   ZINB: pi_i + (1 - pi_i) * dnbinom(0, theta, mu_i), from predict(type = "prob").
# The unconditional Poisson reference is deliberately not in this table.
nmes_comparison <- function(fits, y = fits$nb$y) {
  ll <- lapply(fits, stats::logLik)
  ez <- c(
    poisson = sum(stats::dpois(0, stats::fitted(fits$poisson))),
    nb = sum(stats::dnbinom(0, size = fits$nb$theta,
                            mu = stats::fitted(fits$nb))),
    zinb = sum(stats::predict(fits$zinb, type = "prob")[, 1L]),
    zinb_health = sum(stats::predict(fits$zinb_health, type = "prob")[, 1L]))
  data.frame(
    model = names(fits),
    label = c("Poisson regression", "NB", "ZINB",
              "ZINB, health in zero part")[match(names(fits), names(ez))],
    df = vapply(ll, function(l) as.integer(attr(l, "df")), 1L),
    loglik = vapply(ll, as.numeric, 1),
    AIC = vapply(fits, stats::AIC, 1),
    BIC = vapply(fits, stats::BIC, 1),
    expected_zeros = unname(ez[names(fits)]),
    observed_zeros = sum(y == 0),
    row.names = NULL)
}

# Effects on the multiplicative scale with 95% Wald intervals:
#   count part: incidence rate ratios, IRR = exp(beta);
#   zero part (ZINB only): odds ratios of being a STRUCTURAL zero,
#     OR = exp(gamma). They are not odds ratios of observing a zero.
# The intercepts are left out, and so is the Log(theta) row that the ZINB
# summary lists with the count part. `part` is "count" or "zero".
nmes_effects <- function(fit, part = c("count", "zero"), level = 0.95) {
  part <- match.arg(part)
  s <- summary(fit)$coefficients
  if (inherits(fit, "zinbreg")) {
    tab <- s[[part]]
    prefix <- paste0(part, "_")
  } else {
    if (part == "zero") stop("an nbreg fit has no zero part")
    tab <- s
    prefix <- ""
  }
  keep <- !rownames(tab) %in% c("(Intercept)", "Log(theta)")
  tab <- tab[keep, , drop = FALSE]
  ci <- unclass(stats::confint(fit, parm = paste0(prefix, rownames(tab)),
                               level = level))
  data.frame(
    term = rownames(tab),
    estimate = tab[, "Estimate"],
    ratio = exp(tab[, "Estimate"]),
    lower = exp(ci[, 1L]),
    upper = exp(ci[, 2L]),
    p_value = tab[, "Pr(>|z|)"],
    row.names = NULL)
}

# Observed and expected frequencies of the counts 0, ..., max_count and
# "more than max_count", for any of the four fits, in the layout of the
# package's expected_counts(). The package's function covers the NB and ZINB
# fits; for the Poisson regression (a glm, which expected_counts() does not
# take) the expected frequency of count k is the sum over people of the
# Poisson probability of k at their own fitted mean.
nmes_count_table <- function(fit, max_count = 30L) {
  if (!inherits(fit, "glm")) {
    return(expected_counts(fit, max_count = max_count))
  }
  y <- fit$y
  mu <- stats::fitted(fit)
  k <- 0:max_count
  expected <- vapply(k, function(j) sum(stats::dpois(j, mu)), 1)
  data.frame(
    count = c(as.character(k), paste0(max_count + 1L, "+")),
    observed = as.numeric(c(tabulate(y + 1L, nbins = max_count + 1L),
                            sum(y > max_count))),
    expected = c(expected, sum(stats::ppois(max_count, mu,
                                            lower.tail = FALSE))))
}

# A hanging rootogram drawn from such a table, in the style of the package's
# rootogram(): each bar is the square root of an observed frequency, hung
# from the square root of the expected frequency (the red curve). A bar that
# stops above the dashed zero line marks a count the model over-predicts; a
# bar that goes below it marks a count the model under-predicts. The grouped
# last row is not drawn. `ylim` lets several panels share one scale.
nmes_rootogram_panel <- function(tab, main, ylim = NULL) {
  shown <- seq_len(nrow(tab) - 1L)
  k <- shown - 1L
  top <- sqrt(tab$expected[shown])
  bottom <- top - sqrt(tab$observed[shown])
  if (is.null(ylim)) ylim <- range(0, bottom, top)
  graphics::plot(NA_real_, NA_real_, type = "n", xlab = "Visits",
                 ylab = "sqrt(Frequency)", main = main,
                 xlim = c(-0.5, max(k) + 0.5), ylim = ylim)
  graphics::rect(k - 0.4, bottom, k + 0.4, top, col = "grey80",
                 border = "grey40")
  graphics::abline(h = 0, lty = 2)
  graphics::lines(k, top, type = "b", pch = 19, col = "firebrick", cex = 0.6)
  invisible(tab)
}

# Theta with its Wald interval (on the log scale, transformed back) and its
# profile-likelihood interval.
nmes_theta <- function(fit, level = 0.95) {
  w <- unclass(stats::confint(fit, parm = "theta", level = level,
                              method = "wald"))
  p <- unclass(stats::confint(fit, parm = "theta", level = level,
                              method = "profile"))
  c(theta = fit$theta, se = fit$SE.theta,
    wald_lower = w[1L, 1L], wald_upper = w[1L, 2L],
    profile_lower = p[1L, 1L], profile_upper = p[1L, 2L])
}

# ---------------------------------------------------------------------------
# The bootstrap LR test of NB against ZINB, with checkpoints
# ---------------------------------------------------------------------------
#
# How the run is split and resumed. vuong_test(seed = NULL) draws its
# bootstrap data sets from R's global random number stream and leaves the
# stream advanced, so set.seed(seed) followed by several calls with B = 25
# draws exactly the same data sets, in the same order, as one call with
# B = 999 and seed = seed (the tests check this). After each chunk the run
# saves a checkpoint holding everything produced so far and the state of
# the random number stream (.Random.seed). If the run stops, the same command
# reads the checkpoint, checks that nothing that could change the results
# has changed, restores the stream and continues with the next chunk, so the
# finished replicates are kept unchanged and the remaining ones are the same
# as in an uninterrupted run.

# md5 of a text representation (17 significant digits) of named numbers.
nmes_md5_numbers <- function(x) {
  txt <- paste(names(x), sprintf("%.17g", x), sep = "=", collapse = ";")
  f <- tempfile()
  on.exit(unlink(f))
  writeBin(charToRaw(txt), f)
  unname(tools::md5sum(f))
}

# A fingerprint of the two fitted models the bootstrap uses.
nmes_fit_hash <- function(nb, zinb) {
  nmes_md5_numbers(c(stats::coef(nb), theta_nb = nb$theta,
                     loglik_nb = nb$loglik, stats::coef(zinb),
                     theta_zinb = zinb$theta, loglik_zinb = zinb$loglik,
                     n = nb$nobs))
}

# The settings a checkpoint must share with the run that resumes it.
nmes_boot_header <- function(B, chunk, seed, alpha, maxit, nb, zinb, prov) {
  list(B = as.integer(B), chunk = as.integer(chunk), seed = as.integer(seed),
       alpha = alpha, maxit = as.integer(maxit),
       rng_kind = RNGkind(), r_version = R.version.string,
       fingerprint_installed = prov$fingerprint_installed,
       helper_md5 = prov$helper_md5,
       fit_hash = nmes_fit_hash(nb, zinb))
}

nmes_header_mismatch <- function(saved, now) {
  fields <- names(now)
  fields[!vapply(fields, function(f) identical(saved[[f]], now[[f]]), TRUE)]
}

# Save a checkpoint atomically: write a temporary file, then rename it over
# the old checkpoint, so a crash during the write leaves the previous
# checkpoint intact.
nmes_write_checkpoint <- function(ck, file) {
  tmp <- paste0(file, ".tmp")
  saveRDS(ck, tmp)
  if (!file.rename(tmp, file)) {
    unlink(tmp)
    stop("could not replace the checkpoint ", file)
  }
  invisible(file)
}

nmes_checkpoint_file <- function(out_dir) {
  file.path(out_dir, "nmes_bootstrap_checkpoint.rds")
}

# The bootstrap p-values from the pooled replicates, with exactly the
# formulas of the package's vuong_test() (the tests check they agree). With
# B replicates, B_u usable, F = B - B_u failed and c usable replicates with
# LR* >= the observed LR:
#   p_hat = (1 + c) / (B_u + 1)       (failures dropped)
#   p_lo  = (1 + c) / (B + 1)         (every failure counted below LR)
#   p_hi  = (1 + c + F) / (B + 1)     (every failure counted at or above LR)
# The test is inconclusive if the failures could change the decision at
# alpha (p_lo <= alpha < p_hi), if more than 10% failed, or if fewer than 99
# are usable. An inconclusive test decides nothing (`reject` is NA).
nmes_pool_bootstrap <- function(lr, lr_obs, alpha = 0.05) {
  B <- length(lr)
  usable <- !is.na(lr)
  B_u <- sum(usable)
  F <- B - B_u
  c_ge <- sum(lr[usable] >= lr_obs)
  p_hat <- if (B_u > 0L) (1 + c_ge) / (B_u + 1) else NA_real_
  p_lo <- (1 + c_ge) / (B + 1)
  p_hi <- (1 + c_ge + F) / (B + 1)
  reasons <- c(
    if (p_lo <= alpha && alpha < p_hi) {
      sprintf(paste(
        "the failed replicates could change the decision at alpha = %s",
        "(p between %s and %s)"), format(alpha), format(signif(p_lo, 3)),
        format(signif(p_hi, 3)))
    },
    if (F > 0.1 * B) {
      sprintf("%d of %d replicates (more than 10%%) failed", F, B)
    },
    if (B_u < 99L) {
      sprintf("only %d usable replicates (at least 99 are needed)", B_u)
    })
  list(p.value = p_hat, p.bounds = c(lower = p_lo, upper = p_hi),
       used = B_u, failures = F, exceed = c_ge,
       conclusive = length(reasons) == 0L, reasons = reasons,
       reject = if (length(reasons)) NA else p_hi <= alpha)
}

# Run (or resume) the bootstrap.
#   B, chunk, seed, alpha: the test's settings; maxit: the EM limit of the
#     ZINB fit (the refits inherit it);
#   out_dir: where the checkpoint lives;
#   resume: TRUE continues an existing checkpoint (after checking it
#     matches); FALSE refuses to run if one exists;
#   max_chunks: stop after this many new chunks (for tests and pilots);
#   fits: list(nb =, zinb =) to use instead of the NMES1988 fits (tests);
#   prov: a provenance record to use instead of computing one (tests);
#   maxit_decision: the record of how maxit was chosen (kept in the results).
# Returns the results; `complete` is FALSE when chunks remain.
run_nmes_bootstrap <- function(B = 999L, chunk = 25L, seed = 651L,
                               alpha = 0.05, maxit = 5000L,
                               out_dir = "nmes_boot", source_path = ".",
                               resume = TRUE, max_chunks = Inf, fits = NULL,
                               prov = NULL, maxit_decision = NULL,
                               progress = TRUE) {
  say <- function(...) if (progress) message(...)
  if (is.null(prov)) {
    prov <- provenance(source_path,
                       files = c("benchmark_functions.R", "nmes_functions.R"))
  }
  if (is.null(fits)) {
    f <- nmes_formulas()
    d <- nmes_data()
    fits <- list(nb = nbreg(f$count, data = d),
                 zinb = zinbreg(f$zinb, data = d,
                                control = zinbreg_control(maxit = maxit)))
  }
  header <- nmes_boot_header(B, chunk, seed, alpha, maxit, fits$nb,
                             fits$zinb, prov)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ck_file <- nmes_checkpoint_file(out_dir)
  # A leftover temporary file is an interrupted write: ignore it.
  unlink(paste0(ck_file, ".tmp"))
  genv <- globalenv()
  if (file.exists(ck_file)) {
    if (!resume) {
      stop("a checkpoint already exists in ", out_dir, "; resume it or ",
           "remove the directory first")
    }
    ck <- readRDS(ck_file)
    bad <- nmes_header_mismatch(ck$header, header)
    if (length(bad)) {
      stop("the checkpoint in ", out_dir, " was made with different ",
           paste(bad, collapse = ", "), "; it cannot be resumed by this ",
           "run (the results would not be one bootstrap). Use the same ",
           "settings and installed code, or start in a new directory.")
    }
    ck$resumes <- rbind(ck$resumes, data.frame(
      time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      at_chunk = length(ck$chunks) + 1L))
    say(sprintf("resuming at chunk %d of %d", length(ck$chunks) + 1L,
                ceiling(B / chunk)))
  } else {
    set.seed(seed)
    ck <- list(header = header, provenance = prov,
               rng_start = get(".Random.seed", envir = genv),
               chunks = list(),
               resumes = data.frame(time = character(), at_chunk = integer()),
               started = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    nmes_write_checkpoint(ck, ck_file)
  }
  n_chunks <- as.integer(ceiling(B / chunk))
  done_now <- 0L
  while (length(ck$chunks) < n_chunks && done_now < max_chunks) {
    i <- length(ck$chunks) + 1L
    k <- min(chunk, B - (i - 1L) * chunk)
    # Continue the stream exactly where the last finished chunk left it.
    assign(".Random.seed",
           if (i == 1L) ck$rng_start else ck$chunks[[i - 1L]]$rng_after,
           envir = genv)
    warns <- character()
    t0 <- proc.time()[["elapsed"]]
    v <- withCallingHandlers(
      vuong_test(fits$zinb, fits$nb, B = k, alpha = alpha, seed = NULL),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    b <- v$bootstrap
    ck$chunks[[i]] <- list(
      index = i, from = (i - 1L) * chunk + 1L, to = (i - 1L) * chunk + k,
      replicates = b$replicates, status = b$status, monotone = b$monotone,
      first_failure = b$first_failure, observed = b$statistic,
      loglik = b$loglik, warnings = warns,
      seconds = proc.time()[["elapsed"]] - t0,
      rng_after = get(".Random.seed", envir = genv),
      finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    nmes_write_checkpoint(ck, ck_file)
    done_now <- done_now + 1L
    say(sprintf("chunk %d of %d done: replicates %d-%d, %d failed (%.0f s)",
                i, n_chunks, ck$chunks[[i]]$from, ck$chunks[[i]]$to,
                sum(is.na(b$replicates)), ck$chunks[[i]]$seconds))
  }
  nmes_boot_results(ck, maxit_decision)
}

# Assemble the results from a checkpoint (complete or not).
nmes_boot_results <- function(ck, maxit_decision = NULL) {
  h <- ck$header
  chunks <- ck$chunks
  n_chunks <- as.integer(ceiling(h$B / h$chunk))
  lr <- unlist(lapply(chunks, `[[`, "replicates"))
  reps <- data.frame(
    b = seq_along(lr),
    chunk = rep(seq_along(chunks),
                vapply(chunks, function(x) length(x$replicates), 1L)),
    lr_star = lr, usable = !is.na(lr))
  obs <- if (length(chunks)) chunks[[1L]]$observed else NA_real_
  status <- if (length(chunks)) {
    Reduce(`+`, lapply(chunks, function(x) x$status))
  }
  complete <- length(chunks) == n_chunks && length(lr) == h$B
  list(
    settings = h[c("B", "chunk", "seed", "alpha", "maxit")],
    header = h, complete = complete,
    chunks_done = length(chunks), chunks_total = n_chunks,
    observed = list(statistic = obs,
                    loglik = if (length(chunks)) chunks[[1L]]$loglik),
    observed_consistent = length(unique(vapply(
      chunks, function(x) x$observed, 1))) <= 1L,
    replicates = reps,
    status = status,
    monotone = sum(vapply(chunks, function(x) as.integer(x$monotone), 1L)),
    failure_messages = {
      ff <- vapply(chunks, function(x) as.character(x$first_failure), "")
      unique(ff[!is.na(ff)])
    },
    pool = if (complete) nmes_pool_bootstrap(lr, obs, h$alpha),
    chunk_seconds = vapply(chunks, `[[`, 1, "seconds"),
    chunk_warnings = lapply(chunks, `[[`, "warnings"),
    resumes = ck$resumes, started = ck$started,
    maxit_decision = maxit_decision,
    provenance = ck$provenance)
}

# Two separate verdicts:
#   run_status (computational): COMPLETE when the provenance was verified,
#     every chunk was run, all B replicates were attempted and the observed
#     LR was the same in every chunk; FAILED otherwise. Failed replicates are
#     results of the bootstrap, not a failure of the run.
#   inference (statistical): CONCLUSIVE or INCONCLUSIVE (with the reasons)
#     for a complete run; NOT EVALUATED otherwise.
nmes_bootstrap_status <- function(results) {
  reasons <- c(
    if (results$provenance$status != "verified") {
      paste("provenance unverified:",
            paste(results$provenance$reasons, collapse = "; "))
    },
    if (!results$complete) {
      sprintf(paste("only %d of %d chunks were run; re-run the same command",
                    "to resume from the checkpoint"),
              results$chunks_done, results$chunks_total)
    },
    if (!results$observed_consistent) {
      "the observed LR statistic differed between chunks"
    })
  run_status <- if (length(reasons)) "FAILED" else "COMPLETE"
  inference <- if (run_status != "COMPLETE") "NOT EVALUATED"
               else if (results$pool$conclusive) "CONCLUSIVE"
               else "INCONCLUSIVE"
  list(run_status = run_status, run_reasons = reasons, inference = inference,
       inference_reasons = if (run_status == "COMPLETE") results$pool$reasons)
}

# Save the results: the full object (.rds), one row per replicate (.csv) and
# a status file whose first line is COMPLETE or FAILED and whose second line
# is CONCLUSIVE, INCONCLUSIVE or NOT EVALUATED. Called on every exit.
save_nmes_bootstrap <- function(results, out_dir) {
  st <- nmes_bootstrap_status(results)
  results$run_status <- st$run_status
  results$run_reasons <- st$run_reasons
  results$inference <- st$inference
  results$inference_reasons <- st$inference_reasons
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(results, file.path(out_dir, "nmes_bootstrap_results.rds"))
  utils::write.csv(results$replicates,
                   file.path(out_dir, "nmes_bootstrap_replicates.csv"),
                   row.names = FALSE)
  s <- results$settings
  p <- results$provenance
  pl <- results$pool
  md <- results$maxit_decision
  lines <- c(
    st$run_status,
    if (length(st$run_reasons)) paste("-", st$run_reasons),
    paste0(st$inference,
           if (length(st$inference_reasons))
             paste0(": ", paste(st$inference_reasons, collapse = "; "))),
    "",
    sprintf(paste("Parametric bootstrap LR test of NB against ZINB on",
                  "NMES1988: B = %d in chunks of %d, seed %d, alpha %s,",
                  "ZINB maxit %d."), s$B, s$chunk, s$seed, format(s$alpha),
            s$maxit),
    if (!is.null(md)) sprintf("maxit decision: %s", md$summary),
    sprintf("Observed LR = %.4f; chunks run %d of %d; resumes %d.",
            results$observed$statistic, results$chunks_done,
            results$chunks_total, nrow(results$resumes)),
    if (!is.null(results$status)) sprintf("Replicates: %s; monotone violations %d.",
            paste(names(results$status), as.integer(results$status),
                  sep = " ", collapse = ", "), results$monotone),
    if (!is.null(pl)) sprintf(paste(
      "p = %s (bounds %s to %s counting failures below / at or above);",
      "%d usable, %d failed, %d at or above the observed LR."),
      format(signif(pl$p.value, 4)), format(signif(pl$p.bounds[["lower"]], 4)),
      format(signif(pl$p.bounds[["upper"]], 4)), pl$used, pl$failures,
      pl$exceed),
    sprintf("RNG: %s; compute time %.0f s.",
            paste(results$header$rng_kind, collapse = ", "),
            sum(results$chunk_seconds)),
    "",
    "Provenance:",
    sprintf("- status: %s%s", p$status,
            if (length(p$reasons)) paste0(" (", paste(p$reasons, collapse = "; "), ")") else ""),
    sprintf("- git commit: %s (relevant paths clean: %s)", p$git_commit, p$git_clean),
    sprintf("- code fingerprint (installed / source): %s / %s",
            p$fingerprint_installed, p$fingerprint_source),
    sprintf("- %s on %s", p$r_version, p$platform),
    sprintf("- started %s", results$started))
  writeLines(lines, file.path(out_dir, "nmes_bootstrap_status.txt"),
             useBytes = TRUE)
  invisible(results)
}

# The maxit rule, fixed before the pilot: if more than 5% of the pilot's
# replicates failed at the default maxit, a second pilot with maxit = 20000
# is run, and the setting with the lower failure share is used (the default
# when they are equal). Only failure shares enter the rule, never a p-value
# or an LR statistic, so the choice cannot be steered by the result.
nmes_maxit_decision <- function(pilot_default, pilot_high = NULL,
                                threshold = 0.05, high = 20000L) {
  share <- function(r) mean(!r$replicates$usable)
  s0 <- share(pilot_default)
  if (s0 <= threshold || is.null(pilot_high)) {
    chosen <- pilot_default$settings$maxit
    why <- if (s0 <= threshold) {
      sprintf("default kept: pilot failure share %.3f is at most %.2f",
              s0, threshold)
    } else {
      sprintf("pilot failure share %.3f exceeds %.2f; a maxit = %d pilot is needed",
              s0, threshold, high)
    }
    s1 <- NA_real_
  } else {
    s1 <- share(pilot_high)
    chosen <- if (s1 < s0) pilot_high$settings$maxit else pilot_default$settings$maxit
    why <- sprintf(paste("pilot failure shares %.3f (maxit %d) and %.3f",
                         "(maxit %d); the lower is used, the default on a tie"),
                   s0, pilot_default$settings$maxit, s1,
                   pilot_high$settings$maxit)
  }
  list(maxit = as.integer(chosen), share_default = s0, share_high = s1,
       threshold = threshold,
       needs_high_pilot = s0 > threshold && is.null(pilot_high),
       summary = sprintf("maxit %d (%s)", as.integer(chosen), why))
}

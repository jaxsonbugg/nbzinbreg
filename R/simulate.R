#' Simulate responses from nbreg and zinbreg fits
#'
#' @description
#' Draws new response vectors from a fitted model, keeping the covariates
#' fixed. For an `nbreg` fit each response is a negative binomial draw with
#' the fitted mean and `theta`; for a `zinbreg` fit a structural zero is drawn
#' first with probability \eqn{\pi_i}, and otherwise a negative binomial
#' draw.
#'
#' @details
#' One response is drawn for every row of the fitted data (including rows with
#' weight 0), as `stats::simulate()` does: it simulates the model at each
#' covariate row. Prior weights are not used; a frequency-weighted parametric
#' bootstrap, in which a row of weight `w` stands for `w` independent
#' observations, is what [boot_se()] does.
#'
#' Random numbers follow the conventions of [stats::simulate()]: with
#' `seed = NULL` the current stream is used (and its state recorded); otherwise
#' `set.seed(seed)` is used and the caller's stream is restored afterwards.
#'
#' @section Errors and warnings:
#' `nbzinbreg_error_input`: `nsim` not a whole number of at least 1, a result
#' of more than 1e8 values (`nsim` times the number of rows), a `seed` that
#' `set.seed()` cannot use exactly, or unused arguments.
#'
#' @param object An `nbreg` or `zinbreg` fit.
#' @param nsim Number of response vectors to simulate.
#' @param seed `NULL`, or a whole number for [set.seed()].
#' @param ... Not used (an error if supplied).
#' @return A data frame with one column per simulation (`sim_1`, `sim_2`,
#'   ...) and one row per row of the fitted data, with the attribute `"seed"`
#'   as in [stats::simulate()].
#' @seealso [boot_se()]
#' @examples
#' d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
#' fit <- nbreg(y ~ x, data = d)
#' simulate(fit, nsim = 3, seed = 1)
#' @name simulate-nbzinbreg
NULL

#' @rdname simulate-nbzinbreg
#' @export
simulate.nbreg <- function(object, nsim = 1, seed = NULL, ...) {
  check_no_dots(...length(), ...names(), "simulate",
                c("object", "nsim", "seed"))
  simulate_fit(object, nsim, seed)
}

#' @rdname simulate-nbzinbreg
#' @export
simulate.zinbreg <- function(object, nsim = 1, seed = NULL, ...) {
  check_no_dots(...length(), ...names(), "simulate",
                c("object", "nsim", "seed"))
  simulate_fit(object, nsim, seed)
}

simulate_fit <- function(object, nsim, seed) {
  check_count_arg(nsim, "nsim")
  check_seed(seed)
  n <- length(object$y)
  if (as.numeric(nsim) * n > 1e8) {
    abort_nbz(sprintf(paste(
      "simulate() would create %s values (`nsim` = %s times %d rows), more",
      "than the limit of 1e8 (about 800 MB). Use a smaller `nsim`."),
      format(as.numeric(nsim) * n, big.mark = ",", scientific = FALSE),
      format(nsim, big.mark = ",", scientific = FALSE), n),
      "nbzinbreg_error_input")
  }
  nsim <- as.integer(nsim)
  with_seed(seed, {
    cols <- lapply(seq_len(nsim), function(j) {
      simulate_rows(object, rep(1L, n))$y
    })
    names(cols) <- paste0("sim_", seq_len(nsim))
    out <- as.data.frame(cols)
    rownames(out) <- names(object$y)
    out
  })
}

#' Parametric bootstrap standard errors
#'
#' @description
#' Simulates `R` new data sets from a fitted model, refits the model to each,
#' and reports the standard deviation of the refitted estimates. This is an
#' independent check on the standard errors computed from the information
#' matrix.
#'
#' @details
#' **Bootstrap data.** A row with prior weight `w` stands for `w` independent
#' observations, so `w` responses are drawn for it independently (a row of
#' weight 0 gets none) and the refit uses the distinct values with their
#' multiplicities as weights. This is exactly the bootstrap of the fully
#' expanded data. Weights that are not whole numbers have no such meaning and
#' are refused. The total number of draws per data set is limited to 1e7.
#'
#' **Safeguarded refits.** Each data set is refitted starting from the
#' original estimates. If that refit fails, does not converge, or ends with an
#' estimate on a boundary, the data set is also refitted from the default
#' starting values, and of the refits that converged the one with the higher
#' log-likelihood is used: a boundary estimate reached from one start can be a
#' lower local supremum than an interior maximum reached from another.
#'
#' **The zero-inflation boundary.** For `zinbreg`, a data set may show no
#' zero inflation, so its estimate lies at `pi -> 0`, where EM can stop
#' unconverged (the zero part drifts toward minus infinity and the interior
#' gradient check becomes unavailable) without anything being wrong. If both
#' refits fail, the data set's negative binomial (NB) model is refitted, and
#' the zero-inflated model is refitted once more from that NB fit with `pi`
#' near 0. If that refit also fails, the data set is accepted at the boundary
#' only when all of these hold:
#'
#' * the zero part has an intercept, so `pi -> 0` gives the NB model exactly;
#' * the NB refit converged;
#' * the score at `pi = 0` is not positive, so the likelihood does not rise
#'   into the interior;
#' * no refit reached a log-likelihood above the NB value by more than
#'   `1e-6 (|l| + 1)`, so no interior improvement was found;
#' * the refit from the boundary was monotone and did not fall below the NB
#'   value.
#'
#' The accepted estimate is the boundary point: the count part and `theta`
#' from the NB refit, with every zero coefficient on the boundary. Otherwise
#' the data set is a **numerical failure**.
#'
#' **Boundary data sets are kept.** A simulated data set can have a maximum
#' likelihood estimate on a boundary: a coefficient running off to infinity,
#' a collapsed zero part (no zero inflation, including a simulated response
#' without any zero), or `theta` on a bound. These are valid samples, so the
#' data set is kept and only the affected parameters are left out of that
#' replicate. For each parameter, the bootstrap standard error is the standard
#' deviation over the usable replicates in which it has a finite, identified
#' estimate. If more than 5% of the usable replicates put a parameter on a
#' boundary, its bootstrap distribution has mass at infinity or on a bound,
#' and its standard error is reported as `NA` ("unbounded bootstrap
#' distribution"). Coefficients that already diverge in the original fit get
#' `NA` too.
#'
#' **Numerical failures** are dropped and counted. Up to 10% of `R` gives an
#' `nbzinbreg_warning_bootstrap` warning with the counts (dropping them may
#' bias the standard errors downward, if the failed samples are unusual
#' ones). More than 10% gives an `nbzinbreg_error_fit`.
#'
#' **Parameters.** For `nbreg`: the coefficients and `theta`. For `zinbreg`:
#' the count and zero coefficients (named as by `coef()`) and `log(theta)`.
#' Each `zinbreg` refit is a full EM fit, so `R` refits take some time;
#' `progress = TRUE` reports progress every 10%.
#'
#' @section Errors and warnings:
#' * `nbzinbreg_error_input`: `object` of the wrong class, `R` not a whole
#'   number from 2 to 100,000, an invalid `seed` or `progress`, non-integer or
#'   vanishing frequency weights, a data set of more than 1e7 draws, or
#'   unused arguments.
#' * `nbzinbreg_error_fit`: `object` did not converge (its estimates are not
#'   maximum likelihood estimates, so no bootstrap is run), or more than 10%
#'   of the replicates failed.
#' * `nbzinbreg_warning_bootstrap`: some replicates failed and were dropped.
#'
#' @param object An `nbreg` or `zinbreg` fit.
#' @param R Number of bootstrap replicates.
#' @param seed `NULL`, or a whole number for [set.seed()].
#' @param progress Logical; if `TRUE`, report progress every 10%.
#' @param x An `nbzinbreg_boot` object, for `print()`.
#' @param digits Number of significant digits to print.
#' @param ... Not used (an error if supplied), except for `print()`.
#' @return An object of class `nbzinbreg_boot`, a list with `estimate`,
#'   `se_model` (from the information matrix), `se_boot` (bootstrap), `used`
#'   (replicates used for each parameter), `boundary` (replicates with that
#'   parameter on a boundary), `status` (counts of the replicate outcomes),
#'   `failures`, `boundary_accepted` (replicates accepted at the
#'   zero-inflation boundary after unconverged refits), `R`, `replicates` (the matrix of refitted estimates, `NA`
#'   where left out), `notes`, and the attribute `"seed"`.
#' @seealso [simulate.nbreg()], [confint.nbreg()]
#' @examples
#' d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
#' fit <- nbreg(y ~ x, data = d)
#' boot_se(fit, R = 50, seed = 1)
#' @export
boot_se <- function(object, R = 199, seed = NULL, progress = FALSE, ...) {
  check_fit_class(object, c("nbreg", "zinbreg"))
  check_no_dots(...length(), ...names(), "boot_se",
                c("object", "R", "seed", "progress"))
  check_count_arg(R, "R", min = 2, max = 1e5)
  check_seed(seed)
  check_flag(progress, "progress")
  if (!isTRUE(object$converged)) {
    abort_nbz(sprintf(paste(
      "boot_se() needs a converged fit, but this %s fit did not converge, so",
      "its estimates are not the maximum likelihood estimates: data simulated",
      "from them, and bootstrap standard errors around them, would look",
      "reliable without being so. Refit with a larger `maxit` (see the",
      "control function) or a simpler model first."), class(object)[1L]),
      "nbzinbreg_error_fit")
  }
  w_int <- check_frequency_weights(object$prior.weights, "boot_se()")
  check_bootstrap_size(w_int, "boot_se()")
  R <- as.integer(R)
  zinb <- inherits(object, "zinbreg")
  if (zinb) {
    est <- c(coef(object), `log(theta)` = log(object$theta))
    se_model <- c(sqrt(diag(object$vcov)), `log(theta)` = object$SE.logtheta)
  } else {
    est <- c(object$coefficients, theta = object$theta)
    se_model <- c(sqrt(diag(object$vcov)), theta = object$SE.theta)
  }
  k <- length(est)
  reps <- matrix(NA_real_, R, k, dimnames = list(NULL, names(est)))
  boundary <- matrix(FALSE, R, k, dimnames = list(NULL, names(est)))
  outcome <- character(R)
  accepted <- logical(R)
  first_error <- NA_character_
  step <- max(1L, ceiling(R / 10))

  result <- with_seed(seed, {
    for (b in seq_len(R)) {
      s <- bootstrap_sample(object, w_int)
      r <- boot_refit(object, s)
      outcome[b] <- r$outcome
      accepted[b] <- isTRUE(r$accepted_boundary)
      if (r$outcome == "failure") {
        if (is.na(first_error)) first_error <- r$message
      } else {
        boundary[b, ] <- names(est) %in% r$boundary
        reps[b, ] <- ifelse(boundary[b, ], NA_real_, r$values[names(est)])
      }
      if (progress && (b %% step == 0L || b == R)) {
        message(sprintf("boot_se(): %d of %d replicates done", b, R))
      }
    }
    TRUE
  })

  usable <- outcome != "failure"
  failures <- sum(!usable)
  if (failures > 0.1 * R) {
    abort_nbz(sprintf(paste(
      "boot_se(): %d of %d replicates (%s) failed even after the retries,",
      "more than 10%%, so the bootstrap standard errors would not be",
      "reliable. First failure: %s"), failures, R,
      sprintf("%.0f%%", 100 * failures / R), first_error),
      "nbzinbreg_error_fit")
  }
  if (failures > 0L) {
    warn_nbz(sprintf(paste(
      "boot_se(): %d of %d replicates failed even after the retries and were",
      "dropped. Dropping them can bias the standard errors downward if the",
      "failed samples are unusual. First failure: %s"), failures, R,
      first_error), "nbzinbreg_warning_bootstrap")
  }
  n_usable <- sum(usable)
  n_boundary <- colSums(boundary[usable, , drop = FALSE])
  used <- colSums(!boundary[usable, , drop = FALSE] &
                    is.finite(reps[usable, , drop = FALSE]))
  se_boot <- vapply(seq_len(k), function(j) {
    v <- reps[usable & !boundary[, j], j]
    v <- v[is.finite(v)]
    if (length(v) < 2L) NA_real_ else stats::sd(v)
  }, numeric(1))
  names(se_boot) <- names(est)
  notes <- character()
  unbounded <- n_boundary > 0.05 * n_usable
  if (any(unbounded)) {
    se_boot[unbounded] <- NA_real_
    notes <- c(notes, sprintf(paste(
      "%s: on a boundary in %d of %d usable replicates (more than 5%%), so",
      "the bootstrap distribution is unbounded and no standard error is",
      "given"), names(est)[unbounded], as.integer(n_boundary[unbounded]),
      n_usable))
  }
  diverging <- intersect(names(est), object$diverging)
  if (length(diverging)) {
    se_boot[diverging] <- NA_real_
    notes <- c(notes, sprintf(
      "%s: diverges in the original fit, so no standard error is given",
      diverging))
  }
  n_accepted <- sum(accepted)
  if (n_accepted > 0L) {
    notes <- c(notes, sprintf(paste(
      "%d replicate(s) did not converge in the interior and were accepted at",
      "the zero-inflation boundary (pi -> 0) after the boundary checks"),
      n_accepted))
  }
  structure(list(
    estimate = est, se_model = se_model, se_boot = se_boot, used = used,
    boundary = n_boundary,
    status = table(factor(outcome, levels = c("ok", "boundary", "failure"))),
    failures = failures, boundary_accepted = n_accepted, R = R,
    replicates = reps, notes = notes),
    seed = attr(result, "seed"), class = "nbzinbreg_boot")
}

#' @rdname boot_se
#' @export
print.nbzinbreg_boot <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 ...) {
  cat("\nParametric bootstrap standard errors (", x$R, " replicates: ",
      x$status[["ok"]], " ok, ", x$status[["boundary"]],
      " with a boundary estimate, ", x$failures, " failed)\n\n", sep = "")
  tab <- data.frame(Estimate = x$estimate, `Model SE` = x$se_model,
                    `Bootstrap SE` = x$se_boot,
                    Ratio = x$se_boot / x$se_model, Used = x$used,
                    check.names = FALSE)
  print(tab, digits = digits)
  for (note in x$notes) cat("Note: ", note, ".\n", sep = "")
  invisible(x)
}

# One bootstrap replicate: refit with safeguards, and classify it.
#   outcome "ok": an interior fit; "boundary": a valid sample whose estimate
#   is on a boundary for some parameters (listed in `boundary`); "failure":
#   every refit failed (a numerical failure).
# A simulated response with no zero (zinbreg) or only zeros is a boundary
# sample, not a failure: for no zeros, the count part comes from an NB refit
# and every zero coefficient is on the boundary (pi -> 0); for only zeros,
# every parameter is.
boot_refit <- function(object, s) {
  pos <- s$w > 0
  if (all(s$y[pos] == 0)) {
    names_all <- boot_param_names(object)
    return(list(outcome = "boundary", values = stats::setNames(
      rep(NA_real_, length(names_all)), names_all), boundary = names_all))
  }
  if (inherits(object, "zinbreg") && !any(s$y[pos] == 0)) {
    return(boot_no_zero_replicate(object, s))
  }
  # Starts: the original estimates (fast, and the usual bootstrap practice).
  # If that refit does not end as a regular interior fit (status "ok"), the
  # default starting values are tried too, and of the refits that neither
  # fail nor stop unconverged the one with the higher log-likelihood is used
  # (the first on a tie). A boundary result from one start can be a lower
  # local supremum than an interior maximum found from another, so the first
  # regular refit is not simply taken.
  failed <- function(r) r$status %in% c("error", "not_converged")
  attempts <- list(refit_model(object, s$y, s$w, rows = s$row,
                               start = boot_start(object)))
  if (!identical(attempts[[1L]]$status, "ok")) {
    attempts <- c(attempts, list(
      refit_model(object, s$y, s$w, rows = s$row, start = NULL)))
  }
  regular <- Filter(Negate(failed), attempts)
  if (!length(regular)) {
    if (!inherits(object, "zinbreg")) {
      r <- attempts[[length(attempts)]]
      return(list(outcome = "failure",
                  message = if (!is.na(r$message)) r$message
                  else "the refit did not converge"))
    }
    return(boot_zero_boundary(object, s, attempts))
  }
  ll <- vapply(regular, function(a) a$loglik, numeric(1))
  ll[is.na(ll)] <- -Inf
  r <- regular[[which.max(ll)]]
  b <- boot_boundary(object, r)
  list(outcome = if (length(b)) "boundary" else "ok",
       values = boot_values(object, r), boundary = b)
}

# A zinbreg replicate whose refits all failed or stopped unconverged. Near
# the legitimate boundary pi -> 0 (no zero inflation in this sample), EM can
# stop unconverged without anything being wrong: the zero part drifts toward
# minus infinity, its curvature vanishes, and the interior gradient check
# becomes unavailable (-H not positive definite) or EM runs out of
# iterations. Such a sample has a boundary estimate, and the replicate is
# kept; any other is a numerical failure. The decision, in order:
#   1. The zero part must have an intercept, so that pi -> 0 for every row
#      reproduces the negative binomial (NB) model exactly; otherwise the
#      boundary cannot be identified, and the replicate fails.
#   2. The NB refit of the replicate (same count model) must converge. Its
#      log-likelihood l_NB is the supremum of the ZINB likelihood along
#      pi -> 0, i.e. the value at the boundary.
#   3. A ZINB refit is started at that boundary (count part and theta from
#      the NB refit, pi = 1e-8). If it ends regularly, it is the replicate's
#      fit, as for any other refit.
#   4. Otherwise the boundary is accepted only with all of this evidence:
#      a. first order: the score for moving every pi away from 0 together
#         (the zero intercept, from the boundary) is not positive,
#         D = sum_i w_i e^{o_i} (I(y_i = 0) / f0_i - 1) <= 0, with f0_i the
#         NB probability of a zero and o_i the zero-part offset. D > 0 means
#         the likelihood rises into the interior, so an interior maximum
#         exists that the refits failed to find.
#      b. likelihood: no ZINB refit (any start) found a log-likelihood above
#         l_NB + 1e-6 (|l_NB| + 1). A higher value is an interior improvement
#         that could not be confirmed as a maximum: a failure, not a boundary.
#      c. optimization: the refit from the boundary start was monotone (no
#         recorded violation) and did not fall below l_NB - 1e-6 (|l_NB| + 1);
#         EM from the boundary cannot lose likelihood unless something is
#         numerically wrong.
#   The accepted estimate is the boundary point itself: the count part and
#   theta from the NB refit, every zero coefficient on the boundary (as for a
#   replicate without zeros).
boot_zero_boundary <- function(object, s, attempts) {
  fail <- function(why) list(outcome = "failure", message = paste(
    "the zero-inflated refits did not converge, and the sample is not at",
    "the pi -> 0 boundary:", why))
  pos <- s$w > 0
  Z <- object$z[s$row, , drop = FALSE]
  j <- intercept_index(Z, pos)
  if (j == 0L) return(fail("the zero part has no intercept"))
  X <- object$x[s$row, , drop = FALSE]
  off <- object$offset$count[s$row]
  off_zero <- object$offset$zero[s$row]
  ctl <- object$control$nb_control
  ctl$trace <- FALSE
  nb <- tryCatch(nbreg_fit(X, s$y, off, s$w, ctl, compute_se = FALSE),
                 nbzinbreg_error = function(e) e)
  if (inherits(nb, "condition") || !isTRUE(nb$converged)) {
    return(fail("its negative binomial refit did not converge"))
  }
  g <- stats::setNames(numeric(ncol(Z)), colnames(Z))
  g[j] <- stats::qlogis(1e-8)
  rb <- refit_model(object, s$y, s$w, rows = s$row, start = list(
    count = nb$coefficients, zero = g, theta = nb$theta))
  if (!(rb$status %in% c("error", "not_converged"))) {
    b <- boot_boundary(object, rb)
    return(list(outcome = if (length(b)) "boundary" else "ok",
                values = boot_values(object, rb), boundary = b))
  }
  ev <- boot_boundary_evidence(s, nb, X, off, off_zero, c(attempts, list(rb)),
                               rb)
  if (!ev$accept) return(fail(ev$why))
  out <- boot_nb_boundary_values(object, s, nb, X)
  out$accepted_boundary <- TRUE
  out
}

# The evidence of steps 4a-4c above, for the replicate `s` with NB refit
# `nb`; `attempts` are all ZINB refits, `rb` the one from the boundary start.
boot_boundary_evidence <- function(s, nb, X, off, off_zero, attempts, rb) {
  pos <- s$w > 0
  tol <- 1e-6 * (abs(nb$loglik) + 1)
  log_f0 <- nb_log_f0(lin_pred(X, nb$coefficients, off), nb$theta)
  o <- off_zero - max(off_zero[pos])
  term <- ifelse(s$y == 0, exp(-log_f0), 0) - 1
  D <- sum((s$w * exp(o) * term)[pos])
  if (!is.finite(D) || D > 0) {
    return(list(accept = FALSE, why = sprintf(paste(
      "the score at pi = 0 is positive (D = %.3g), so the likelihood rises",
      "into the interior"), D)))
  }
  ll <- vapply(attempts, function(a) a$loglik, numeric(1))
  ll <- ll[is.finite(ll)]
  if (length(ll) && max(ll) > nb$loglik + tol) {
    return(list(accept = FALSE, why = sprintf(paste(
      "a refit reached a log-likelihood %.3g above the boundary value but",
      "did not converge"), max(ll) - nb$loglik)))
  }
  if (!is.finite(rb$loglik) || rb$loglik < nb$loglik - tol ||
      is.null(rb$engine) || NROW(rb$engine$violations) > 0L) {
    return(list(accept = FALSE, why = paste(
      "the refit from the boundary start failed, lost likelihood, or was",
      "not monotone")))
  }
  list(accept = TRUE, why = NA_character_, D = D)
}

# Replicate values at the pi -> 0 boundary, from an NB fit `fit` of the
# replicate's count model: the count part and theta from the NB fit, every
# zero coefficient on the boundary (NA).
boot_nb_boundary_values <- function(object, s, fit, X) {
  zero_names <- paste0("zero_", colnames(object$z))
  values <- c(stats::setNames(fit$coefficients,
                              paste0("count_", names(fit$coefficients))),
              stats::setNames(rep(NA_real_, length(zero_names)), zero_names),
              `log(theta)` = log(fit$theta))
  div <- union(fit$diverging, structural_divergence(X, s$y, s$w))
  b <- c(zero_names, if (length(div)) paste0("count_", div),
         if (!identical(fit$theta_bound, "none")) "log(theta)")
  list(outcome = "boundary", values = values, boundary = b)
}

boot_param_names <- function(object) {
  if (inherits(object, "zinbreg")) {
    c(names(coef(object)), "log(theta)")
  } else {
    c(names(object$coefficients), "theta")
  }
}

# The first start for a refit: the original estimates.
boot_start <- function(object) {
  if (!inherits(object, "zinbreg")) {
    return(list(beta = object$coefficients, theta = object$theta))
  }
  list(count = object$coefficients$count, zero = object$coefficients$zero,
       theta = object$theta)
}

boot_values <- function(object, r) {
  if (inherits(object, "zinbreg")) {
    cf <- r$coefficients
    c(stats::setNames(cf$count, paste0("count_", names(cf$count))),
      stats::setNames(cf$zero, paste0("zero_", names(cf$zero))),
      `log(theta)` = log(r$theta))
  } else {
    c(r$coefficients, theta = r$theta)
  }
}

# Parameters on a boundary in a refit: diverging coefficients, every zero
# coefficient when the zero part collapsed or is completely separated (the
# engine's flags; see zero_part_separated()), and theta on a bound.
boot_boundary <- function(object, r) {
  eng <- r$engine
  out <- r$diverging
  if (inherits(object, "zinbreg")) {
    if (isTRUE(eng$zero_collapsed) || isTRUE(eng$zero_separated)) {
      out <- union(out, paste0("zero_", colnames(object$z)))
    }
    if (!identical(eng$theta_bound, "none")) out <- c(out, "log(theta)")
  } else if (!identical(eng$theta_bound, "none")) {
    out <- c(out, "theta")
  }
  out
}

# A zinbreg replicate without any zero: the zero part is on its boundary
# (pi -> 0), so the zero coefficients are left out, and the count part and
# theta come from an NB fit of the same count model.
boot_no_zero_replicate <- function(object, s) {
  X <- object$x[s$row, , drop = FALSE]
  ctl <- object$control$nb_control
  ctl$trace <- FALSE
  fit <- tryCatch(
    nbreg_fit(X, s$y, object$offset$count[s$row], s$w, ctl,
              compute_se = FALSE),
    nbzinbreg_error = function(e) e)
  if (inherits(fit, "condition") || !isTRUE(fit$converged)) {
    return(list(outcome = "failure", message = paste(
      "a replicate without zeros, whose negative binomial refit failed")))
  }
  boot_nb_boundary_values(object, s, fit, X)
}

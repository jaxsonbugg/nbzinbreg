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
#' original estimates. If that refit fails or does not converge, it is retried
#' from the default starting values, and for `zinbreg` (with a zero-part
#' intercept) once more from the NB-boundary start (count part at the
#' estimates, `pi` near 0). Only a data set whose refits all fail is a
#' **numerical failure**.
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
#' * `nbzinbreg_error_fit`: more than 10% of the replicates failed.
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
#'   `failures`, `R`, `replicates` (the matrix of refitted estimates, `NA`
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
  first_error <- NA_character_
  step <- max(1L, ceiling(R / 10))

  result <- with_seed(seed, {
    for (b in seq_len(R)) {
      s <- bootstrap_sample(object, w_int)
      r <- boot_refit(object, s)
      outcome[b] <- r$outcome
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
  structure(list(
    estimate = est, se_model = se_model, se_boot = se_boot, used = used,
    boundary = n_boundary,
    status = table(factor(outcome, levels = c("ok", "boundary", "failure"))),
    failures = failures, R = R, replicates = reps, notes = notes),
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
  # Starts, in order: the original estimates (fast, and the usual bootstrap
  # practice), the default starting values, and for zinbreg the NB-boundary
  # start. The first refit that neither fails nor stops unconverged is used.
  starts <- boot_starts(object)
  r <- refit_model(object, s$y, s$w, rows = s$row, start = starts[[1L]])
  for (st in c(list(NULL), starts[-1L])) {
    if (!(r$status %in% c("error", "not_converged"))) break
    r <- refit_model(object, s$y, s$w, rows = s$row, start = st)
  }
  if (r$status %in% c("error", "not_converged")) {
    return(list(outcome = "failure",
                message = if (!is.na(r$message)) r$message
                else "the refit did not converge"))
  }
  b <- boot_boundary(object, r)
  list(outcome = if (length(b)) "boundary" else "ok",
       values = boot_values(object, r), boundary = b)
}

boot_param_names <- function(object) {
  if (inherits(object, "zinbreg")) {
    c(names(coef(object)), "log(theta)")
  } else {
    c(names(object$coefficients), "theta")
  }
}

# Alternative starts for a failed refit: the original estimates, and for
# zinbreg with a zero-part intercept also the NB-boundary start (count part
# at the estimates, zero part at logit(1e-8) on the intercept).
boot_starts <- function(object) {
  if (!inherits(object, "zinbreg")) {
    return(list(list(beta = object$coefficients, theta = object$theta)))
  }
  out <- list(list(count = object$coefficients$count,
                   zero = object$coefficients$zero, theta = object$theta))
  j <- intercept_index(object$z, object$prior.weights > 0)
  if (j > 0L) {
    g <- stats::setNames(numeric(ncol(object$z)), colnames(object$z))
    g[j] <- stats::qlogis(1e-8)
    out <- c(out, list(list(count = object$coefficients$count, zero = g,
                            theta = object$theta)))
  }
  out
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
  zero_names <- paste0("zero_", colnames(object$z))
  values <- c(stats::setNames(fit$coefficients,
                              paste0("count_", names(fit$coefficients))),
              stats::setNames(rep(NA_real_, length(zero_names)), zero_names),
              `log(theta)` = log(fit$theta))
  b <- c(zero_names,
         paste0("count_", union(fit$diverging,
                                structural_divergence(X, s$y, s$w))),
         if (!identical(fit$theta_bound, "none")) "log(theta)")
  list(outcome = "boundary", values = values, boundary = b)
}

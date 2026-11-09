#' Control settings for nbreg()
#'
#' @description
#' Collects and checks the tolerances, iteration limits and other settings
#' used by [nbreg()]. Every argument is validated, and an invalid value stops
#' with an error of class `nbzinbreg_error_input` that names the argument.
#'
#' @param tol Convergence tolerance for the outer loop (one IRLS stage plus
#'   one Newton-Raphson stage per pass). The fit has converged when the
#'   relative change in the log-likelihood,
#'   `abs(new - old) / (abs(old) + 0.1)`, is below `tol` **and** the changes
#'   in the coefficients and in `log(theta)` between passes are below
#'   `sqrt(tol)`.
#' @param maxit Maximum number of outer passes. Reaching it gives a warning
#'   and `converged = FALSE`.
#' @param trace Logical; if `TRUE`, print one line of progress per pass.
#' @param irls_tol Tolerance for the IRLS stage (with `theta` held fixed): it
#'   stops when the relative change in the log-likelihood is below `irls_tol`
#'   and the relative change in the coefficients is below `sqrt(irls_tol)`.
#' @param irls_maxit Maximum number of rounds in one IRLS stage.
#' @param theta_tol Tolerance for the Newton-Raphson stage: it stops when the
#'   step in `log(theta)` is smaller than `theta_tol`.
#' @param theta_maxit Maximum number of rounds in one Newton-Raphson stage.
#' @param theta_start Optional starting value for `theta`. By default a
#'   moment-type estimate from an initial Poisson fit is used.
#' @param theta_min,theta_max Bounds for `theta`. They constrain the
#'   parameter; the likelihood itself is never altered. A fit that ends on a
#'   bound is reported with a warning.
#' @param theta_warn Estimates of `theta` above this value trigger a warning
#'   that the data show little or no overdispersion (a Poisson model may be
#'   adequate). Must satisfy `theta_min < theta_warn <= theta_max`.
#' @param max_halving Maximum number of times a step is halved (backtracking)
#'   when it would lower the log-likelihood or lead to invalid values.
#'
#' @section Provisional:
#' `tol`, `maxit` and `trace` are part of the required interface. The other
#' arguments are provisional and may change in later versions.
#'
#' @return A list of class `nbreg_control` holding the validated settings.
#' @seealso [nbreg()], [zinbreg_control()]
#' @examples
#' nbreg_control()
#' nbreg_control(tol = 1e-10, maxit = 200)
#' @export
nbreg_control <- function(tol = 1e-8, maxit = 100, trace = FALSE,
                          irls_tol = 1e-10, irls_maxit = 50,
                          theta_tol = 1e-10, theta_maxit = 50,
                          theta_start = NULL, theta_min = 1e-8,
                          theta_max = 1e8, theta_warn = 1e4,
                          max_halving = 30) {
  check_tolerance(tol, "tol")
  check_whole(maxit, "maxit")
  check_flag(trace, "trace")
  check_tolerance(irls_tol, "irls_tol")
  check_whole(irls_maxit, "irls_maxit")
  check_tolerance(theta_tol, "theta_tol")
  check_whole(theta_maxit, "theta_maxit")
  check_whole(max_halving, "max_halving")

  check_positive(theta_min, "theta_min")
  check_positive(theta_max, "theta_max")
  if (theta_min >= theta_max) {
    abort_nbz(sprintf(
      "`theta_min` (%s) must be smaller than `theta_max` (%s).",
      format(theta_min), format(theta_max)), "nbzinbreg_error_input")
  }
  check_positive(theta_warn, "theta_warn")
  if (!(theta_warn > theta_min && theta_warn <= theta_max)) {
    abort_nbz(sprintf(paste(
      "`theta_warn` (%s) must lie above `theta_min` (%s) and at most",
      "`theta_max` (%s)."),
      format(theta_warn), format(theta_min), format(theta_max)),
      "nbzinbreg_error_input")
  }
  if (!is.null(theta_start)) {
    check_positive(theta_start, "theta_start")
    if (theta_start < theta_min || theta_start > theta_max) {
      abort_nbz(sprintf(paste(
        "`theta_start` (%s) must lie between `theta_min` (%s) and",
        "`theta_max` (%s)."),
        format(theta_start), format(theta_min), format(theta_max)),
        "nbzinbreg_error_input")
    }
  }

  structure(
    list(tol = tol, maxit = as.integer(round(maxit)), trace = trace,
         irls_tol = irls_tol, irls_maxit = as.integer(round(irls_maxit)),
         theta_tol = theta_tol,
         theta_maxit = as.integer(round(theta_maxit)),
         theta_start = theta_start, theta_min = theta_min,
         theta_max = theta_max, theta_warn = theta_warn,
         max_halving = as.integer(round(max_halving))),
    class = "nbreg_control"
  )
}

#' Control settings for zinbreg()
#'
#' @description
#' Collects and checks the tolerances, iteration limits and other settings
#' used by the EM algorithm in [zinbreg()]. Every argument is validated, and
#' an invalid value stops with an error of class `nbzinbreg_error_input` that
#' names the argument.
#'
#' @param tol Convergence tolerance for the EM algorithm. EM stops when the
#'   relative change in the observed-data log-likelihood,
#'   `abs(new - old) / (abs(old) + 0.1)`, is below `tol` **and** the relative
#'   changes in the count coefficients, the zero coefficients and the change
#'   in `log(theta)` are below `sqrt(tol)`, **and** the gradient check passes
#'   (see `grad_tol`). EM converges slowly when the zero part is poorly
#'   determined, so small changes between iterations do not by themselves
#'   show that the maximum has been reached; the gradient check does.
#' @param maxit Maximum number of EM iterations. Reaching it gives a warning
#'   and `converged = FALSE`.
#' @param trace Logical; if `TRUE`, print one line of progress per EM
#'   iteration.
#' @param start Optional starting values: a list with exactly the elements
#'   `count` (the count coefficients), `zero` (the zero coefficients) and
#'   `theta` (a positive number). Their lengths are checked against the model
#'   in [zinbreg()]. By default the count part starts from a negative
#'   binomial fit and the zero part from a logistic fit of `I(y == 0)`.
#' @param hessian How the standard errors are computed: `"analytic"` (second
#'   derivatives worked out by hand, the default) or `"numeric"` (numerical
#'   differentiation of the log-likelihood with [numDeriv::hessian()]).
#' @param hessian_check Logical; if `TRUE`, also compute the numeric Hessian
#'   and warn if it differs from the analytic one by more than a relative
#'   `1e-4`.
#' @param monotone_tol EM must never lower the log-likelihood, and each
#'   M-step must never lower its own part of the EM objective. A decrease is
#'   reported (with a warning of class `nbzinbreg_warning_monotone`) only if it
#'   is larger than `monotone_tol * (abs(value) + 1)`, which allows for
#'   rounding error.
#' @param grad_tol Threshold for the gradient check. The Newton decrement
#'   \eqn{\delta = g^\top (-H)^{-1} g} of the observed-data log-likelihood
#'   (with `g` its gradient and `H` its Hessian, both analytic) must be at
#'   most `grad_tol`: half of it approximates the log-likelihood that could
#'   still be gained, and every estimate is within \eqn{\sqrt{\delta}}
#'   standard errors of the maximum. EM keeps iterating while the check
#'   fails. The default `1e-8` means within 0.0001 standard errors.
#' @param zero_tol,zero_maxit Tolerance and maximum number of rounds for the
#'   logistic IRLS M-step of the zero part. It stops when the relative change
#'   in its objective is below `zero_tol` and the relative change in the zero
#'   coefficients is below `sqrt(zero_tol)`.
#' @param nb_control Settings for the negative binomial M-step of the count
#'   part, from [nbreg_control()]. Its `theta_min`, `theta_max` and
#'   `theta_warn` also apply to the final estimate of `theta`.
#'
#' @section Provisional:
#' `tol`, `maxit` and `trace` are part of the required interface. The other
#' arguments are provisional and may change in later versions.
#'
#' @return A list of class `zinbreg_control` holding the validated settings.
#' @seealso [zinbreg()], [nbreg_control()]
#' @examples
#' zinbreg_control()
#' zinbreg_control(tol = 1e-10, hessian = "numeric")
#' zinbreg_control(start = list(count = c(0.5, 0.1), zero = -1, theta = 2))
#' @export
zinbreg_control <- function(tol = 1e-8, maxit = 5000, trace = FALSE,
                            start = NULL,
                            hessian = c("analytic", "numeric"),
                            hessian_check = FALSE, monotone_tol = 1e-8,
                            grad_tol = 1e-8, zero_tol = 1e-10,
                            zero_maxit = 50,
                            nb_control = nbreg_control(tol = 1e-10)) {
  check_tolerance(tol, "tol")
  check_whole(maxit, "maxit")
  check_flag(trace, "trace")
  hessian <- check_option(hessian, c("analytic", "numeric"), "hessian")
  check_flag(hessian_check, "hessian_check")
  check_tolerance(monotone_tol, "monotone_tol")
  check_tolerance(grad_tol, "grad_tol")
  check_tolerance(zero_tol, "zero_tol")
  check_whole(zero_maxit, "zero_maxit")
  nb_control <- resolve_control(nb_control, nbreg_control, "nbreg_control",
                                arg = "nb_control")
  start <- check_start(start, nb_control)

  structure(
    list(tol = tol, maxit = as.integer(round(maxit)), trace = trace,
         start = start, hessian = hessian, hessian_check = hessian_check,
         monotone_tol = monotone_tol, grad_tol = grad_tol,
         zero_tol = zero_tol, zero_maxit = as.integer(round(zero_maxit)),
         nb_control = nb_control),
    class = "zinbreg_control"
  )
}

# One value from a fixed set of options, for the control functions. The
# default (the whole set) means the first; partial matching is not used, so a
# misspelling is an error rather than a silent guess.
check_option <- function(x, choices, name) {
  if (identical(x, choices)) return(choices[1L])
  if (!is.character(x) || length(x) != 1L || is.na(x) || !(x %in% choices)) {
    abort_nbz(sprintf("`%s` must be one of %s, not %s.", name,
                      paste0("\"", choices, "\"", collapse = " or "),
                      describe_value(x)), "nbzinbreg_error_input")
  }
  x
}

# Starting values for zinbreg(): NULL, or a list with exactly `count`, `zero`
# and `theta`. Only the structure and the values are checked here; the
# lengths are checked against the model matrices in zinbreg().
check_start <- function(start, nb_control) {
  if (is.null(start)) return(NULL)
  fields <- c("count", "zero", "theta")
  if (!is.list(start) || is.data.frame(start)) {
    abort_nbz(sprintf(paste(
      "`start` must be NULL or a list with elements `count`, `zero` and",
      "`theta`, not %s."), describe_value(start)), "nbzinbreg_error_input")
  }
  nms <- names(start)
  if (length(start) && (is.null(nms) || anyNA(nms) || any(nms == ""))) {
    abort_nbz(paste(
      "Every element of `start` must be named: use",
      "`list(count = ..., zero = ..., theta = ...)`."), "nbzinbreg_error_input")
  }
  dups <- unique(nms[duplicated(nms)])
  if (length(dups)) {
    abort_nbz(sprintf("`start` gives %s more than once.",
                      paste0("`", dups, "`", collapse = ", ")),
              "nbzinbreg_error_input")
  }
  missing_fields <- setdiff(fields, nms)
  extra_fields <- setdiff(nms, fields)
  if (length(missing_fields) || length(extra_fields)) {
    problems <- c(
      if (length(missing_fields)) {
        sprintf("is missing %s", paste0("`", missing_fields, "`",
                                        collapse = ", "))
      },
      if (length(extra_fields)) {
        sprintf("has unknown element(s) %s",
                paste0("`", extra_fields, "`", collapse = ", "))
      })
    abort_nbz(sprintf(paste(
      "`start` %s. It must be a list with exactly the elements `count`,",
      "`zero` and `theta`."), paste(problems, collapse = " and ")),
      "nbzinbreg_error_input")
  }
  for (part in c("count", "zero")) {
    v <- start[[part]]
    if (!is.numeric(v) || length(v) == 0L || is.matrix(v) ||
        any(!is.finite(v))) {
      abort_nbz(sprintf(paste(
        "`start$%s` must be a non-empty numeric vector of finite values, not",
        "%s."), part, describe_value(v)), "nbzinbreg_error_input")
    }
  }
  check_positive(start$theta, "start$theta")
  if (start$theta < nb_control$theta_min ||
      start$theta > nb_control$theta_max) {
    abort_nbz(sprintf(paste(
      "`start$theta` (%s) must lie between `theta_min` (%s) and `theta_max`",
      "(%s) of `nb_control`."), format(start$theta),
      format(nb_control$theta_min), format(nb_control$theta_max)),
      "nbzinbreg_error_input")
  }
  list(count = stats::setNames(as.numeric(start$count), names(start$count)),
       zero = stats::setNames(as.numeric(start$zero), names(start$zero)),
       theta = as.numeric(start$theta))
}

# Accept a control argument as an object made by its constructor (class
# `cls`) or as a named list, and always revalidate it: an object may have been
# edited after it was created (for example ctrl$tol <- -1), so it is never
# trusted as is. `arg` is the argument name used in messages ("control", or
# "nb_control" inside zinbreg_control()). The other model's control object
# (an easy mix-up) gets its own message rather than a list of unknown
# settings.
resolve_control <- function(control, constructor, cls, arg = "control") {
  fun <- cls  # the constructor has the same name as the class
  other <- setdiff(c("nbreg_control", "zinbreg_control"), cls)
  if (!inherits(control, cls) && inherits(control, other)) {
    abort_nbz(sprintf(
      "`%s` was created by `%s()`, but `%s()` is needed here.%s", arg,
      other, fun,
      if (other == "nbreg_control") {
        paste(
          " Settings for the count-part M-step go in",
          "`zinbreg_control(nb_control = nbreg_control(...))`.")
      } else {
        ""
      }), "nbzinbreg_error_input")
  }
  if (inherits(control, cls)) {
    if (!is.list(control)) {
      abort_nbz(sprintf(paste(
        "`%s` has class \"%s\" but is not a list. Create it with `%s(...)`."),
        arg, cls, fun), "nbzinbreg_error_input")
    }
    expected <- names(formals(constructor))
    nms <- names(control)
    if (is.null(nms) || anyNA(nms) || any(nms == "") || anyDuplicated(nms)) {
      abort_nbz(sprintf(paste(
        "`%s` is a malformed %s object: its elements must have unique names.",
        "Create it with `%s(...)`."), arg, cls, fun), "nbzinbreg_error_input")
    }
    missing_fields <- setdiff(expected, nms)
    extra_fields <- setdiff(nms, expected)
    if (length(missing_fields) || length(extra_fields)) {
      problems <- c(
        if (length(missing_fields)) {
          sprintf("missing field(s) %s",
                  paste0("`", missing_fields, "`", collapse = ", "))
        },
        if (length(extra_fields)) {
          sprintf("unknown field(s) %s",
                  paste0("`", extra_fields, "`", collapse = ", "))
        })
      abort_nbz(sprintf(paste(
        "`%s` is a modified %s object with %s. Create it with `%s(...)`",
        "instead of editing it."), arg, cls,
        paste(problems, collapse = " and "), fun), "nbzinbreg_error_input")
    }
    # Re-run every check on the current values.
    return(do.call(constructor, unclass(control)[expected]))
  }
  if (is.list(control)) {
    nms <- names(control)
    if (length(control) && (is.null(nms) || anyNA(nms) || any(nms == ""))) {
      abort_nbz(sprintf(paste(
        "Every element of a `%s` list must be named, for example",
        "`list(tol = 1e-10)`. Prefer `%s = %s(...)`."), arg, arg, fun),
        "nbzinbreg_error_input")
    }
    dups <- unique(nms[duplicated(nms)])
    if (length(dups)) {
      abort_nbz(sprintf(paste(
        "The `%s` list gives setting(s) %s more than once. Give each",
        "setting only once."), arg, paste0("`", dups, "`", collapse = ", ")),
        "nbzinbreg_error_input")
    }
    unknown <- setdiff(nms, names(formals(constructor)))
    if (length(unknown)) {
      abort_nbz(sprintf(paste(
        "Unknown setting(s) in `%s`: %s. See `?%s` for the available",
        "settings."), arg, paste0("`", unknown, "`", collapse = ", "), fun),
        "nbzinbreg_error_input")
    }
    return(do.call(constructor, control))
  }
  abort_nbz(sprintf(
    "`%s` must be created by `%s()` (or be a list), not %s.",
    arg, fun, describe_value(control)), "nbzinbreg_error_input")
}

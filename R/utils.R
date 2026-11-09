# Internal helpers shared across the package.

# Signal that an exported function is still a documented stub.
#
# Every planned function exists with its final (or provisional) signature and
# documentation, but until it is implemented its body calls this helper. The
# error has its own class, `nbzinbreg_not_implemented`, so tests can tell
# "not written yet" apart from any other error. When a function is
# implemented, its call to not_implemented() disappears and its name is
# removed from the `stubbed` list in tests/testthat/test-stubs.R.
not_implemented <- function(what) {
  cnd <- structure(
    class = c("nbzinbreg_not_implemented", "error", "condition"),
    list(
      message = sprintf(
        "`%s()` is not yet implemented in this development version of nbzinbreg.",
        what
      ),
      call = NULL
    )
  )
  stop(cnd)
}

# ---------------------------------------------------------------------------
# Classed errors and warnings
#
# Every error and warning the package raises has its own class, so that tests
# can check exactly which problem was detected and users can catch one kind
# of problem without catching all of them. The classes are:
#
#   errors:   nbzinbreg_error_input  (an argument or the data are invalid)
#             nbzinbreg_error_rank   (collinear model-matrix columns)
#             nbzinbreg_error_fit    (an unrecoverable numerical failure)
#   warnings: nbzinbreg_warning_convergence (an iteration limit was reached)
#             nbzinbreg_warning_theta       (theta very large or on a bound)
#             nbzinbreg_warning_divergence  (coefficients running to -Inf)
#             nbzinbreg_warning_se          (standard errors unavailable)
#
# Each also carries the parent class nbzinbreg_error or nbzinbreg_warning.
# The call is left out (call = NULL) because the messages already say which
# argument or quantity is at fault.
# ---------------------------------------------------------------------------

abort_nbz <- function(message, class) {
  stop(structure(
    class = c(class, "nbzinbreg_error", "error", "condition"),
    list(message = message, call = NULL)
  ))
}

warn_nbz <- function(message, class) {
  warning(structure(
    class = c(class, "nbzinbreg_warning", "warning", "condition"),
    list(message = message, call = NULL)
  ))
}

# A short description of a value for error messages, for example
# `"abc"`, `-1`, `NULL`, or "a numeric vector of length 3".
describe_value <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.function(x)) return("a function")
  if (is.list(x) && !is.data.frame(x)) {
    return(sprintf("a list of length %d", length(x)))
  }
  if (length(x) == 0L) return(sprintf("an empty %s vector", class(x)[1L]))
  if (length(x) > 1L) {
    return(sprintf("a %s vector of length %d", class(x)[1L], length(x)))
  }
  if (is.character(x)) return(sprintf('"%s"', x))
  if (is.factor(x)) return(sprintf('the factor value "%s"', as.character(x)))
  format(x)
}

# The first few values of x, for messages such as
# "Found 3 negative values, for example -1, -2 and -5."
format_values <- function(x, n = 3L) {
  shown <- format(x[seq_len(min(n, length(x)))], digits = 6, trim = TRUE)
  if (length(shown) == 1L) return(shown)
  out <- paste0(paste(shown[-length(shown)], collapse = ", "), " and ",
                shown[length(shown)])
  if (length(x) > n) out <- paste0(out, " (and ", length(x) - n, " more)")
  out
}

# ---------------------------------------------------------------------------
# Validators for single settings (used by the control functions)
# ---------------------------------------------------------------------------

check_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    abort_nbz(sprintf("`%s` must be a single finite number, not %s.",
                      name, describe_value(x)),
              "nbzinbreg_error_input")
  }
  invisible(x)
}

# A tolerance: a single number strictly between 0 and 1.
check_tolerance <- function(x, name) {
  check_number(x, name)
  if (x <= 0 || x >= 1) {
    abort_nbz(sprintf(
      "`%s` must be a number strictly between 0 and 1 (for example 1e-8), not %s.",
      name, describe_value(x)), "nbzinbreg_error_input")
  }
  invisible(x)
}

# A positive number (used for bounds on theta).
check_positive <- function(x, name) {
  check_number(x, name)
  if (x <= 0) {
    abort_nbz(sprintf("`%s` must be a positive number, not %s.",
                      name, describe_value(x)),
              "nbzinbreg_error_input")
  }
  invisible(x)
}

# A whole number of at least `min` (iteration limits).
check_whole <- function(x, name, min = 1) {
  check_number(x, name)
  if (abs(x - round(x)) > sqrt(.Machine$double.eps) || x < min) {
    abort_nbz(sprintf("`%s` must be a whole number of at least %d, not %s.",
                      name, as.integer(min), describe_value(x)),
              "nbzinbreg_error_input")
  }
  invisible(x)
}

# A single TRUE or FALSE.
check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    abort_nbz(sprintf("`%s` must be TRUE or FALSE, not %s.",
                      name, describe_value(x)),
              "nbzinbreg_error_input")
  }
  invisible(x)
}

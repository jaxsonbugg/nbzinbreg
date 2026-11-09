#' Negative binomial regression by IRLS and Newton-Raphson
#'
#' @description
#' Fits a negative binomial (NB2) regression model for count data, with the
#' estimation written from scratch: iteratively reweighted least squares
#' (IRLS) for the regression coefficients, alternated with Newton-Raphson for
#' the dispersion parameter `theta`. No general-purpose optimizer is used.
#' The parameterization is the one used by `MASS::glm.nb()`.
#'
#' @param formula A model formula such as `y ~ x1 + x2`. The response must
#'   contain non-negative whole numbers (counts). `offset()` terms are
#'   allowed.
#' @param data A data frame (or list or environment) containing the variables
#'   in `formula`. If missing, variables are taken from the environment in
#'   which `formula` was created.
#' @param offset Optional numeric vector added to the linear predictor, for
#'   example the log of an exposure time, with one value per row of `data`.
#'   It is added to any `offset()` terms in `formula`. `NULL` (the default)
#'   means no offset argument.
#' @param control A list of settings from [nbreg_control()]. A plain list is
#'   also accepted and is checked by passing it to [nbreg_control()].
#' @param weights Optional non-negative prior weights, one per observation.
#'   Each observation's contribution to the log-likelihood is multiplied by
#'   its weight, so a weight of 2 counts an observation twice and a weight of
#'   0 leaves it out. Fractional weights are allowed.
#' @param subset Optional expression selecting the rows of `data` to use.
#' @param na.action How to handle missing values; defaults to
#'   `getOption("na.action")` (usually [stats::na.omit()]).
#'   [stats::na.exclude()] pads fitted values and residuals with `NA` to the
#'   original length.
#' @param model Logical; if `TRUE` (the default), the model frame is kept in
#'   the result.
#' @param ... Not used. Supplying anything here is an error, so that a
#'   misspelled argument is not silently ignored.
#'
#' @details
#' The model is \eqn{Y_i \sim NB(\mu_i, \theta)} with
#' \eqn{\log \mu_i = x_i^\top \beta + o_i}, so that
#' \eqn{E(Y_i) = \mu_i} and \eqn{Var(Y_i) = \mu_i + \mu_i^2 / \theta}. A small
#' \eqn{\theta} means much more spread than a Poisson model; as
#' \eqn{\theta \to \infty} the model becomes Poisson.
#'
#' The fit alternates two stages until the log-likelihood, the coefficients
#' and \eqn{\log\theta} all stop changing (see [nbreg_control()]):
#' 1. **IRLS for the coefficients, with `theta` held fixed.** Each round
#'    computes the weights \eqn{w_i = a_i \theta \mu_i / (\theta + \mu_i)}
#'    and the working response \eqn{z_i = (\eta_i - o_i) + (y_i - \mu_i)/\mu_i}
#'    and solves a weighted least-squares problem.
#' 2. **Newton-Raphson for `theta`, with the coefficients held fixed.** Steps
#'    are taken on \eqn{\log\theta}, using the slope and curvature of the
#'    log-likelihood (built from `digamma()` and `trigamma()`), with a line
#'    search that guarantees the log-likelihood increases.
#'
#' Every step is checked: a step that would lower the log-likelihood, or lead
#' to infinite or invalid values, is halved until it is acceptable. The
#' log-likelihood itself is always evaluated exactly.
#'
#' Starting values come from a Poisson IRLS fit and a moment-type estimate of
#' `theta`. Standard errors of the coefficients come from
#' \eqn{(X^\top W X)^{-1}} at convergence, treating `theta` as known, as
#' `MASS::glm.nb()` does; the standard error of `theta` comes from the
#' curvature of the log-likelihood in `theta`.
#'
#' Response values within `1e-8` of a whole number are accepted as rounding
#' error and rounded before fitting.
#'
#' @section Errors and warnings:
#' Problems are signalled with classed conditions, so they can be caught
#' selectively. Errors:
#' * `nbzinbreg_error_input`: an invalid argument or invalid data. This covers a
#'   missing, non-formula, one-sided or two-part (`|`) formula; `data` that is
#'   not a data frame, list or environment, or has no rows; unused arguments
#'   in `...`; an invalid `model` flag or `control`; variables that cannot be
#'   found; `offset` or `weights` of the wrong length; missing values with
#'   `na.action = na.fail`; no rows left after `subset` and `na.action`; a
#'   response that is not numeric, has more than one column, or contains
#'   negative, non-integer, missing or infinite values; weights that are not
#'   numeric, negative, infinite or all zero; an offset that is not numeric or
#'   not finite; a response that is zero for every observation; a model with
#'   no coefficients; covariates with infinite values; and too few
#'   observations for the number of coefficients.
#' * `nbzinbreg_error_rank`: columns of the model matrix are linear
#'   combinations of other columns (the message names them).
#' * `nbzinbreg_error_fit`: an unrecoverable numerical failure during fitting.
#'
#' Warnings (the fit is still returned):
#' * `nbzinbreg_warning_convergence`: an iteration or step-halving limit was
#'   reached; `converged` is then `FALSE`.
#' * `nbzinbreg_warning_theta`: `theta` exceeds `theta_warn` (little or no
#'   overdispersion), or ended on one of its bounds.
#' * `nbzinbreg_warning_divergence`: some coefficients run off toward plus or
#'   minus infinity, typically because a group has only zero counts. Their
#'   estimates, standard errors and p-values are not meaningful.
#' * `nbzinbreg_warning_se`: standard errors could not be computed.
#'
#' @return An object of class `nbreg`, a list containing:
#' \describe{
#'   \item{`coefficients`}{Named vector of estimated coefficients.}
#'   \item{`vcov`}{Variance-covariance matrix of the coefficients.}
#'   \item{`theta`, `SE.theta`}{The dispersion estimate and its standard
#'     error.}
#'   \item{`loglik`}{The maximized log-likelihood.}
#'   \item{`df`}{Number of estimated parameters (coefficients plus `theta`).}
#'   \item{`iter`, `converged`}{Number of outer passes, and whether the fit
#'     converged.}
#'   \item{`theta_bound`}{`"none"`, or `"upper"`/`"lower"` if `theta` ended
#'     on a bound.}
#'   \item{`diverging`}{Names of coefficients flagged as diverging.}
#'   \item{`fitted.values`, `linear.predictors`, `residuals`}{Fitted means,
#'     linear predictors (including offsets), and response residuals.}
#'   \item{`prior.weights`, `offset`, `y`, `x`}{The weights, the combined
#'     offset, the (rounded) response and the model matrix.}
#'   \item{`nobs`, `df.residual`}{Number of observations with positive
#'     weight, and that number minus the number of coefficients.}
#'   \item{`loglik_trace`}{Log-likelihood after each outer pass.}
#'   \item{`terms`, `model`, `xlevels`, `contrasts`, `na.action`}{Model
#'     frame information used by the methods.}
#'   \item{`offset_arg_used`}{Whether the `offset` argument was supplied.}
#'   \item{`formula`, `call`, `control`}{The formula, the matched call and
#'     the settings used.}
#' }
#' @seealso [nbreg_control()], [zinbreg()], [nbreg-methods]
#' @examples
#' # The 8-observation example used in the package's design notes
#' d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
#' fit <- nbreg(y ~ x, data = d)
#' fit
#' summary(fit)
#'
#' # Article counts of biochemistry PhD students
#' if (requireNamespace("pscl", quietly = TRUE)) {
#'   data("bioChemists", package = "pscl")
#'   fit2 <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
#'   summary(fit2)
#' }
#' @export
nbreg <- function(formula, data, offset = NULL, control = nbreg_control(),
                  weights = NULL, subset, na.action, model = TRUE, ...) {
  cl <- match.call()

  # ---- arguments --------------------------------------------------------
  if (!missing(formula)) formula <- force_arg(formula, "formula")
  check_formula_arg(formula, missing(formula))
  if (!missing(data)) check_data_arg(force_arg(data, "data"))
  check_dots(...length(), ...names())
  if (!is.logical(model) || length(model) != 1L || is.na(model)) {
    abort_nbz(sprintf("`model` must be TRUE or FALSE, not %s.",
                      describe_value(model)), "nbzinbreg_error_input")
  }
  control <- resolve_nbreg_control(control)

  # ---- model frame ------------------------------------------------------
  mf <- cl[c(1L, match(c("formula", "data", "subset", "weights",
                         "na.action", "offset"), names(cl), 0L))]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- tryCatch(
    eval(mf, parent.frame()),
    error = function(e) {
      abort_nbz(explain_model_frame_error(conditionMessage(e)),
                "nbzinbreg_error_input")
    }
  )
  if (nrow(mf) == 0L) {
    abort_nbz(paste(
      "No observations are left after applying `subset` and `na.action`.",
      "Check the subset condition and the missing values in the data."),
      "nbzinbreg_error_input")
  }
  mt <- attr(mf, "terms")

  # ---- response, weights, offset ---------------------------------------
  y <- check_response(stats::model.response(mf, "any"))
  n <- length(y)
  w <- check_weights(stats::model.weights(mf), n)
  off <- check_offset(
    tryCatch(stats::model.offset(mf), error = function(e) {
      abort_nbz(sprintf("The offset must be a numeric vector. (Original error: %s)",
                        conditionMessage(e)), "nbzinbreg_error_input")
    }),
    n)
  pos <- w > 0
  if (all(y[pos] == 0)) {
    abort_nbz(paste(
      "Every response value (with positive weight) is zero, so the mean",
      "cannot be estimated: its estimate would be zero and the coefficients",
      "would run off to minus infinity."), "nbzinbreg_error_input")
  }

  # ---- model matrix -----------------------------------------------------
  X <- tryCatch(
    stats::model.matrix(mt, mf),
    error = function(e) {
      abort_nbz(explain_model_matrix_error(conditionMessage(e), mf),
                "nbzinbreg_error_input")
    }
  )
  if (ncol(X) == 0L) {
    abort_nbz(paste(
      "The model has no coefficients. Include at least an intercept (for",
      "example `y ~ 1`)."), "nbzinbreg_error_input")
  }
  bad_cols <- colnames(X)[colSums(!is.finite(X)) > 0]
  if (length(bad_cols)) {
    abort_nbz(sprintf(paste(
      "Covariate column(s) %s contain infinite values. Remove or transform",
      "them before fitting."), paste0("`", bad_cols, "`", collapse = ", ")),
      "nbzinbreg_error_input")
  }
  aliased <- aliased_columns(X[pos, , drop = FALSE])
  if (length(aliased)) {
    abort_nbz(sprintf(paste(
      "The model matrix is rank-deficient: %s %s a linear combination of",
      "other columns. Remove %s from the model."),
      paste0("`", aliased, "`", collapse = ", "),
      if (length(aliased) == 1L) "is" else "are",
      if (length(aliased) == 1L) "it" else "them"), "nbzinbreg_error_rank")
  }
  if (sum(pos) <= ncol(X)) {
    abort_nbz(sprintf(paste(
      "Not enough data: %d observation(s) with positive weight for %d",
      "coefficient(s). The negative binomial model also estimates theta, so",
      "more observations than coefficients are needed."),
      sum(pos), ncol(X)), "nbzinbreg_error_input")
  }

  # ---- fit and report ---------------------------------------------------
  fit <- nbreg_fit(X, y, off, w, control)
  diverging <- union(fit$diverging, structural_divergence(X, y, w))
  report_nbreg_fit(fit, diverging, control)

  names(y) <- rownames(mf)
  mu <- stats::setNames(fit$fitted.values, rownames(mf))
  eta <- stats::setNames(fit$linear.predictors, rownames(mf))
  structure(
    list(
      coefficients = fit$coefficients,
      vcov = fit$vcov,
      theta = fit$theta,
      SE.theta = fit$SE.theta,
      loglik = fit$loglik,
      df = ncol(X) + 1L,
      iter = fit$iter,
      converged = fit$converged,
      theta_bound = fit$theta_bound,
      diverging = diverging,
      fitted.values = mu,
      linear.predictors = eta,
      residuals = y - mu,
      prior.weights = w,
      offset = off,
      offset_arg_used = !is.null(cl$offset),
      y = y,
      x = X,
      nobs = sum(pos),
      df.residual = sum(pos) - ncol(X),
      loglik_trace = fit$loglik_trace,
      terms = mt,
      model = if (model) mf,
      xlevels = stats::.getXlevels(mt, mf),
      contrasts = attr(X, "contrasts"),
      na.action = attr(mf, "na.action"),
      formula = formula,
      call = cl,
      control = control
    ),
    class = "nbreg"
  )
}

# ---------------------------------------------------------------------------
# Input checks used by nbreg() (and later by zinbreg()). Each one raises a
# classed error with a message that says what is wrong and how to fix it.
# ---------------------------------------------------------------------------

# Evaluate an argument, turning an evaluation error (typically an object that
# does not exist) into a classed input error.
force_arg <- function(x, name) {
  tryCatch(x, error = function(e) {
    abort_nbz(sprintf("`%s` could not be evaluated: %s", name,
                      conditionMessage(e)), "nbzinbreg_error_input")
  })
}

check_formula_arg <- function(formula, is_missing, fun = "nbreg") {
  if (is_missing) {
    abort_nbz(sprintf(
      "`formula` is missing. Supply a model formula such as `y ~ x1 + x2`."),
      "nbzinbreg_error_input")
  }
  if (is.character(formula)) {
    abort_nbz(sprintf(paste(
      "`formula` must be a formula such as `y ~ x1 + x2`, not a character",
      "string. Use `as.formula(%s)` to convert it."),
      encodeString(formula[1L], quote = "\"")), "nbzinbreg_error_input")
  }
  if (!inherits(formula, "formula")) {
    abort_nbz(sprintf(
      "`formula` must be a formula such as `y ~ x1 + x2`, not %s.",
      describe_value(formula)), "nbzinbreg_error_input")
  }
  if (length(formula) != 3L) {
    abort_nbz(paste(
      "`formula` has no response. Put the count variable on the left of",
      "`~`, as in `y ~ x1 + x2`."), "nbzinbreg_error_input")
  }
  if (fun == "nbreg" && "|" %in% all.names(formula)) {
    abort_nbz(paste(
      "`formula` contains `|`. Two-part formulas (`y ~ x | z`) are for",
      "`zinbreg()`; `nbreg()` takes a one-part formula such as `y ~ x`."),
      "nbzinbreg_error_input")
  }
  invisible(formula)
}

check_data_arg <- function(data) {
  if (is.matrix(data)) {
    abort_nbz(paste(
      "`data` must be a data frame, not a matrix. Convert it with",
      "`as.data.frame()`."), "nbzinbreg_error_input")
  }
  if (!is.data.frame(data) && !is.list(data) && !is.environment(data)) {
    abort_nbz(sprintf(
      "`data` must be a data frame (or a list or environment), not %s.",
      describe_value(data)), "nbzinbreg_error_input")
  }
  if (is.data.frame(data) && nrow(data) == 0L) {
    abort_nbz("`data` has no rows.", "nbzinbreg_error_input")
  }
  invisible(data)
}

check_dots <- function(n_dots, dot_names, fun = "nbreg") {
  if (n_dots == 0L) return(invisible())
  if (is.null(dot_names)) dot_names <- rep("", n_dots)
  dot_names[is.na(dot_names)] <- ""
  named <- dot_names[dot_names != ""]
  msg <- if (length(named)) {
    sprintf("Unused argument(s): %s.", paste0("`", named, "`", collapse = ", "))
  } else {
    sprintf("%d unnamed extra argument(s) were supplied.", n_dots)
  }
  hints <- character()
  if (any(named %in% c("family", "link"))) {
    hints <- c(hints, sprintf(paste(
      "`%s()` always fits a negative binomial (NB2) model with a log link;",
      "it has no `family` or `link` argument."), fun))
  }
  if ("init.theta" %in% named) {
    hints <- c(hints, sprintf(paste(
      "To choose a starting value for theta, use",
      "`control = %s_control(theta_start = ...)`."), fun))
  }
  hints <- c(hints, sprintf(
    "Fitting settings go in `control = %s_control(...)`.", fun))
  abort_nbz(paste(c(msg, hints), collapse = " "), "nbzinbreg_error_input")
}

# Accept `control` as an nbreg_control object or a named list, and always
# revalidate it: an nbreg_control object may have been edited after it was
# created (for example ctrl$tol <- -1), so it is never trusted as is.
resolve_nbreg_control <- function(control) {
  if (inherits(control, "nbreg_control")) {
    if (!is.list(control)) {
      abort_nbz(paste(
        "`control` has class \"nbreg_control\" but is not a list. Create it",
        "with `nbreg_control(...)`."), "nbzinbreg_error_input")
    }
    expected <- names(formals(nbreg_control))
    nms <- names(control)
    if (is.null(nms) || anyNA(nms) || any(nms == "") || anyDuplicated(nms)) {
      abort_nbz(paste(
        "`control` is a malformed nbreg_control object: its elements must",
        "have unique names. Create it with `nbreg_control(...)`."),
        "nbzinbreg_error_input")
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
        "`control` is a modified nbreg_control object with %s. Create it",
        "with `nbreg_control(...)` instead of editing it."),
        paste(problems, collapse = " and ")), "nbzinbreg_error_input")
    }
    # Re-run every check on the current values.
    return(do.call(nbreg_control, unclass(control)[expected]))
  }
  if (is.list(control)) {
    nms <- names(control)
    if (length(control) && (is.null(nms) || anyNA(nms) || any(nms == ""))) {
      abort_nbz(paste(
        "Every element of a `control` list must be named, for example",
        "`list(tol = 1e-10)`. Prefer `control = nbreg_control(...)`."),
        "nbzinbreg_error_input")
    }
    dups <- unique(nms[duplicated(nms)])
    if (length(dups)) {
      abort_nbz(sprintf(paste(
        "The `control` list gives setting(s) %s more than once. Give each",
        "setting only once."), paste0("`", dups, "`", collapse = ", ")),
        "nbzinbreg_error_input")
    }
    unknown <- setdiff(nms, names(formals(nbreg_control)))
    if (length(unknown)) {
      abort_nbz(sprintf(paste(
        "Unknown setting(s) in `control`: %s. See `?nbreg_control` for the",
        "available settings."), paste0("`", unknown, "`", collapse = ", ")),
        "nbzinbreg_error_input")
    }
    return(do.call(nbreg_control, control))
  }
  abort_nbz(sprintf(
    "`control` must be created by `nbreg_control()` (or be a list), not %s.",
    describe_value(control)), "nbzinbreg_error_input")
}

# Turn an error from stats::model.frame() into a message that explains it,
# keeping the original text.
explain_model_frame_error <- function(msg) {
  main <- if (grepl("object '.*' not found", msg)) {
    var <- sub(".*object '([^']*)' not found.*", "\\1", msg)
    sprintf(paste(
      "Variable `%s` was not found in `data` or in the environment where",
      "`formula` was created."), var)
  } else if (grepl("variable lengths differ", msg)) {
    var <- sub(".*found for '([^']*)'.*", "\\1", msg)
    if (identical(var, "(offset)")) {
      "`offset` must have one value per row of `data`."
    } else if (identical(var, "(weights)")) {
      "`weights` must have one value per row of `data`."
    } else {
      sprintf(paste(
        "Variable `%s` has a different length from the other variables in",
        "the model; every variable must have one value per row of `data`."),
        var)
    }
  } else if (grepl("missing values in object", msg)) {
    paste(
      "The data contain missing values, and `na.action = na.fail` does not",
      "allow them. Remove them, or use `na.action = na.omit`.")
  } else if (grepl("weights", msg) && grepl("numeric", msg)) {
    "`weights` must be a numeric vector."
  } else if (grepl("offset", msg) && grepl("numeric", msg)) {
    "`offset` must be a numeric vector."
  } else {
    paste(
      "Could not build the model frame from `formula`, `data` and the other",
      "arguments (for example `subset`).")
  }
  sprintf("%s (Original error: %s)", main, msg)
}

# Turn an error from stats::model.matrix() into a message that explains it.
# The most common one, "contrasts can be applied only to factors with 2 or
# more levels", does not say which variable is at fault, so the factor,
# character and logical covariates with fewer than two distinct values are
# looked up and named.
explain_model_matrix_error <- function(msg, mf) {
  main <- if (grepl("contrasts can be applied only to factors", msg)) {
    covariates <- setdiff(names(mf),
                          c(names(mf)[1L], "(weights)", "(offset)"))
    single <- covariates[vapply(covariates, function(v) {
      x <- mf[[v]]
      (is.factor(x) || is.character(x) || is.logical(x)) &&
        length(unique(x[!is.na(x)])) < 2L
    }, logical(1))]
    if (length(single)) {
      sprintf(paste(
        "Factor(s) %s %s only one level in the data used for fitting (after",
        "`subset` and `na.action`), so %s no contrast to estimate. Remove",
        "%s from the formula, or check `subset`."),
        paste0("`", single, "`", collapse = ", "),
        if (length(single) == 1L) "has" else "have",
        if (length(single) == 1L) "it has" else "they have",
        if (length(single) == 1L) "it" else "them")
    } else {
      paste(
        "A factor in the model has only one level in the data used for",
        "fitting, so it has no contrast to estimate.")
    }
  } else {
    "Could not build the model matrix from `formula`."
  }
  sprintf("%s (Original error: %s)", main, msg)
}

# The response must be a numeric vector of non-negative whole numbers.
# Values within 1e-8 of a whole number are rounded.
check_response <- function(y) {
  if (is.null(y)) {
    abort_nbz("The model has no response variable.", "nbzinbreg_error_input")
  }
  if (is.factor(y) || is.character(y) || is.logical(y)) {
    abort_nbz(sprintf(paste(
      "The response must be numeric counts, not a %s. If it holds counts,",
      "convert it with `as.numeric(as.character(...))`."),
      if (is.factor(y)) "factor" else class(y)[1L]), "nbzinbreg_error_input")
  }
  if (is.matrix(y) || is.data.frame(y)) {
    if (NCOL(y) != 1L) {
      abort_nbz(sprintf(paste(
        "The response must be a single column of counts, not a matrix with",
        "%d columns."), NCOL(y)), "nbzinbreg_error_input")
    }
    y <- as.vector(as.matrix(y))
  }
  if (!is.numeric(y)) {
    abort_nbz(sprintf("The response must be numeric counts, not %s.",
                      describe_value(y)), "nbzinbreg_error_input")
  }
  y <- as.vector(y)
  if (anyNA(y)) {
    abort_nbz(sprintf(paste(
      "The response contains %d missing value(s). Use `na.action = na.omit`",
      "(the default) to leave those rows out."), sum(is.na(y))),
      "nbzinbreg_error_input")
  }
  if (any(is.infinite(y))) {
    abort_nbz(sprintf("The response contains %d infinite value(s).",
                      sum(is.infinite(y))), "nbzinbreg_error_input")
  }
  # Negative beyond rounding error (a value such as -1e-10 is a zero that
  # picked up rounding noise, and is accepted below).
  if (any(y < -1e-8)) {
    neg <- y[y < -1e-8]
    abort_nbz(sprintf(paste(
      "The response must contain non-negative counts. Found %d negative",
      "value(s), for example %s."), length(neg), format_values(neg)),
      "nbzinbreg_error_input")
  }
  frac <- abs(y - round(y)) > 1e-8
  if (any(frac)) {
    abort_nbz(sprintf(paste(
      "The response must contain whole numbers (counts). Found %d",
      "non-integer value(s), for example %s."), sum(frac),
      format_values(y[frac])), "nbzinbreg_error_input")
  }
  round(y) + 0  # "+ 0" turns a rounded -0 into 0
}

check_weights <- function(w, n) {
  if (is.null(w)) return(rep(1, n))
  if (!is.numeric(w)) {
    abort_nbz(sprintf("`weights` must be numeric, not %s.", describe_value(w)),
              "nbzinbreg_error_input")
  }
  w <- as.vector(w)
  if (anyNA(w)) {
    abort_nbz(sprintf(paste(
      "`weights` contains %d missing value(s). Use `na.action = na.omit`",
      "(the default) to leave those rows out."), sum(is.na(w))),
      "nbzinbreg_error_input")
  }
  if (any(is.infinite(w))) {
    abort_nbz(sprintf("`weights` contains %d infinite value(s).",
                      sum(is.infinite(w))), "nbzinbreg_error_input")
  }
  if (any(w < 0)) {
    abort_nbz(sprintf(paste(
      "`weights` must be non-negative. Found %d negative value(s), for",
      "example %s."), sum(w < 0), format_values(w[w < 0])),
      "nbzinbreg_error_input")
  }
  if (all(w == 0)) {
    abort_nbz("All `weights` are zero, so there is nothing to fit.",
              "nbzinbreg_error_input")
  }
  w
}

check_offset <- function(off, n) {
  if (is.null(off)) return(rep(0, n))
  if (!is.numeric(off)) {
    abort_nbz(sprintf("The offset must be numeric, not %s.",
                      describe_value(off)), "nbzinbreg_error_input")
  }
  off <- as.vector(off)
  if (any(!is.finite(off))) {
    abort_nbz(sprintf(paste(
      "The offset contains %d missing or infinite value(s). An offset such",
      "as `log(exposure)` must be finite, so exposures must be positive."),
      sum(!is.finite(off))), "nbzinbreg_error_input")
  }
  off
}

# Turn the engine's flags into classed warnings.
report_nbreg_fit <- function(fit, diverging, control) {
  if (any(fit$issues)) {
    reasons <- c(
      maxit = sprintf("the outer loop reached `maxit` = %d passes",
                      control$maxit),
      irls_maxit = sprintf(
        "the IRLS stage reached `irls_maxit` = %d rounds in the final pass",
        control$irls_maxit),
      theta_maxit = sprintf(paste(
        "the Newton-Raphson stage for theta reached `theta_maxit` = %d",
        "rounds in the final pass"), control$theta_maxit),
      halving = sprintf(paste(
        "a step could not be accepted even after `max_halving` = %d",
        "halvings"), control$max_halving)
    )
    warn_nbz(sprintf(paste(
      "nbreg() did not converge: %s. The estimates may be unreliable;",
      "try raising the limit(s) in `nbreg_control()`."),
      paste(reasons[names(which(fit$issues))], collapse = "; ")),
      "nbzinbreg_warning_convergence")
  }
  if (fit$theta_bound == "upper") {
    warn_nbz(sprintf(paste(
      "theta reached its upper bound `theta_max` = %s: the data show no",
      "overdispersion relative to a Poisson model, so a Poisson model is",
      "adequate. The standard error of theta is not available (NA)."),
      format(control$theta_max)), "nbzinbreg_warning_theta")
  } else if (fit$theta_bound == "lower") {
    warn_nbz(sprintf(paste(
      "theta reached its lower bound `theta_min` = %s, which indicates",
      "extreme overdispersion; the fit may be unreliable. The standard error",
      "of theta is not available (NA)."), format(control$theta_min)),
      "nbzinbreg_warning_theta")
  } else if (fit$theta > control$theta_warn) {
    warn_nbz(sprintf(paste(
      "The estimate of theta (%s) exceeds `theta_warn` = %s: the data show",
      "little or no overdispersion, and a Poisson model may be adequate."),
      format(signif(fit$theta, 6)), format(control$theta_warn)),
      "nbzinbreg_warning_theta")
  }
  if (length(diverging)) {
    warn_nbz(sprintf(paste(
      "Coefficient(s) %s diverge toward plus or minus infinity. The",
      "log-likelihood has converged, but these estimates keep growing because",
      "some observations have fitted means that are effectively zero (for",
      "example a group in which every count is 0). Their estimates, standard",
      "errors and p-values are not meaningful; theta and the other",
      "coefficients can still be interpreted. Consider merging or dropping",
      "that group."), paste0("`", diverging, "`", collapse = ", ")),
      "nbzinbreg_warning_divergence")
  }
  if (fit$vcov_failed) {
    warn_nbz(paste(
      "The covariance matrix of the coefficients could not be computed (the",
      "information matrix is not positive definite); standard errors are NA."),
      "nbzinbreg_warning_se")
  }
  if (fit$se_theta_failed) {
    warn_nbz(paste(
      "The standard error of theta could not be computed (the",
      "log-likelihood is not curved downward in theta at the estimate); it",
      "is NA."), "nbzinbreg_warning_se")
  }
  invisible()
}

#' Methods for nbreg fits
#'
#' @description
#' Standard methods for objects returned by [nbreg()]: printing, summaries
#' with a coefficient table, extracting coefficients, the variance-covariance
#' matrix, the log-likelihood (with a `df` attribute, so `AIC()` and `BIC()`
#' work), fitted values, residuals, the number of observations, and
#' predictions.
#'
#' @param x,object An `nbreg` fit (or, for `print.summary.nbreg()`, its
#'   summary).
#' @param digits Number of significant digits to print.
#' @param type For `residuals()`: `"response"` (\eqn{y - \hat\mu}) or
#'   `"pearson"` (the response residual times the square root of the prior
#'   weight, divided by the estimated standard deviation
#'   \eqn{\sqrt{\hat\mu + \hat\mu^2/\hat\theta}}). For `predict()`: `"link"`
#'   (the linear predictor, on the log scale) or `"response"` (the expected
#'   count).
#' @param newdata Optional data frame of covariate values for `predict()`. If
#'   omitted, predictions are for the data used in the fit.
#' @param se.fit Logical; if `TRUE`, `predict()` also returns standard errors.
#'   On the link scale they are \eqn{\sqrt{x^\top V x}}; on the response scale
#'   they are \eqn{\hat\mu} times that (the delta method). They treat `theta`
#'   as known, like the coefficient standard errors.
#' @param offset For `predict()` with `newdata`, only when the model was
#'   fitted with the `offset` argument: the offset values for the new rows.
#'   Needed when the original offset cannot be recomputed from `newdata`
#'   (for example when it was a stand-alone numeric vector). Offsets written
#'   as `offset()` terms in the formula are always taken from `newdata`.
#' @param ... Not used.
#'
#' @details
#' Fitted values and residuals follow `na.action`: with
#' [stats::na.exclude()] they are padded with `NA` to the length of the
#' original data.
#'
#' For `predict()` with `newdata` and a model fitted with the `offset`
#' argument, the original offset expression is re-evaluated in `newdata` when
#' it refers to at least one of its columns (for example `offset = log(t)`
#' with a column `t`). Otherwise the offset for the new rows must be supplied
#' with `offset =`; it is never silently reused from the fit.
#'
#' @return `coef()` a named vector; `vcov()` a matrix; `logLik()` an object
#'   of class `logLik`; `fitted()` and `residuals()` numeric vectors; `nobs()`
#'   the number of observations with positive weight; `predict()` a numeric
#'   vector, or with `se.fit = TRUE` a list with elements `fit` and `se.fit`;
#'   `summary()` an object of class `summary.nbreg`; the print methods return
#'   their argument invisibly.
#' @name nbreg-methods
#' @seealso [nbreg()]
#' @examples
#' d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
#' fit <- nbreg(y ~ x, data = d)
#' fit
#' summary(fit)
#' coef(fit)
#' AIC(fit)
#' residuals(fit, type = "pearson")
#' predict(fit, newdata = data.frame(x = c(2.5, 9)), type = "response",
#'         se.fit = TRUE)
NULL

#' @rdname nbreg-methods
#' @export
print.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Negative binomial regression (NB2, log link)\n\n")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  ll <- logLik(x)
  cat("\nTheta: ", format(x$theta, digits = digits),
      " (SE ", format(x$SE.theta, digits = digits), ")\n", sep = "")
  cat("Log-likelihood: ", format(as.numeric(ll), digits = max(digits, 6L)),
      " on ", attr(ll, "df"), " df;  AIC: ",
      format(stats::AIC(ll), digits = max(digits, 6L)), "\n", sep = "")
  cat("Observations: ", x$nobs, "\n", sep = "")
  print_fit_notes(x)
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
summary.nbreg <- function(object, ...) {
  est <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z <- est / se
  coefs <- cbind(Estimate = est, `Std. Error` = se, `z value` = z,
                 `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  ll <- logLik(object)
  structure(
    list(call = object$call, coefficients = coefs, theta = object$theta,
         SE.theta = object$SE.theta, loglik = as.numeric(ll),
         df = attr(ll, "df"), aic = stats::AIC(ll), bic = stats::BIC(ll),
         nobs = object$nobs, df.residual = object$df.residual,
         iter = object$iter, converged = object$converged,
         theta_bound = object$theta_bound, diverging = object$diverging,
         pearson = stats::quantile(residuals(object, type = "pearson"),
                                   na.rm = TRUE)),
    class = "summary.nbreg"
  )
}

#' @rdname nbreg-methods
#' @export
print.summary.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Pearson residuals:\n")
  q <- x$pearson
  names(q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(q, digits = digits)
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits, na.print = "NA")
  cat("\nTheta: ", format(x$theta, digits = digits),
      " (SE ", format(x$SE.theta, digits = digits), ")\n", sep = "")
  cat("Log-likelihood: ", format(x$loglik, digits = max(digits, 6L)),
      " on ", x$df, " df\n", sep = "")
  cat("AIC: ", format(x$aic, digits = max(digits, 6L)),
      ";  BIC: ", format(x$bic, digits = max(digits, 6L)), "\n", sep = "")
  cat("Observations: ", x$nobs, ";  residual df: ", x$df.residual, "\n",
      sep = "")
  cat("Outer passes: ", x$iter, ";  converged: ",
      if (x$converged) "yes" else "NO", "\n", sep = "")
  print_fit_notes(x)
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
coef.nbreg <- function(object, ...) {
  object$coefficients
}

#' @rdname nbreg-methods
#' @export
vcov.nbreg <- function(object, ...) {
  object$vcov
}

#' @rdname nbreg-methods
#' @export
logLik.nbreg <- function(object, ...) {
  structure(object$loglik, df = object$df, nobs = object$nobs,
            class = "logLik")
}

#' @rdname nbreg-methods
#' @export
fitted.nbreg <- function(object, ...) {
  stats::napredict(object$na.action, object$fitted.values)
}

#' @rdname nbreg-methods
#' @export
residuals.nbreg <- function(object, type = c("response", "pearson"), ...) {
  type <- check_choice(type, c("response", "pearson"), "type")
  mu <- object$fitted.values
  r <- object$y - mu
  if (type == "pearson") {
    r <- r * sqrt(object$prior.weights) / sqrt(mu + mu^2 / object$theta)
  }
  stats::naresid(object$na.action, r)
}

#' @rdname nbreg-methods
#' @export
nobs.nbreg <- function(object, ...) {
  object$nobs
}

#' @rdname nbreg-methods
#' @export
predict.nbreg <- function(object, newdata, type = c("link", "response"),
                          se.fit = FALSE, offset = NULL, ...) {
  if (...length() > 0L) {
    abort_nbz(sprintf(paste(
      "Unused argument(s) in predict(): %d extra argument(s). The",
      "arguments are `newdata`, `type`, `se.fit` and `offset`."),
      ...length()), "nbzinbreg_error_input")
  }
  type <- check_choice(type, c("link", "response"), "type")
  check_flag(se.fit, "se.fit")

  if (missing(newdata) || is.null(newdata)) {
    if (!is.null(offset)) {
      abort_nbz(paste(
        "`offset` is only used together with `newdata`. Predictions for the",
        "fitted data already include the offset used in the fit."),
        "nbzinbreg_error_input")
    }
    eta <- object$linear.predictors
    X <- object$x
    padded <- function(v) stats::napredict(object$na.action, v)
  } else {
    X <- newdata_matrix(object, newdata)
    eta <- drop(X %*% object$coefficients) +
      newdata_offset(object, newdata, X, offset)
    names(eta) <- rownames(X)
    padded <- function(v) v
  }

  fit <- if (type == "link") eta else exp(eta)
  if (!se.fit) return(padded(fit))
  se_link <- sqrt(rowSums((X %*% object$vcov) * X))
  se <- if (type == "link") se_link else exp(eta) * se_link
  names(se) <- names(fit)
  list(fit = padded(fit), se.fit = padded(se))
}

# ---------------------------------------------------------------------------
# Helpers for the methods
# ---------------------------------------------------------------------------

# Notes printed at the end of print() and print(summary()).
print_fit_notes <- function(x) {
  if (!x$converged) {
    cat("Note: the fit did not converge; see the warning given when fitting.\n")
  }
  if (x$theta_bound != "none") {
    cat("Note: theta is on its ", x$theta_bound,
        " bound; its standard error is not available.\n", sep = "")
  }
  if (length(x$diverging)) {
    cat("Note: diverging coefficient(s) ",
        paste0("`", x$diverging, "`", collapse = ", "),
        "; their estimates and standard errors are not meaningful.\n",
        sep = "")
  }
  invisible()
}

# One value from a fixed set; the default (the whole set) means the first.
check_choice <- function(x, choices, name) {
  if (identical(x, choices)) return(choices[1L])
  if (!is.character(x) || length(x) != 1L || !(x %in% choices)) {
    abort_nbz(sprintf("`%s` must be one of %s, not %s.", name,
                      paste0("\"", choices, "\"", collapse = " or "),
                      describe_value(x)), "nbzinbreg_error_input")
  }
  x
}

# Model matrix for new data, using the factor levels and contrasts of the fit.
newdata_matrix <- function(object, newdata) {
  if (!is.data.frame(newdata)) {
    abort_nbz(sprintf("`newdata` must be a data frame, not %s.",
                      describe_value(newdata)), "nbzinbreg_error_input")
  }
  tt <- stats::delete.response(object$terms)
  mf <- tryCatch(
    stats::model.frame(tt, newdata, na.action = stats::na.pass,
                       xlev = object$xlevels),
    error = function(e) {
      msg <- conditionMessage(e)
      main <- if (grepl("object '.*' not found", msg)) {
        sprintf("Variable `%s` is missing from `newdata`.",
                sub(".*object '([^']*)' not found.*", "\\1", msg))
      } else if (grepl("has new level", msg)) {
        sprintf(paste(
          "`newdata` has factor level(s) that were not present when the",
          "model was fitted: %s."), sub(".*factor ", "factor ", msg))
      } else {
        "Could not build the model frame for `newdata`."
      }
      abort_nbz(sprintf("%s (Original error: %s)", main, msg),
                "nbzinbreg_error_input")
    }
  )
  X <- tryCatch(
    stats::model.matrix(tt, mf, contrasts.arg = object$contrasts),
    error = function(e) {
      abort_nbz(sprintf(
        "Could not build the model matrix for `newdata`. (Original error: %s)",
        conditionMessage(e)), "nbzinbreg_error_input")
    }
  )
  attr(X, "model_frame") <- mf
  X
}

# Offset for new data: offset() terms in the formula are evaluated in
# newdata. If the fit used the `offset` argument, its values for the new rows
# come from `offset =` in predict(), or from re-evaluating the original
# expression in newdata when it refers to a newdata column; otherwise this is
# an error, so the fitted offset is never silently reused.
newdata_offset <- function(object, newdata, X, offset) {
  n_new <- nrow(newdata)
  mf <- attr(X, "model_frame")
  off <- tryCatch(stats::model.offset(mf), error = function(e) {
    abort_nbz(sprintf(
      "The offset() term could not be evaluated in `newdata`. (Original error: %s)",
      conditionMessage(e)), "nbzinbreg_error_input")
  })
  if (is.null(off)) off <- rep(0, n_new)

  if (!isTRUE(object$offset_arg_used)) {
    if (!is.null(offset)) {
      abort_nbz(paste(
        "The model was fitted without an `offset` argument, so `offset`",
        "cannot be used in predict(). Offsets written as `offset()` terms in",
        "the formula are taken from `newdata` automatically."),
        "nbzinbreg_error_input")
    }
    return(off)
  }

  if (!is.null(offset)) {
    if (!is.numeric(offset) || length(offset) != n_new ||
        any(!is.finite(offset))) {
      abort_nbz(sprintf(paste(
        "`offset` must be a finite numeric vector with one value per row of",
        "`newdata` (%d), not %s."), n_new, describe_value(offset)),
        "nbzinbreg_error_input")
    }
    return(off + offset)
  }

  expr <- object$call$offset
  label <- paste(deparse(expr), collapse = " ")
  if (!length(intersect(all.vars(expr), names(newdata)))) {
    abort_nbz(sprintf(paste(
      "The model was fitted with `offset = %s`, which does not refer to any",
      "column of `newdata`, so it cannot be recomputed for the new rows.",
      "Supply the offset for the new rows with `predict(..., offset = )`."),
      label), "nbzinbreg_error_input")
  }
  value <- tryCatch(
    eval(expr, newdata, environment(object$terms)),
    error = function(e) {
      abort_nbz(sprintf(paste(
        "The offset expression `%s` could not be evaluated in `newdata`.",
        "(Original error: %s)"), label, conditionMessage(e)),
        "nbzinbreg_error_input")
    }
  )
  if (!is.numeric(value) || length(value) != n_new) {
    abort_nbz(sprintf(paste(
      "Re-evaluating the offset `%s` in `newdata` gave %s, but one numeric",
      "value per row of `newdata` (%d) is needed. Supply it with",
      "`predict(..., offset = )`."), label, describe_value(value), n_new),
      "nbzinbreg_error_input")
  }
  off + value
}

#' Summaries, Predictions and Other Methods for zinbreg Fits
#'
#' @description
#' Standard methods for objects returned by [zinbreg()]: printing, summaries
#' with coefficient tables for the count and zero parts, extracting
#' coefficients, the variance-covariance matrix, the log-likelihood (with a
#' `df` attribute, so `AIC()` and `BIC()` work), fitted values, residuals, the
#' number of observations, and predictions.
#'
#' @param x,object A `zinbreg` fit (or, for `print.summary.zinbreg()`, its
#'   summary).
#' @param digits Number of significant digits to print.
#' @param model For `coef()` and `vcov()`: `"full"` (both parts, with names
#'   `count_*` and `zero_*`), `"count"` or `"zero"` (one part, without the
#'   prefix), as in `pscl::zeroinfl()`.
#' @param type For `residuals()`: `"response"` (\eqn{y - \hat f}, with the
#'   fitted mean \eqn{\hat f = (1 - \hat\pi)\hat\mu}) or `"pearson"` (the
#'   response residual times the square root of the prior weight, divided by
#'   the estimated standard deviation
#'   \eqn{\sqrt{\hat f (1 + (\hat\pi + 1/\hat\theta)\hat\mu)}}). For
#'   `predict()`: `"response"` (the mean \eqn{(1 - \pi)\mu}), `"count"` (the
#'   count-part mean \eqn{\mu}), `"zero"` (the structural-zero probability
#'   \eqn{\pi}), or `"prob"` (the probabilities \eqn{P(Y = k)} for the counts
#'   in `at`).
#' @param newdata Optional data frame of covariate values for `predict()`. If
#'   omitted, predictions are for the data used in the fit.
#' @param se.fit Logical; if `TRUE`, `predict()` also returns standard errors
#'   by the delta method (not for `type = "prob"`). They use the covariance
#'   matrix of the count and zero coefficients; `theta` does not enter these
#'   predictions. They are `NA` for rows that involve a coefficient without a
#'   finite estimate.
#' @param at For `type = "prob"`: the counts \eqn{k} (non-negative whole
#'   numbers, without duplicates, in any order) at which to evaluate
#'   \eqn{P(Y = k)}; the columns of the result follow this order. The default
#'   is `0:max(y)`, with `y` the response of the fit.
#' @param offset For `predict()` with `newdata`, only when the model was
#'   fitted with the `offset` argument: the count-part offset for the new
#'   rows. Offsets written as `offset()` terms in either part of the formula
#'   are always taken from `newdata`.
#' @param ... Not used.
#'
#' @details
#' Fitted values, residuals and predictions for the fitted data follow
#' `na.action`: with [stats::na.exclude()] they are padded with `NA` to the
#' length of the original data.
#'
#' With `type = "prob"`, the result has one row per observation and one
#' column per value of `at`. To protect against an enormous allocation (for
#' example when the largest count in the data is huge), a result with more
#' than 10 million cells (about 80 MB) is refused with an error that suggests
#' an explicit `at`.
#'
#' For new data, the count-part offset follows the rules of
#' [predict.nbreg()][nbreg-methods]: `offset()` terms are evaluated in
#' `newdata`; an `offset` argument used in the fit is re-evaluated in
#' `newdata` when its expression refers to a `newdata` column, and must
#' otherwise be supplied with `offset =`. All offsets must be finite.
#'
#' @return `coef()` a named vector; `vcov()` a matrix; `logLik()` an object
#'   of class `logLik`; `fitted()` and `residuals()` numeric vectors; `nobs()`
#'   the number of observations with positive weight; `predict()` a numeric
#'   vector (a matrix for `type = "prob"`), or with `se.fit = TRUE` a list
#'   with elements `fit` and `se.fit`; `summary()` an object of class
#'   `summary.zinbreg`; the print methods return their argument invisibly.
#' @name zinbreg-methods
#' @seealso [zinbreg()]
#' @examples
#' if (requireNamespace("pscl", quietly = TRUE)) {
#'   data("bioChemists", package = "pscl")
#'   fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
#'                  data = bioChemists)
#'   fit
#'   summary(fit)
#'   coef(fit, model = "zero")
#'   AIC(fit)
#'   nd <- data.frame(fem = "Women", mar = "Married", kid5 = 0, phd = 3,
#'                    ment = c(0, 10))
#'   predict(fit, newdata = nd, type = "response", se.fit = TRUE)
#'   predict(fit, newdata = nd, type = "prob", at = 0:3)
#' }
NULL

#' @rdname zinbreg-methods
#' @export
print.zinbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Zero-inflated negative binomial regression\n",
      "(count part: NB2 with log link; zero part: logit link)\n\n", sep = "")
  cat("Count coefficients:\n")
  print.default(format(x$coefficients$count, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\nZero coefficients:\n")
  print.default(format(x$coefficients$zero, digits = digits),
                print.gap = 2L, quote = FALSE)
  ll <- logLik(x)
  cat("\nTheta: ", format(x$theta, digits = digits),
      " (SE ", format(x$SE.theta, digits = digits), ")\n", sep = "")
  cat("Log-likelihood: ", format(as.numeric(ll), digits = max(digits, 6L)),
      " on ", attr(ll, "df"), " df;  AIC: ",
      format(stats::AIC(ll), digits = max(digits, 6L)), "\n", sep = "")
  cat("Observations: ", x$nobs, "\n", sep = "")
  cat("EM iterations: ", x$iter, ";  converged: ",
      if (x$converged) "yes" else "NO", "\n", sep = "")
  print_fit_notes(x)
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
summary.zinbreg <- function(object, ...) {
  p <- length(object$coefficients$count)
  q <- length(object$coefficients$zero)
  est <- c(object$coefficients$count, object$coefficients$zero,
           log(object$theta))
  se <- sqrt(diag(object$vcov_full))
  table <- function(i, names) {
    z <- est[i] / se[i]
    out <- cbind(Estimate = est[i], `Std. Error` = se[i], `z value` = z,
                 `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
    rownames(out) <- names
    out
  }
  ll <- logLik(object)
  status <- object$mstep_status
  structure(
    list(call = object$call,
         coefficients = list(
           count = table(c(seq_len(p), p + q + 1L),
                         c(names(object$coefficients$count), "Log(theta)")),
           zero = table(p + seq_len(q), names(object$coefficients$zero))),
         theta = object$theta, SE.theta = object$SE.theta,
         loglik = as.numeric(ll), df = attr(ll, "df"),
         aic = stats::AIC(ll), bic = stats::BIC(ll), nobs = object$nobs,
         df.residual = object$df.residual, iter = object$iter,
         converged = object$converged, theta_bound = object$theta_bound,
         diverging = object$diverging,
         gradient_check = object$gradient_check,
         grad_tol = object$control$grad_tol,
         mstep_status = status[status > 0L],
         em_violations = nrow(object$em_violations),
         start_method = object$start_method,
         pearson = stats::quantile(residuals(object, type = "pearson"),
                                   na.rm = TRUE)),
    class = "summary.zinbreg"
  )
}

#' @rdname zinbreg-methods
#' @export
print.summary.zinbreg <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Pearson residuals:\n")
  q <- x$pearson
  names(q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(q, digits = digits)
  cat("\nCount model coefficients (negative binomial with log link):\n")
  stats::printCoefmat(x$coefficients$count, digits = digits, na.print = "NA")
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  stats::printCoefmat(x$coefficients$zero, digits = digits, na.print = "NA")
  cat("\nTheta: ", format(x$theta, digits = digits),
      " (SE ", format(x$SE.theta, digits = digits), ")\n", sep = "")
  cat("Log-likelihood: ", format(x$loglik, digits = max(digits, 6L)),
      " on ", x$df, " df\n", sep = "")
  cat("AIC: ", format(x$aic, digits = max(digits, 6L)),
      ";  BIC: ", format(x$bic, digits = max(digits, 6L)), "\n", sep = "")
  cat("Observations: ", x$nobs, ";  residual df: ", x$df.residual, "\n",
      sep = "")
  cat("EM iterations: ", x$iter, ";  converged: ",
      if (x$converged) "yes" else "NO", "\n", sep = "")
  gc <- x$gradient_check
  cat("Gradient check: ",
      if (!isTRUE(gc$available)) {
        "unavailable (the negative Hessian is not positive definite)"
      } else {
        sprintf("Newton decrement %s (grad_tol %s)%s",
                format(signif(gc$decrement, 3)), format(x$grad_tol),
                if (isTRUE(gc$ok)) "" else ", NOT passed")
      }, "\n", sep = "")
  cat("Starting values: count part from ", x$start_method[["count"]],
      "; zero part from ", x$start_method[["zero"]], "\n", sep = "")
  if (length(x$mstep_status)) {
    cat("M-step notes (number of EM iterations): ",
        paste(sprintf("%s %d", gsub("_", " ", names(x$mstep_status)),
                      x$mstep_status), collapse = ", "), "\n", sep = "")
  }
  if (x$em_violations) {
    cat("Note: ", x$em_violations, " monotonicity violation(s) were recorded;",
        " see the warning given when fitting.\n", sep = "")
  }
  print_fit_notes(x)
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
coef.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- check_choice(model, c("full", "count", "zero"), "model")
  cf <- object$coefficients
  switch(model,
         full = c(stats::setNames(cf$count, paste0("count_", names(cf$count))),
                  stats::setNames(cf$zero, paste0("zero_", names(cf$zero)))),
         count = cf$count,
         zero = cf$zero)
}

#' @rdname zinbreg-methods
#' @export
vcov.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- check_choice(model, c("full", "count", "zero"), "model")
  vc <- object$vcov
  if (model == "full") return(vc)
  p <- length(object$coefficients$count)
  i <- if (model == "count") seq_len(p) else -seq_len(p)
  out <- vc[i, i, drop = FALSE]
  nms <- names(object$coefficients[[model]])
  dimnames(out) <- list(nms, nms)
  out
}

#' @rdname zinbreg-methods
#' @export
logLik.zinbreg <- function(object, ...) {
  structure(object$loglik, df = object$df, nobs = object$nobs,
            class = "logLik")
}

#' @rdname zinbreg-methods
#' @export
fitted.zinbreg <- function(object, ...) {
  stats::napredict(object$na.action, object$fitted.values)
}

#' @rdname zinbreg-methods
#' @export
residuals.zinbreg <- function(object, type = c("response", "pearson"), ...) {
  type <- check_choice(type, c("response", "pearson"), "type")
  f <- object$fitted.values
  r <- object$y - f
  if (type == "pearson") {
    k <- (object$pi + 1 / object$theta) * object$mu
    v <- f * (1 + k)
    r <- r * sqrt(object$prior.weights) / sqrt(v)
    # The variance underflows to exactly 0 when mu or 1 - pi does
    # (f = (1 - pi) mu), and the formula above would give 0/0 = NaN for a
    # zero count. For those rows only, the residual is computed on the log
    # scale, log f = eta + log(1 - pi) and log v = log f + log1p(k), so
    # nothing is divided by an underflowed variance: a zero count gets
    # -sqrt(f / (1 + k)) (0 to machine precision), and a positive count a
    # large but finite value. The other rows are unchanged.
    under <- !(v > 0)
    if (any(under)) {
      log_f <- object$linear.predictors$count[under] +
        stats::plogis(-object$linear.predictors$zero[under], log.p = TRUE)
      log_v <- log_f + log1p(k[under])
      y_u <- object$y[under]
      r[under] <- sqrt(object$prior.weights[under]) *
        ifelse(y_u == 0, -exp(0.5 * (log_f - log1p(k[under]))),
               (y_u - f[under]) * exp(-0.5 * log_v))
    }
  }
  stats::naresid(object$na.action, r)
}

#' @rdname zinbreg-methods
#' @export
nobs.zinbreg <- function(object, ...) {
  object$nobs
}

#' @rdname zinbreg-methods
#' @export
predict.zinbreg <- function(object, newdata,
                            type = c("response", "count", "zero", "prob"),
                            se.fit = FALSE, at = NULL, offset = NULL, ...) {
  if (...length() > 0L) {
    abort_nbz(sprintf(paste(
      "Unused argument(s) in predict(): %d extra argument(s). The arguments",
      "are `newdata`, `type`, `se.fit`, `at` and `offset`."), ...length()),
      "nbzinbreg_error_input")
  }
  type <- check_choice(type, c("response", "count", "zero", "prob"), "type")
  check_flag(se.fit, "se.fit")
  if (type == "prob" && se.fit) {
    abort_nbz(paste(
      "`se.fit = TRUE` is not available for `type = \"prob\"`. Standard",
      "errors are given for types \"response\", \"count\" and \"zero\"."),
      "nbzinbreg_error_input")
  }
  if (type != "prob" && !is.null(at)) {
    abort_nbz(sprintf(
      "`at` is only used with `type = \"prob\"`, not with `type = \"%s\"`.",
      type), "nbzinbreg_error_input")
  }

  beta <- object$coefficients$count
  gamma <- object$coefficients$zero
  if (missing(newdata) || is.null(newdata)) {
    if (!is.null(offset)) {
      abort_nbz(paste(
        "`offset` is only used together with `newdata`. Predictions for the",
        "fitted data already include the offset used in the fit."),
        "nbzinbreg_error_input")
    }
    X <- object$x
    Z <- object$z
    eta <- object$linear.predictors$count
    zeta <- object$linear.predictors$zero
    padded <- function(v) stats::napredict(object$na.action, v)
  } else {
    X <- newdata_matrix(object, newdata, object$terms$count,
                        object$xlevels$count, object$contrasts$count)
    Z <- newdata_matrix(object, newdata, object$terms$zero,
                        object$xlevels$zero, object$contrasts$zero)
    eta <- drop(X %*% beta) + newdata_offset(
      object, newdata, X, offset, env = environment(object$terms$count),
      terms_label = "The count-part `offset()` term(s) in the formula")
    zeta <- drop(Z %*% gamma) + zero_newdata_offset(Z, nrow(newdata))
    names(eta) <- names(zeta) <- rownames(X)
    padded <- function(v) v
  }

  if (type == "prob") {
    return(padded(zinb_prob_matrix(object, eta, zeta, at)))
  }
  mu <- exp(eta)
  pi <- stats::plogis(zeta)
  fit <- switch(type,
                response = mu * stats::plogis(-zeta),
                count = mu,
                zero = pi)
  if (!se.fit) return(padded(fit))
  se <- zinb_prediction_se(object, X, Z, mu, zeta, type)
  names(se) <- names(fit)
  list(fit = padded(fit), se.fit = padded(se))
}

# ---------------------------------------------------------------------------
# Helpers for the zinbreg methods
# ---------------------------------------------------------------------------

# Zero-part offset for new rows: the offset() terms of the zero part,
# evaluated in newdata (the `offset` argument belongs to the count part).
zero_newdata_offset <- function(Z, n_new) {
  off <- tryCatch(stats::model.offset(attr(Z, "model_frame")),
                  error = function(e) {
                    abort_nbz(sprintf(paste(
                      "The zero-part offset() term could not be evaluated in",
                      "`newdata`. (Original error: %s)"), conditionMessage(e)),
                      "nbzinbreg_error_input")
                  })
  if (is.null(off)) return(rep(0, n_new))
  check_newdata_offset(off, "The zero-part `offset()` term(s) in the formula")
  off
}

# P(Y = k) for every row and every k in `at`:
#   P(Y = 0) = pi + (1 - pi) f_NB(0),   P(Y = k) = (1 - pi) f_NB(k), k >= 1.
# The size of the result is checked before anything is allocated (the
# default at = 0:max(y) is not even formed when it would be too large).
zinb_prob_matrix <- function(object, eta, zeta, at) {
  limit <- 1e7
  n_rows <- length(eta)
  size_error <- function(n_cols, why) {
    abort_nbz(sprintf(paste(
      "`type = \"prob\"` would create a %s x %s matrix (%s cells), more than",
      "the limit of 1e7 cells (about 80 MB). %s Supply only the counts you",
      "need with `at` (for example `at = 0:10`), or predict for fewer rows",
      "with `newdata`."),
      format(n_rows, big.mark = ",", scientific = FALSE),
      format(n_cols, big.mark = ",", scientific = FALSE),
      format(n_rows * n_cols, big.mark = ",", scientific = FALSE), why),
      "nbzinbreg_error_input")
  }
  if (is.null(at)) {
    y_max <- max(object$y)
    if (n_rows * (y_max + 1) > limit) {
      size_error(y_max + 1, sprintf(paste(
        "The default `at = 0:max(y)` uses the largest count in the data",
        "(%s)."), format(y_max, big.mark = ",", scientific = FALSE)))
    }
    at <- 0:y_max
  } else {
    if (!is.numeric(at) || length(at) == 0L || anyNA(at) ||
        any(!is.finite(at))) {
      abort_nbz(sprintf(paste(
        "`at` must be a non-empty numeric vector of counts without missing or",
        "infinite values, not %s."), describe_value(at)),
        "nbzinbreg_error_input")
    }
    if (any(at < 0) || any(abs(at - round(at)) > 1e-8)) {
      bad <- at[at < 0 | abs(at - round(at)) > 1e-8]
      abort_nbz(sprintf(paste(
        "`at` must contain non-negative whole numbers (counts). Found %s."),
        format_values(bad)), "nbzinbreg_error_input")
    }
    at <- round(at)
    if (anyDuplicated(at)) {
      abort_nbz(sprintf("`at` contains duplicated values: %s.",
                        format_values(unique(at[duplicated(at)]))),
                "nbzinbreg_error_input")
    }
    if (n_rows * length(at) > limit) size_error(length(at), "")
  }
  mu <- exp(eta)
  not_pi <- stats::plogis(-zeta)
  pi <- stats::plogis(zeta)
  out <- matrix(0, n_rows, length(at))
  for (j in seq_along(at)) {
    out[, j] <- not_pi * stats::dnbinom(at[j], size = object$theta, mu = mu)
    if (at[j] == 0) out[, j] <- out[, j] + pi
  }
  dimnames(out) <- list(names(eta), as.character(at))
  out
}

# Delta-method standard errors of the predictions, from the covariance
# matrix of the count and zero coefficients:
#   count:    d mu / d beta = mu x
#   zero:     d pi / d gamma = pi (1 - pi) z
#   response: d f / d beta = f x,  d f / d gamma = -f pi z,  f = (1 - pi) mu.
# Coefficients without a finite estimate have NA variances. The prediction
# standard error is NA for exactly the rows that involve them, and is
# computed from the other coefficients for every other row. Whether a row
# involves a coefficient is decided from the model matrices (a nonzero entry
# in its column, for the parts that the prediction type uses), not from the
# derivatives: a derivative that underflows to 0 (for example a fitted mean
# of exp(-800)) must not turn an unknown standard error into a zero. If the
# covariance matrix is entirely unavailable, every standard error is NA.
zinb_prediction_se <- function(object, X, Z, mu, zeta, type) {
  V <- object$vcov
  p <- ncol(X)
  q <- ncol(Z)
  missing_var <- is.na(diag(V))
  if (all(missing_var)) return(rep(NA_real_, nrow(X)))
  V[is.na(V)] <- 0
  pi <- stats::plogis(zeta)
  G <- switch(type,
              count = cbind(mu * X, matrix(0, nrow(X), q)),
              zero = cbind(matrix(0, nrow(Z), p), logit_var(zeta) * Z),
              response = {
                f <- mu * stats::plogis(-zeta)
                cbind(f * X, -f * pi * Z)
              })
  se <- sqrt(pmax(rowSums((G %*% V) * G), 0))
  if (any(missing_var)) {
    used <- switch(type,
                   count = cbind(X, matrix(0, nrow(X), q)),
                   zero = cbind(matrix(0, nrow(Z), p), Z),
                   response = cbind(X, Z))
    involved <- rowSums(used[, missing_var, drop = FALSE] != 0) > 0
    se[involved] <- NA_real_
  }
  se
}

#' Zero-inflated negative binomial regression by EM with IRLS M-steps
#'
#' @description
#' Fits a zero-inflated negative binomial (ZINB) model: each observation is
#' either a "structural zero" (with probability `pi`) or a draw from a
#' negative binomial distribution, which can also produce zeros. The model is
#' fitted with an EM algorithm whose M-steps are IRLS fits, written from
#' scratch. No general-purpose optimizer is used. The parameterization and
#' the coefficient names are those of `pscl::zeroinfl(dist = "negbin")`.
#'
#' @param formula A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model
#'   covariates before the bar, zero-model covariates after it. Use `| 1` for
#'   a zero part with only an intercept. Without a bar (`y ~ x1 + x2`) the
#'   whole right-hand side, including any `offset()` terms, is used for both
#'   parts, as in `pscl::zeroinfl()`; to keep an offset in the count part
#'   only, write it before the bar (`y ~ x + offset(log(t)) | x`). `offset()`
#'   terms are allowed in either part. The response must contain
#'   non-negative whole numbers, with at least one zero.
#' @param data A data frame (or list or environment) containing the variables
#'   in `formula`. If missing, variables are taken from the environment in
#'   which `formula` was created.
#' @param offset Optional numeric vector added to the linear predictor of the
#'   count part, for example the log of an exposure time, with one value per
#'   row of `data`. It is added to any `offset()` terms of the count part.
#'   `NULL` (the default) means no offset argument.
#' @param control A list of settings from [zinbreg_control()]. A plain list is
#'   also accepted and is checked by passing it to [zinbreg_control()].
#' @param weights Optional non-negative prior weights, one per observation.
#'   Each observation's contribution to the log-likelihood is multiplied by
#'   its weight; a weight of 0 leaves it out.
#' @param subset Optional expression selecting the rows of `data` to use.
#' @param na.action How to handle missing values; defaults to
#'   `getOption("na.action")`.
#' @param model Logical; if `TRUE` (the default), the model frame is kept in
#'   the result.
#' @param ... Not used. Supplying anything here is an error, so that a
#'   misspelled argument is not silently ignored.
#'
#' @details
#' The model is
#' \deqn{P(Y_i = 0) = \pi_i + (1 - \pi_i) f_{NB}(0; \mu_i, \theta),}
#' \deqn{P(Y_i = y) = (1 - \pi_i) f_{NB}(y; \mu_i, \theta), \quad y \ge 1,}
#' with \eqn{\log \mu_i = x_i^\top \beta + o_i} and
#' \eqn{\mathrm{logit}\, \pi_i = z_i^\top \gamma + o^z_i}. Because it is a
#' mixture, it is not a GLM, so IRLS cannot be applied to it directly.
#'
#' The EM algorithm treats "is a structural zero" as missing data and
#' repeats:
#' * **E-step:** for each zero, compute the probability \eqn{\tau_i} that it is
#'   a structural zero (it is 0 for every positive count).
#' * **Zero M-step:** a logistic IRLS fit of the fractional response
#'   \eqn{\tau_i} on the zero-part covariates. It never divides by
#'   \eqn{\pi_i(1 - \pi_i)}, so it stays exact when probabilities are
#'   extremely close to 0 or 1, and it is damped when its weighted
#'   cross-product matrix is singular.
#' * **Count M-step:** [nbreg()]'s algorithm with prior weights
#'   \eqn{1 - \tau_i}.
#'
#' Every iteration is checked: EM must never lower the log-likelihood, and
#' neither M-step may lower its own part of the EM objective (generalized
#' EM). The flags of every M-step are recorded.
#'
#' EM converges slowly when the zero part is poorly determined, so small
#' changes between iterations are not enough to stop. EM stops when the
#' changes are small **and** the gradient check passes: the Newton decrement
#' \eqn{g^\top (-H)^{-1} g} of the log-likelihood (analytic gradient and
#' Hessian) must be at most `grad_tol`. The fit is reported as converged only
#' if, in addition, the M-steps of the final iteration converged and no
#' monotonicity violation occurred.
#'
#' Standard errors come from the inverse of the negative Hessian of the
#' observed-data log-likelihood in \eqn{(\beta, \gamma, \log\theta)}: analytic
#' by default, numeric on request (see [zinbreg_control()]).
#'
#' Starting values come from a negative binomial fit for the count part and a
#' logistic fit of `I(y == 0)` for the zero part (with fallbacks if that
#' logistic fit separates), unless supplied with
#' `zinbreg_control(start = )`.
#'
#' Coefficients that run off to plus or minus infinity are detected per part
#' and then held fixed; their standard errors are `NA`. See the section
#' below.
#'
#' @section Errors and warnings:
#' Errors:
#' * `nbzinbreg_error_input`: an invalid argument or invalid data. Besides the
#'   problems listed in [nbreg()] (formula, data, unused arguments, control,
#'   variables, offsets, weights, response values, model matrices), this
#'   covers a formula with more than one response or more than two
#'   right-hand parts, a count or zero part without coefficients, a response
#'   without any zero among the observations with positive weight (there is
#'   then nothing to model with the zero part; use [nbreg()]), a non-finite
#'   zero-part offset, too few observations for the
#'   \eqn{p + q + 1} parameters, and `start` values whose lengths or names do
#'   not match the model.
#' * `nbzinbreg_error_rank`: columns of the count or zero model matrix are
#'   linear combinations of other columns (the message names the part and
#'   the columns).
#' * `nbzinbreg_error_fit`: an unrecoverable numerical failure (the message
#'   gives the EM iteration).
#'
#' Warnings (the fit is still returned):
#' * `nbzinbreg_warning_convergence`: the fit did not converge, with the
#'   reason(s): EM reached `maxit`; the gradient check failed or could not be
#'   carried out; or an M-step of the final iteration did not converge.
#'   `converged` is then `FALSE`.
#' * `nbzinbreg_warning_monotone`: EM lowered the log-likelihood, or an
#'   M-step lowered its part of the EM objective, beyond rounding error.
#' * `nbzinbreg_warning_theta`: `theta` exceeds `theta_warn`, or ended on a
#'   bound (settings of `nb_control`).
#' * `nbzinbreg_warning_divergence`: coefficients without a finite estimate,
#'   explained per part: a zero-part group with no zeros (its structural-zero
#'   probability goes to 0); a group whose counts are all zero (either part
#'   can absorb it, so it is not identified); other coefficients detected
#'   running off to infinity; or no zero inflation at all (every fitted
#'   probability of a structural zero below `1e-8`; [nbreg()] may be
#'   adequate). Their estimates are not meaningful and their standard errors
#'   are `NA`.
#' * `nbzinbreg_warning_se`: standard errors unavailable or possibly
#'   inaccurate (non-finite, asymmetric or ill-conditioned Hessian, a
#'   negative Hessian that is not positive definite, or analytic and numeric
#'   Hessians that disagree with `hessian_check = TRUE`).
#'
#' @return An object of class `zinbreg`, a list containing:
#' \describe{
#'   \item{`coefficients`}{List with the `count` and `zero` coefficient
#'     vectors.}
#'   \item{`vcov`}{Covariance matrix of the count and zero coefficients, with
#'     names `count_*` and `zero_*`.}
#'   \item{`vcov_full`}{Covariance matrix including `Log(theta)`.}
#'   \item{`theta`, `SE.theta`, `SE.logtheta`}{The dispersion estimate, its
#'     standard error (delta method) and the standard error of
#'     `log(theta)`.}
#'   \item{`loglik`, `df`}{The maximized log-likelihood and the number of
#'     estimated parameters, \eqn{p + q + 1}.}
#'   \item{`iter`, `converged`}{Number of EM iterations and whether the fit
#'     converged (see Details).}
#'   \item{`gradient`, `newton_decrement`, `gradient_check`}{The analytic
#'     gradient at the estimates, the Newton decrement, and the result of the
#'     gradient check.}
#'   \item{`mstep_status`}{Counts of M-steps that did not converge, failed to
#'     halve a step, were damped, were restarted and so on, over all
#'     iterations.}
#'   \item{`em_violations`}{Data frame of monotonicity violations (usually
#'     empty).}
#'   \item{`start_method`}{How the starting values were obtained.}
#'   \item{`theta_bound`}{`"none"`, or `"upper"`/`"lower"` if `theta` ended
#'     on a bound.}
#'   \item{`diverging`}{Names (`count_*`, `zero_*`) of coefficients without a
#'     finite estimate.}
#'   \item{`fitted.values`}{The fitted means \eqn{(1 - \pi)\mu}.}
#'   \item{`mu`, `pi`, `linear.predictors`}{Count-part means, structural-zero
#'     probabilities, and the linear predictors of both parts (including
#'     offsets).}
#'   \item{`residuals`}{Response residuals \eqn{y - (1 - \pi)\mu}.}
#'   \item{`prior.weights`, `offset`, `y`, `x`, `z`}{The weights, the offsets
#'     of both parts, the (rounded) response and the two model matrices.}
#'   \item{`nobs`, `df.residual`}{Number of observations with positive
#'     weight, and that number minus `df`.}
#'   \item{`loglik_trace`}{Log-likelihood after each EM iteration.}
#'   \item{`terms`, `model`, `xlevels`, `contrasts`, `na.action`}{Model
#'     frame information used by the methods.}
#'   \item{`hessian`, `hessian_difference`}{The Hessian method used and, with
#'     `hessian_check = TRUE`, the relative difference between the analytic
#'     and numeric Hessians.}
#'   \item{`offset_arg_used`, `formula`, `call`, `control`}{Whether the
#'     `offset` argument was supplied, the two-part formula, the matched call
#'     and the settings used.}
#' }
#' @seealso [zinbreg_control()], [nbreg()], [zinbreg-methods]
#' @examples
#' if (requireNamespace("pscl", quietly = TRUE)) {
#'   data("bioChemists", package = "pscl")
#'   fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
#'                  data = bioChemists)
#'   fit
#'   summary(fit)
#' }
#' @export
zinbreg <- function(formula, data, offset = NULL, control = zinbreg_control(),
                    weights = NULL, subset, na.action, model = TRUE, ...) {
  cl <- match.call()

  # ---- arguments --------------------------------------------------------
  if (!missing(formula)) {
    formula <- force_arg(formula, "formula")
    if (inherits(formula, "Formula")) formula <- stats::formula(formula)
  }
  check_formula_arg(formula, missing(formula), fun = "zinbreg")
  if (!missing(data)) check_data_arg(force_arg(data, "data"))
  check_dots(...length(), ...names(), fun = "zinbreg")
  if (!is.logical(model) || length(model) != 1L || is.na(model)) {
    abort_nbz(sprintf("`model` must be TRUE or FALSE, not %s.",
                      describe_value(model)), "nbzinbreg_error_input")
  }
  control <- resolve_control(control, zinbreg_control, "zinbreg_control")
  ff <- zinb_formula(formula)
  tdata <- if (missing(data)) environment(formula) else data

  # ---- model frame ------------------------------------------------------
  mf <- cl[c(1L, match(c("formula", "data", "subset", "weights",
                         "na.action", "offset"), names(cl), 0L))]
  mf$formula <- ff
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
  part_terms <- function(rhs) {
    tryCatch(stats::terms(ff, data = tdata, rhs = rhs),
             error = function(e) {
               abort_nbz(sprintf(
                 "Could not build the terms of `formula`. (Original error: %s)",
                 conditionMessage(e)), "nbzinbreg_error_input")
             })
  }
  mt_count <- part_terms(1L)
  mt_zero <- stats::delete.response(part_terms(2L))

  # ---- response, weights, offsets ---------------------------------------
  y <- check_response(stats::model.response(mf, "any"))
  n <- length(y)
  w <- check_weights(stats::model.weights(mf), n)
  arg_off <- mf[["(offset)"]]
  if (!is.null(arg_off) && !is.numeric(arg_off)) {
    abort_nbz(sprintf("`offset` must be numeric, not %s.",
                      describe_value(arg_off)), "nbzinbreg_error_input")
  }
  count_terms_off <- terms_offset(mt_count, mf)
  off_count <- check_offset(
    if (is.null(arg_off)) count_terms_off
    else if (is.null(count_terms_off)) arg_off
    else arg_off + count_terms_off, n)
  off_zero <- tryCatch(
    check_offset(terms_offset(mt_zero, mf), n),
    nbzinbreg_error_input = function(e) {
      abort_nbz(sub("^The offset", "The zero-part offset", conditionMessage(e)),
                "nbzinbreg_error_input")
    })

  pos <- w > 0
  if (all(y[pos] == 0)) {
    abort_nbz(paste(
      "Every response value (with positive weight) is zero, so the count",
      "part cannot be estimated: its mean would be zero and the coefficients",
      "would run off to minus infinity."), "nbzinbreg_error_input")
  }
  if (!any(y[pos] == 0)) {
    hidden <- sum(y[!pos] == 0)
    abort_nbz(sprintf(paste(
      "The response has no zeros among the observations with positive",
      "weight%s, so there is no zero inflation to model: the zero-part",
      "coefficients would run off to minus infinity. Use `nbreg()` for a",
      "model without zero inflation."),
      if (hidden) sprintf(" (%d zero(s) have weight 0 and do not count)",
                          hidden) else ""), "nbzinbreg_error_input")
  }

  # ---- model matrices ---------------------------------------------------
  build <- function(mt) {
    tryCatch(stats::model.matrix(mt, mf), error = function(e) {
      abort_nbz(explain_model_matrix_error(conditionMessage(e), mf),
                "nbzinbreg_error_input")
    })
  }
  X <- build(mt_count)
  Z <- build(mt_zero)
  if (ncol(X) == 0L) {
    abort_nbz(paste(
      "The count part has no coefficients. Include at least an intercept",
      "before the bar (for example `y ~ 1 | z`)."), "nbzinbreg_error_input")
  }
  if (ncol(Z) == 0L) {
    abort_nbz(paste(
      "The zero part has no coefficients. Use `| 1` for a zero part with",
      "only an intercept, or `nbreg()` for a model without zero inflation."),
      "nbzinbreg_error_input")
  }
  check_part_matrix(X, pos, "count")
  check_part_matrix(Z, pos, "zero")
  p <- ncol(X)
  q <- ncol(Z)
  if (sum(pos) <= p + q + 1L) {
    abort_nbz(sprintf(paste(
      "Not enough data: %d observation(s) with positive weight for %d",
      "parameters (%d count coefficient(s), %d zero coefficient(s) and",
      "theta). More observations than parameters are needed."),
      sum(pos), p + q + 1L, p, q), "nbzinbreg_error_input")
  }
  if (!is.null(control$start)) {
    control$start$count <- match_start(control$start$count, colnames(X),
                                       "count")
    control$start$zero <- match_start(control$start$zero, colnames(Z), "zero")
  }

  # ---- fit, standard errors, report -------------------------------------
  em <- zinb_em(X, Z, y, off_count, off_zero, w, control)
  se <- zinb_standard_errors(em, X, Z, y, off_count, off_zero, w, control)
  converged <- report_zinbreg_fit(em, se, X, Z, control)

  nms <- zinb_param_names(X, Z)
  k <- p + q
  vcov_full <- se$vcov
  vcov <- vcov_full[seq_len(k), seq_len(k), drop = FALSE]
  se_logtheta <- sqrt(vcov_full[k + 1L, k + 1L])
  rn <- rownames(mf)
  mu <- stats::setNames(exp(em$eta), rn)
  pi <- stats::setNames(stats::plogis(em$zeta), rn)
  fitted <- stats::setNames(exp(em$eta) * stats::plogis(-em$zeta), rn)
  names(y) <- rn
  structure(
    list(
      coefficients = list(count = em$beta, zero = em$gamma),
      vcov = vcov,
      vcov_full = vcov_full,
      theta = em$theta,
      SE.theta = em$theta * se_logtheta,
      SE.logtheta = se_logtheta,
      loglik = em$loglik,
      df = k + 1L,
      iter = em$iter,
      converged = converged,
      gradient = em$gradient$gradient,
      newton_decrement = em$gradient$decrement,
      gradient_check = em$gradient[c("available", "ok", "decrement",
                                     "max_scaled")],
      mstep_status = em$mstep_status,
      em_violations = em$violations,
      start_method = em$start_method,
      theta_bound = em$theta_bound,
      diverging = setdiff(em$gradient$excluded, "Log(theta)"),
      fitted.values = fitted,
      mu = mu,
      pi = pi,
      linear.predictors = list(count = stats::setNames(em$eta, rn),
                               zero = stats::setNames(em$zeta, rn)),
      residuals = y - fitted,
      prior.weights = w,
      offset = list(count = off_count, zero = off_zero),
      offset_arg_used = !is.null(cl$offset),
      y = y,
      x = X,
      z = Z,
      nobs = sum(pos),
      df.residual = sum(pos) - (k + 1L),
      loglik_trace = em$loglik_trace,
      terms = list(count = mt_count, zero = mt_zero,
                   full = attr(mf, "terms")),
      model = if (model) mf,
      xlevels = list(count = stats::.getXlevels(mt_count, mf),
                     zero = stats::.getXlevels(mt_zero, mf)),
      contrasts = list(count = attr(X, "contrasts"),
                       zero = attr(Z, "contrasts")),
      na.action = attr(mf, "na.action"),
      hessian = control$hessian,
      hessian_difference = se$difference,
      formula = ff,
      call = cl,
      control = control,
      param_names = nms
    ),
    class = "zinbreg"
  )
}

# ---------------------------------------------------------------------------
# Helpers for zinbreg()
# ---------------------------------------------------------------------------

# Parse the formula into a two-part Formula. One part (y ~ x) means the same
# right-hand side, offset() terms included, for both parts, as in
# pscl::zeroinfl().
zinb_formula <- function(formula) {
  ff <- tryCatch(Formula::as.Formula(formula), error = function(e) {
    abort_nbz(sprintf(
      "`formula` could not be read as a two-part formula. (Original error: %s)",
      conditionMessage(e)), "nbzinbreg_error_input")
  })
  len <- length(ff)
  if (len[1L] != 1L) {
    abort_nbz(sprintf(paste(
      "`formula` must have exactly one response on the left of `~`, not %d",
      "parts separated by `|`."), len[1L]), "nbzinbreg_error_input")
  }
  if (len[2L] > 2L) {
    abort_nbz(sprintf(paste(
      "`formula` has %d parts on the right of `~`; zinbreg() takes at most",
      "two: `y ~ count covariates | zero covariates`."), len[2L]),
      "nbzinbreg_error_input")
  }
  if (len[2L] == 1L) {
    ff <- Formula::as.Formula(stats::formula(ff),
                              stats::formula(ff, lhs = 0L))
  }
  ff
}

# Sum of the offset() terms of one part (NULL if it has none), taken from the
# model frame, whose columns are named like the deparsed terms.
terms_offset <- function(mt, mf) {
  idx <- attr(mt, "offset")
  if (is.null(idx)) return(NULL)
  vars <- as.list(attr(mt, "variables"))[-1L]
  out <- NULL
  for (i in idx) {
    v <- vars[[i]]
    nm <- paste(deparse(v, width.cutoff = 500L,
                        backtick = !is.symbol(v) && is.language(v)),
                collapse = " ")
    value <- mf[[nm]]
    if (is.null(value)) {
      abort_nbz(sprintf("The offset term `%s` could not be found.", nm),
                "nbzinbreg_error_input")
    }
    if (!is.numeric(value)) {
      abort_nbz(sprintf("The offset must be numeric, not %s.",
                        describe_value(value)), "nbzinbreg_error_input")
    }
    value <- as.vector(value)
    out <- if (is.null(out)) value else out + value
  }
  out
}

# Infinite covariate values and collinear columns, per part.
check_part_matrix <- function(M, pos, part) {
  bad_cols <- colnames(M)[colSums(!is.finite(M)) > 0]
  if (length(bad_cols)) {
    abort_nbz(sprintf(paste(
      "Covariate column(s) %s of the %s part contain infinite values. Remove",
      "or transform them before fitting."),
      paste0("`", bad_cols, "`", collapse = ", "), part),
      "nbzinbreg_error_input")
  }
  aliased <- aliased_columns(M[pos, , drop = FALSE])
  if (length(aliased)) {
    abort_nbz(sprintf(paste(
      "The %s-part model matrix is rank-deficient: %s %s a linear",
      "combination of other columns of that part. Remove %s from the %s",
      "part of the formula."), part,
      paste0("`", aliased, "`", collapse = ", "),
      if (length(aliased) == 1L) "is" else "are",
      if (length(aliased) == 1L) "it" else "them", part),
      "nbzinbreg_error_rank")
  }
  invisible(M)
}

# Check user starting values against the model matrix of one part. Unnamed
# values must have the right length; named values (with or without the
# count_/zero_ prefix) must name exactly the columns, in any order.
match_start <- function(values, cols, part) {
  if (length(values) != length(cols)) {
    abort_nbz(sprintf(paste(
      "`start$%s` has %d value(s), but the %s part has %d coefficient(s): %s."),
      part, length(values), part, length(cols),
      paste0("`", cols, "`", collapse = ", ")), "nbzinbreg_error_input")
  }
  nms <- names(values)
  if (is.null(nms)) return(stats::setNames(values, cols))
  nms <- sub(paste0("^", part, "_"), "", nms)
  if (anyNA(nms) || !setequal(nms, cols) || anyDuplicated(nms)) {
    abort_nbz(sprintf(paste(
      "The names of `start$%s` (%s) do not match the %s-part coefficients",
      "(%s). Give the values in the order of the coefficients without",
      "names, or name every coefficient once."), part,
      paste0("`", names(values), "`", collapse = ", "), part,
      paste0("`", cols, "`", collapse = ", ")), "nbzinbreg_error_input")
  }
  stats::setNames(unname(values)[match(cols, nms)], cols)
}

# Covariance matrix from the analytic or numeric Hessian at the estimates,
# after validate_hessian(). Returns the matrix (NA for excluded parameters,
# or everywhere if it cannot be computed), the validation result, whether
# the inversion failed, and with hessian_check the relative difference
# between the two Hessians.
zinb_standard_errors <- function(em, X, Z, y, offset, offset_zero, w,
                                 control) {
  analytic <- em$gradient$hessian
  keep <- !(colnames(analytic) %in% em$gradient$excluded)
  numeric_h <- function() {
    zinb_numeric_hessian(X, Z, y, offset, offset_zero, w, em$beta, em$gamma,
                         em$theta)
  }
  H <- if (control$hessian == "analytic") analytic else numeric_h()
  difference <- NA_real_
  if (control$hessian_check) {
    other <- if (control$hessian == "analytic") numeric_h() else analytic
    difference <- if (all(is.finite(H)) && all(is.finite(other))) {
      hessian_difference(H[keep, keep, drop = FALSE],
                         other[keep, keep, drop = FALSE])
    } else {
      Inf
    }
  }
  val <- validate_hessian(H, keep)
  vc <- if (val$finite) {
    zinb_vcov(val$H, keep)
  } else {
    list(vcov = matrix(NA_real_, nrow(H), ncol(H), dimnames = dimnames(H)),
         failed = TRUE)
  }
  list(vcov = vc$vcov, failed = vc$failed, validation = val,
       difference = difference)
}

# Turn the EM flags and the standard-error checks into classed warnings, and
# decide whether the fit counts as converged.
report_zinbreg_fit <- function(em, se, X, Z, control) {
  nbc <- control$nb_control
  g <- em$gradient
  fin <- em$final
  reasons <- character()
  if (identical(em$stop_reason, "maxit")) {
    reasons <- c(reasons, sprintf(paste(
      "EM reached `maxit` = %d iterations before its convergence criteria",
      "were met"), control$maxit))
  }
  if (!g$available) {
    reasons <- c(reasons, paste(
      "the gradient check could not be carried out because the negative",
      "Hessian is not positive definite, so convergence to a maximum cannot",
      "be confirmed"))
  } else if (!g$ok) {
    reasons <- c(reasons, sprintf(paste(
      "the gradient check failed: the Newton decrement %s exceeds `grad_tol`",
      "= %s (largest gradient component, in standard errors: %s)"),
      format(signif(g$decrement, 3)), format(control$grad_tol),
      format(signif(g$max_scaled, 3))))
  }
  if (!fin$zero_converged) {
    reasons <- c(reasons, sprintf(
      "the zero-part M-step of the final EM iteration did not converge (%s)",
      if (fin$zero_halving_failed) "a step could not be accepted by halving"
      else sprintf("it reached `zero_maxit` = %d rounds", control$zero_maxit)))
  }
  if (!fin$count_converged) {
    why <- c(
      if (fin$count_rank_drop) {
        "its weighted model matrix became rank-deficient"
      },
      if (fin$count_issues[["maxit"]]) {
        sprintf("it reached `maxit` = %d passes of `nb_control`", nbc$maxit)
      },
      if (fin$count_issues[["irls_maxit"]]) {
        sprintf("its IRLS stage reached `irls_maxit` = %d", nbc$irls_maxit)
      },
      if (fin$count_issues[["theta_maxit"]]) {
        sprintf("its theta stage reached `theta_maxit` = %d",
                nbc$theta_maxit)
      },
      if (fin$count_issues[["halving"]]) {
        "a step could not be accepted by halving"
      })
    if (!length(why)) why <- "it did not meet the tolerances of `nb_control`"
    reasons <- c(reasons, sprintf(
      "the count-part M-step of the final EM iteration did not converge (%s)",
      paste(why, collapse = "; ")))
  }
  v <- em$violations
  if (nrow(v)) {
    reasons <- c(reasons,
                 "EM monotonicity was violated (see the separate warning)")
  }
  # The engine's flag already requires every condition (in particular an
  # available, passed gradient check); a fit is never reported as converged
  # while any reason above applies.
  converged <- isTRUE(em$converged) && length(reasons) == 0L
  if (!converged && !length(reasons)) {
    reasons <- "the convergence conditions were not all met"
  }
  if (!converged) {
    warn_nbz(sprintf(paste(
      "zinbreg() did not converge: %s. The estimates may be unreliable; see",
      "`?zinbreg_control` for the settings."),
      paste(reasons, collapse = "; ")), "nbzinbreg_warning_convergence")
  }

  if (nrow(v)) {
    shown <- utils::head(v, 3L)
    warn_nbz(sprintf(paste(
      "EM monotonicity was violated %d time(s) in %d iteration(s), for",
      "example %s. EM and its M-steps must never lower the log-likelihood or",
      "their part of the EM objective (beyond `monotone_tol`), so this points",
      "to a numerical problem; the fit is reported as not converged."),
      nrow(v), length(unique(v$iteration)),
      paste(sprintf("iteration %d (%s fell by %s)", shown$iteration,
                    shown$check, format(signif(shown$decrease, 3))),
            collapse = ", ")), "nbzinbreg_warning_monotone")
  }

  if (em$theta_bound == "upper") {
    warn_nbz(sprintf(paste(
      "theta reached its upper bound `theta_max` = %s: the count part shows",
      "no overdispersion, so a zero-inflated Poisson model would be adequate.",
      "The standard error of theta is not available (NA)."),
      format(nbc$theta_max)), "nbzinbreg_warning_theta")
  } else if (em$theta_bound == "lower") {
    warn_nbz(sprintf(paste(
      "theta reached its lower bound `theta_min` = %s, which indicates",
      "extreme overdispersion; the fit may be unreliable. The standard error",
      "of theta is not available (NA)."), format(nbc$theta_min)),
      "nbzinbreg_warning_theta")
  } else if (em$theta > nbc$theta_warn) {
    warn_nbz(sprintf(paste(
      "The estimate of theta (%s) exceeds `theta_warn` = %s: the count part",
      "shows little or no overdispersion, and a zero-inflated Poisson model",
      "may be adequate."), format(signif(em$theta, 6)),
      format(nbc$theta_warn)), "nbzinbreg_warning_theta")
  }

  report_zinb_divergence(em, Z)
  report_zinb_se(se, control)
  converged
}

# Component-specific divergence warnings.
report_zinb_divergence <- function(em, Z) {
  st <- em$structural
  fr <- em$frozen
  quoted <- function(prefix, x) paste0("`", prefix, x, "`", collapse = ", ")
  if (length(st$zero_to_zero)) {
    warn_nbz(sprintf(paste(
      "Zero-part coefficient(s) %s run off to minus infinity: every count in",
      "the corresponding group is positive, so its probability of a",
      "structural zero goes to 0. Their estimates and p-values are not",
      "meaningful and their standard errors are NA; the other estimates are",
      "unaffected."), quoted("zero_", st$zero_to_zero)),
      "nbzinbreg_warning_divergence")
  }
  if (length(st$count) || length(st$zero_to_one)) {
    named <- c(if (length(st$count)) quoted("count_", st$count),
               if (length(st$zero_to_one)) quoted("zero_", st$zero_to_one))
    warn_nbz(sprintf(paste(
      "Every count is 0 in the group(s) of %s. Either the count part (mean",
      "going to 0) or the zero part (probability of a structural zero going",
      "to 1) can explain such a group, so these coefficient(s) have no",
      "finite estimate and are not identified. Their estimates and p-values",
      "are not meaningful and their standard errors are NA; consider merging",
      "or dropping the group."), paste(named, collapse = ", ")),
      "nbzinbreg_warning_divergence")
  }
  other_count <- setdiff(fr$count, st$count)
  if (length(other_count)) {
    warn_nbz(sprintf(paste(
      "Count-part coefficient(s) %s diverge toward minus infinity: the",
      "fitted means of some zero counts go to 0. Their estimates and p-values",
      "are not meaningful and their standard errors are NA."),
      quoted("count_", other_count)), "nbzinbreg_warning_divergence")
  }
  if (em$zero_collapsed) {
    warn_nbz(sprintf(paste(
      "Every fitted probability of a structural zero is below 1e-8: the data",
      "show no zero inflation beyond what the negative binomial already",
      "predicts, and `nbreg()` may be adequate. The zero-part coefficients",
      "(%s) are not meaningful and their standard errors are NA."),
      quoted("zero_", colnames(Z))), "nbzinbreg_warning_divergence")
  } else if (isTRUE(em$zero_separated)) {
    warn_nbz(sprintf(paste(
      "The zero part is completely separated: the zero-part covariates split",
      "the observations exactly into structural zeros and non-structural",
      "observations (for example, every count above some value of a",
      "covariate is 0 and every count below it is positive). Every fitted",
      "probability of a structural zero is within 1e-8 of 0 or 1, and",
      "scaling the zero-part coefficients to infinity stays inside the 95%%",
      "likelihood region, so the data put no bound on their size. The",
      "zero-part coefficients (%s) have no finite estimates; they are not",
      "meaningful and their standard errors are NA. Consider a simpler zero",
      "part."), quoted("zero_", em$gradient$separated)),
      "nbzinbreg_warning_divergence")
  } else {
    other_zero <- setdiff(fr$zero, st$zero)
    if (length(other_zero)) {
      warn_nbz(sprintf(paste(
        "Zero-part coefficient(s) %s diverge: the fitted probabilities of a",
        "structural zero go to 0 or 1 for the observations they affect.",
        "Their estimates and p-values are not meaningful and their standard",
        "errors are NA."), quoted("zero_", other_zero)),
        "nbzinbreg_warning_divergence")
    }
  }
  invisible()
}

# Standard-error warnings.
report_zinb_se <- function(se, control) {
  val <- se$validation
  if (!val$finite) {
    warn_nbz(sprintf(paste(
      "The %s Hessian contains non-finite values, so the standard errors are",
      "not available (NA)."), control$hessian), "nbzinbreg_warning_se")
    return(invisible())
  }
  if (!val$symmetric) {
    warn_nbz(sprintf(paste(
      "The %s Hessian is not symmetric to within 1e-8 of its largest entry;",
      "it was made symmetric before inversion, but the standard errors may",
      "be inaccurate."), control$hessian), "nbzinbreg_warning_se")
  }
  if (se$failed) {
    warn_nbz(paste(
      "The negative Hessian is not positive definite, so the covariance",
      "matrix cannot be computed; the standard errors are NA. The fit may not",
      "be at a maximum."), "nbzinbreg_warning_se")
  } else if (is.finite(val$kappa) && val$kappa > 1e12) {
    warn_nbz(sprintf(paste(
      "The Hessian is ill-conditioned (condition number %s), so the standard",
      "errors may be inaccurate."), format(signif(val$kappa, 3))),
      "nbzinbreg_warning_se")
  }
  if (isTRUE(control$hessian_check) && !(se$difference <= 1e-4)) {
    warn_nbz(sprintf(paste(
      "The analytic and numeric Hessians differ by a relative %s (more than",
      "1e-4), so the standard errors may be inaccurate."),
      format(signif(se$difference, 3))), "nbzinbreg_warning_se")
  }
  invisible()
}

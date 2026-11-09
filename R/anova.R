# Internals for the likelihood-ratio tests of nested fits (anova()) and the
# overdispersion test. The exported functions are in R/inference.R.

# Do the columns of S lie in the column space of L? Each column of S is
# projected onto span(L) (QR of L, rows with positive weight); its residual
# must be below 1e-8 times the column's own length.
columns_nested <- function(S, L, pos) {
  S <- S[pos, , drop = FALSE]
  L <- L[pos, , drop = FALSE]
  if (ncol(S) == 0L) return(TRUE)
  res <- qr.resid(qr(L), S)
  lens <- sqrt(colSums(S^2))
  all(sqrt(colSums(res^2)) <= 1e-8 * pmax(lens, .Machine$double.xmin))
}

# Coefficients of L that reproduce the linear predictor X_S b_S exactly
# (possible because span(X_S) is inside span(X_L)).
embed_coefficients <- function(S, b, L, pos) {
  target <- drop(S[pos, , drop = FALSE] %*% b)
  out <- qr.coef(qr(L[pos, , drop = FALSE]), target)
  out[is.na(out)] <- 0
  stats::setNames(out, colnames(L))
}

# Why a fit is not regular for an ordinary chi-squared reference (an empty
# vector if it is): it did not converge, has coefficients without a finite
# estimate, has theta on a bound, or (zinbreg) its zero part has collapsed.
irregular_reasons <- function(fit) {
  pos <- fit$prior.weights > 0
  c(if (!isTRUE(fit$converged)) "did not converge",
    if (length(fit$diverging)) "has coefficients without a finite estimate",
    if (!identical(fit$theta_bound, "none")) "has theta on a bound",
    if (inherits(fit, "zinbreg") && all(fit$pi[pos] < 1e-8)) {
      "has a collapsed zero part (no zero inflation)"
    })
}

# Model label: the right-hand side of the formula.
model_label <- function(fit) {
  f <- if (inherits(fit, "zinbreg")) stats::formula(fit$formula) else fit$formula
  paste(deparse(f[[3L]], width.cutoff = 500L), collapse = " ")
}

# Refit the larger model L from the embedding of the smaller fit S. Both
# engines never lower the log-likelihood, so the refit is at least l_S; if it
# beats the stored l_L, the stored larger fit was not at its maximum.
# Returns the refit's log-likelihood (NA if the refit fails).
embedding_refit_loglik <- function(S, L) {
  pos <- L$prior.weights > 0
  y <- unname(L$y)
  if (inherits(L, "zinbreg")) {
    b <- embed_coefficients(S$x, S$coefficients$count, L$x, pos)
    g <- embed_coefficients(S$z, S$coefficients$zero, L$z, pos)
    ctl <- L$control
    ctl$trace <- FALSE
    ctl$nb_control$trace <- FALSE
    ctl$start <- list(count = b, zero = g, theta = S$theta)
    fit <- tryCatch(zinb_em(L$x, L$z, y, L$offset$count, L$offset$zero,
                            L$prior.weights, ctl),
                    nbzinbreg_error = function(e) NULL)
  } else {
    b <- embed_coefficients(S$x, S$coefficients, L$x, pos)
    ctl <- L$control
    ctl$trace <- FALSE
    fit <- tryCatch(nbreg_fit(L$x, y, L$offset, L$prior.weights, ctl,
                              start = list(beta = b, theta = S$theta),
                              compute_se = FALSE),
                    nbzinbreg_error = function(e) NULL)
  }
  if (is.null(fit)) NA_real_ else fit$loglik
}

# The table of likelihood-ratio tests for a list of fits of class `cls`.
anova_fits <- function(fits, cls) {
  other <- setdiff(c("nbreg", "zinbreg"), cls)
  if (length(fits) < 2L) {
    abort_nbz(sprintf(paste(
      "anova() for %s fits needs at least two fits of nested models, for",
      "example anova(fit_small, fit_large). For one fit, use summary()."),
      cls), "nbzinbreg_error_input")
  }
  for (i in seq_along(fits)) {
    f <- fits[[i]]
    if (inherits(f, other)) {
      abort_nbz(sprintf(paste(
        "anova() compares fits of the same model class; fit %d is a `%s`",
        "fit. An nbreg fit is the boundary pi = 0 of a zinbreg fit, which a",
        "chi-squared likelihood-ratio test does not cover: use vuong_test()."),
        i, other), "nbzinbreg_error_input")
    }
    if (!inherits(f, cls)) {
      abort_nbz(sprintf("Argument %d of anova() must be a `%s` fit, not %s.",
                        i, cls, describe_value(f)), "nbzinbreg_error_input")
    }
  }
  first <- fits[[1L]]
  same <- function(a, b) isTRUE(all.equal(unname(a), unname(b), tolerance = 0))
  for (i in seq_along(fits)[-1L]) {
    f <- fits[[i]]
    problem <- if (!identical(names(f$y), names(first$y)) ||
                   !same(f$y, first$y)) {
      "the response (or the rows used)"
    } else if (!same(f$prior.weights, first$prior.weights)) {
      "the weights"
    } else if (cls == "nbreg" && !same(f$offset, first$offset)) {
      "the offset"
    } else if (cls == "zinbreg" &&
               (!same(f$offset$count, first$offset$count) ||
                !same(f$offset$zero, first$offset$zero))) {
      "the offsets"
    }
    if (!is.null(problem)) {
      abort_nbz(sprintf(paste(
        "Fits 1 and %d differ in %s. A likelihood-ratio test needs every",
        "fit made on exactly the same data (same rows, response, weights and",
        "offsets); check `subset`, missing values and `na.action`."),
        i, problem), "nbzinbreg_error_input")
    }
  }
  df <- vapply(fits, function(f) as.numeric(f$df), numeric(1))
  if (anyDuplicated(df)) {
    abort_nbz(sprintf(paste(
      "Fits %s have the same number of parameters (%s), so they cannot be",
      "nested in each other. anova() needs strictly nested models."),
      paste(which(df %in% df[duplicated(df)]), collapse = " and "),
      format(df[duplicated(df)][1L])), "nbzinbreg_error_input")
  }
  o <- order(df)
  fits <- fits[o]
  pos <- first$prior.weights > 0
  k <- length(fits)
  lr <- p <- rep(NA_real_, k)
  ddf <- rep(NA_integer_, k)
  test <- note <- rep("", k)
  for (j in seq_len(k)[-1L]) {
    S <- fits[[j - 1L]]
    L <- fits[[j]]
    nested <- columns_nested(S$x, L$x, pos) &&
      (cls == "nbreg" || columns_nested(S$z, L$z, pos))
    if (!nested) {
      abort_nbz(sprintf(paste(
        "The models are not nested: the %s of model %d (%s) is not contained",
        "in that of model %d (%s). Every smaller model must be a special",
        "case of the next larger one."),
        if (cls == "zinbreg") "count or zero part" else "model matrix",
        o[j - 1L], model_label(S), o[j], model_label(L)),
        "nbzinbreg_error_input")
    }
    tol <- 1e-6 * (abs(L$loglik) + 1)
    ll_ref <- embedding_refit_loglik(S, L)
    if (is.finite(ll_ref) && ll_ref > L$loglik + tol) {
      abort_nbz(sprintf(paste(
        "The larger fit %d (%s) is not at its maximum: refitting it from the",
        "smaller fit gives log-likelihood %s, above its stored %s. Refit it",
        "(for example with a smaller `tol`) before testing."), o[j],
        model_label(L), format(ll_ref, digits = 10),
        format(L$loglik, digits = 10)), "nbzinbreg_error_fit")
    }
    stat <- 2 * (L$loglik - S$loglik)
    if (stat < -2 * tol) {
      abort_nbz(sprintf(paste(
        "The larger fit %d has a lower log-likelihood than the smaller fit",
        "%d (%s against %s), which is impossible for nested models at their",
        "maxima."), o[j], o[j - 1L], format(L$loglik, digits = 10),
        format(S$loglik, digits = 10)), "nbzinbreg_error_fit")
    }
    lr[j] <- max(stat, 0)
    ddf[j] <- as.integer(df[o[j]] - df[o[j - 1L]])
    test[j] <- sprintf("%d vs %d", j - 1L, j)
    reasons <- c(
      if (!is.finite(ll_ref)) "the verification refit failed",
      if (length(r <- irregular_reasons(S))) {
        sprintf("model %d %s", j - 1L, paste(r, collapse = ", "))
      },
      if (length(r <- irregular_reasons(L))) {
        sprintf("model %d %s", j, paste(r, collapse = ", "))
      })
    if (length(reasons)) {
      note[j] <- paste0(paste(reasons, collapse = "; "),
                        ": chi-squared reference not valid")
      cls_w <- if (any(grepl("theta on a bound", reasons))) {
        "nbzinbreg_warning_theta"
      } else if (any(grepl("finite estimate|collapsed", reasons))) {
        "nbzinbreg_warning_divergence"
      } else {
        "nbzinbreg_warning_convergence"
      }
      warn_nbz(sprintf(paste(
        "No p-value for test %s: %s, so the likelihood-ratio statistic does",
        "not follow a chi-squared distribution reliably."), test[j],
        paste(reasons, collapse = "; ")), cls_w)
    } else {
      p[j] <- stats::pchisq(lr[j], ddf[j], lower.tail = FALSE)
    }
  }
  resp <- paste(deparse(stats::formula(first$formula)[[2L]]), collapse = " ")
  out <- data.frame(
    Model = vapply(fits, model_label, character(1)),
    theta = vapply(fits, function(f) f$theta, numeric(1)),
    `Resid. df` = vapply(fits, function(f) as.integer(f$df.residual),
                         integer(1)),
    `2 x log-lik.` = vapply(fits, function(f) 2 * f$loglik, numeric(1)),
    Test = test, df = ddf, `LR stat.` = lr, `Pr(Chi)` = p, Note = note,
    check.names = FALSE, stringsAsFactors = FALSE)
  structure(out, heading = c(sprintf(
    "Likelihood-ratio tests of %s models\n",
    if (cls == "zinbreg") "zero-inflated negative binomial" else
      "negative binomial"), sprintf("Response: %s", resp)),
    class = c("nbzinbreg_anova", "anova", "data.frame"))
}

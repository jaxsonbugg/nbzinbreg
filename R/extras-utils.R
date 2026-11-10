# Internal helpers shared by the extras (confint, anova, overdispersion and
# Vuong tests, simulate, bootstraps): argument validators, reproducible
# random numbers, simulation from a fitted model, the frequency-weight
# bootstrap sample, refits with a single status, and a Poisson fit.
# None of them is exported.

# ---------------------------------------------------------------------------
# Argument validators
# ---------------------------------------------------------------------------

# A probability strictly between 0 and 1 (confidence level or test level).
check_probability <- function(x, name, example) {
  check_number(x, name)
  if (x <= 0 || x >= 1) {
    abort_nbz(sprintf(
      "`%s` must be a number strictly between 0 and 1 (for example %s), not %s.",
      name, example, describe_value(x)), "nbzinbreg_error_input")
  }
  invisible(x)
}

check_level <- function(x, name = "level") check_probability(x, name, "0.95")

check_alpha <- function(x, name = "alpha") check_probability(x, name, "0.05")

# A seed for set.seed(): NULL, or a single finite whole number that
# as.integer() can hold. set.seed() coerces its argument with as.integer(),
# so 1.5 would be silently truncated and 2^31 would become NA.
check_seed <- function(seed) {
  if (is.null(seed)) return(invisible(NULL))
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed)) {
    abort_nbz(sprintf(
      "`seed` must be NULL or a single finite whole number, not %s.",
      describe_value(seed)), "nbzinbreg_error_input")
  }
  if (abs(seed - round(seed)) > 0) {
    abort_nbz(sprintf(paste(
      "`seed` must be a whole number (set.seed() would silently truncate",
      "%s)."), describe_value(seed)), "nbzinbreg_error_input")
  }
  if (abs(seed) > .Machine$integer.max) {
    abort_nbz(sprintf(paste(
      "`seed` must be between -%d and %d (set.seed() needs an integer), not",
      "%s."), .Machine$integer.max, .Machine$integer.max,
      describe_value(seed)), "nbzinbreg_error_input")
  }
  invisible(seed)
}

# A whole number in [min, max] (numbers of replicates or simulations).
check_count_arg <- function(x, name, min = 1, max = .Machine$integer.max) {
  check_whole(x, name, min = min)
  if (x > max) {
    abort_nbz(sprintf("`%s` must be at most %s, not %s.", name,
                      format(max, big.mark = ",", scientific = FALSE),
                      describe_value(x)), "nbzinbreg_error_input")
  }
  invisible(x)
}

# Extra arguments in `...` are an error, so a misspelled argument is never
# silently ignored. `fun` names the function, `args` its arguments.
check_no_dots <- function(n_dots, dot_names, fun, args) {
  if (n_dots == 0L) return(invisible())
  if (is.null(dot_names)) dot_names <- rep("", n_dots)
  dot_names[is.na(dot_names)] <- ""
  named <- dot_names[dot_names != ""]
  what <- if (length(named)) {
    sprintf("Unused argument(s) in %s(): %s.", fun,
            paste0("`", named, "`", collapse = ", "))
  } else {
    sprintf("%d unnamed extra argument(s) in %s().", n_dots, fun)
  }
  abort_nbz(sprintf("%s The arguments are %s.", what,
                    paste0("`", args, "`", collapse = ", ")),
            "nbzinbreg_error_input")
}

# Which parameters: missing means all; otherwise distinct names from
# `names`, or distinct whole-number indices into them. Returns the names.
check_parm <- function(parm, names, missing_parm) {
  if (missing_parm) return(names)
  valid <- paste0("`", names, "`", collapse = ", ")
  if (is.character(parm)) {
    if (length(parm) == 0L || anyNA(parm)) {
      abort_nbz(sprintf(paste(
        "`parm` must name at least one parameter, without missing values.",
        "The parameters are %s."), valid), "nbzinbreg_error_input")
    }
    unknown <- setdiff(parm, names)
    if (length(unknown)) {
      abort_nbz(sprintf("Unknown parameter(s) in `parm`: %s. The parameters are %s.",
                        paste0("`", unknown, "`", collapse = ", "), valid),
                "nbzinbreg_error_input")
    }
  } else if (is.numeric(parm) && !is.factor(parm)) {
    # Finiteness is checked first: Inf - round(Inf) is NaN, and the later
    # comparisons would then give NA instead of TRUE/FALSE.
    if (length(parm) == 0L || anyNA(parm) || any(!is.finite(parm)) ||
        any(abs(parm - round(parm)) > 0) || any(parm < 1) ||
        any(parm > length(names))) {
      abort_nbz(sprintf(paste(
        "`parm` given as numbers must be whole numbers from 1 to %d (the",
        "parameters are %s), not %s."), length(names), valid,
        describe_value(parm)), "nbzinbreg_error_input")
    }
    parm <- names[parm]
  } else {
    abort_nbz(sprintf(paste(
      "`parm` must be parameter names or indices, not %s. The parameters",
      "are %s."), describe_value(parm), valid), "nbzinbreg_error_input")
  }
  if (anyDuplicated(parm)) {
    abort_nbz(sprintf("`parm` names %s more than once.",
                      paste0("`", unique(parm[duplicated(parm)]), "`",
                             collapse = ", ")), "nbzinbreg_error_input")
  }
  parm
}

# The fit must be an nbreg or zinbreg object (or one of `classes`).
check_fit_class <- function(object, classes, arg = "object") {
  if (!inherits(object, classes)) {
    abort_nbz(sprintf("`%s` must be %s fit, not %s.", arg,
                      paste0("a `", classes, "`", collapse = " or "),
                      if (is.object(object)) {
                        sprintf("an object of class \"%s\"",
                                paste(class(object), collapse = "\", \""))
                      } else {
                        describe_value(object)
                      }), "nbzinbreg_error_input")
  }
  invisible(object)
}

# ---------------------------------------------------------------------------
# Reproducible random numbers
# ---------------------------------------------------------------------------

# Evaluate `expr` with the random-number conventions of stats::simulate():
#   * if the RNG has not been used yet (no .Random.seed in the global
#     environment), it is initialized with one runif() draw;
#   * seed = NULL: use the current stream, and record its state;
#   * otherwise: save the current state, set.seed(seed), and restore the
#     saved state on exit, so the caller's stream is not disturbed.
# The value gets a "seed" attribute, as in stats::simulate().
with_seed <- function(seed, expr) {
  check_seed(seed)
  genv <- globalenv()
  if (!exists(".Random.seed", envir = genv, inherits = FALSE)) stats::runif(1)
  if (is.null(seed)) {
    rng_state <- get(".Random.seed", envir = genv, inherits = FALSE)
  } else {
    saved <- get(".Random.seed", envir = genv, inherits = FALSE)
    set.seed(seed)
    rng_state <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", saved, envir = genv))
  }
  value <- expr
  attr(value, "seed") <- rng_state
  value
}

# ---------------------------------------------------------------------------
# Simulation from a fitted model
# ---------------------------------------------------------------------------

# Fitted count-part means of a fit (nbreg: fitted values; zinbreg: mu).
fit_mu <- function(object) {
  unname(if (inherits(object, "zinbreg")) object$mu else object$fitted.values)
}

# counts[i] independent responses for row i of the fit: negative binomial
# with the fitted mean and theta, and for zinbreg a structural zero first
# with probability pi_i. Returns the row index of every draw and the draws.
# This is the only place in the package that calls rnbinom()/rbinom().
simulate_rows <- function(object, counts) {
  idx <- rep.int(seq_along(counts), counts)
  y <- stats::rnbinom(length(idx), size = object$theta,
                      mu = fit_mu(object)[idx])
  if (inherits(object, "zinbreg")) {
    structural <- stats::rbinom(length(idx), 1L, unname(object$pi)[idx]) == 1L
    y[structural] <- 0
  }
  list(row = idx, y = y)
}

# Frequency weights for a parametric bootstrap: whole numbers. Row i stands
# for w_i independent observations. Weights that are not whole numbers have
# no such meaning, so a bootstrap is refused for them.
#   * Rounding noise is tolerated: a positive weight within 1e-8 of a
#     positive whole number (for example 2 + 1e-10) is rounded to it.
#   * A positive weight that would round to 0 (for example 1e-9) is refused:
#     rounding it would silently remove an observation from the fit.
#   * Weights of exactly 0 stay 0 (those rows are not in the fit).
check_frequency_weights <- function(w, what) {
  pos <- w > 0
  frac <- pos & abs(w - round(w)) > 1e-8
  vanishing <- pos & !frac & round(w) < 1
  if (any(frac)) {
    abort_nbz(sprintf(paste(
      "%s needs frequency weights (whole numbers): each row then stands for",
      "that many independent observations, which are simulated separately.",
      "This fit has %d non-integer weight(s), for example %s, which have no",
      "such meaning."), what, sum(frac), format_values(w[frac])),
      "nbzinbreg_error_input")
  }
  if (any(vanishing)) {
    abort_nbz(sprintf(paste(
      "%s needs frequency weights (whole numbers). This fit has %d positive",
      "weight(s) smaller than 1, for example %s, which would round to 0 and",
      "silently drop those observations. Use whole-number weights, or weight",
      "0 to leave a row out."), what, sum(vanishing),
      format_values(w[vanishing])), "nbzinbreg_error_input")
  }
  round(w)
}

# Upper limit on the number of individual draws in one bootstrap data set,
# sum(w): every replicate simulates that many responses and sorts them, which
# needs about 40 bytes of working memory per draw (1e7 draws is roughly
# 400 MB). The sum is formed in double precision, so it cannot overflow
# (huge weights give a huge or infinite sum, which fails the check), and it
# is checked before any vector of that length is allocated.
max_bootstrap_draws <- 1e7

check_bootstrap_size <- function(w_int, what) {
  total <- sum(as.numeric(w_int))
  if (!is.finite(total) || total > max_bootstrap_draws) {
    abort_nbz(sprintf(paste(
      "%s would simulate %s individual observations per replicate (the sum",
      "of the frequency weights), more than the limit of %s (about 400 MB of",
      "working memory). Fit the model to fewer or smaller-weighted rows."),
      what, format(total, big.mark = ",", scientific = FALSE),
      format(max_bootstrap_draws, big.mark = ",", scientific = FALSE)),
      "nbzinbreg_error_input")
  }
  invisible(total)
}

# One bootstrap data set under frequency weights: w_i independent draws for
# row i (none for weight 0), aggregated to distinct (row, value) pairs with
# their multiplicities as weights. Refitting this is exactly refitting the
# fully expanded data, without forming it. (Drawing once per row and keeping
# weight w_i would treat the w_i observations as identical copies.)
# Callers (boot_se() and the bootstrap LR test) run check_bootstrap_size()
# once before their replicate loop; it is repeated here as a safeguard, so
# no oversized allocation can happen even if a caller forgets.
bootstrap_sample <- function(object, w_int) {
  check_bootstrap_size(w_int, "The bootstrap")
  d <- simulate_rows(object, w_int)
  o <- order(d$row, d$y)
  r <- d$row[o]
  y <- d$y[o]
  n <- length(r)
  first <- c(TRUE, r[-1L] != r[-n] | y[-1L] != y[-n])
  list(row = r[first], y = y[first], w = tabulate(cumsum(first)))
}

# ---------------------------------------------------------------------------
# Refits with a single status
# ---------------------------------------------------------------------------

# Refit the model of `object` to a new response: the same model matrices,
# offsets and control (rows taken by `rows`), weights `w`, optionally from
# `start` (nbreg: list(beta, theta); zinbreg: list(count, zero, theta)).
# The engines never warn; any other warning is captured and counted.
# Progress messages are switched off.
#
# Status: the first that applies, in this order.
#   error          an nbzinbreg error, or a response the model cannot fit
#                  (all zero; for zinbreg also no zero)
#   not_converged  the engine's converged flag is FALSE (for zinbreg this
#                  covers maxit, the gradient check, the final M-steps and
#                  monotonicity)
#   zero_boundary  zinbreg only: every pi < 1e-8 (no zero inflation left)
#   diverging      coefficients without a finite estimate (held fixed by
#                  the engine, or structural)
#   theta_bound    theta ended on theta_min or theta_max
#   ok             otherwise
# `maxit` (zinbreg only) overrides the EM iteration limit of the fit's
# control; NULL keeps it.
refit_model <- function(object, y, w, rows = seq_along(y), start = NULL,
                        maxit = NULL) {
  zinb <- inherits(object, "zinbreg")
  pos <- w > 0
  out <- list(status = "error", message = NA_character_, loglik = NA_real_,
              coefficients = NULL, theta = NA_real_, diverging = character(),
              warnings = 0L)
  if (all(y[pos] == 0)) {
    out$message <- "every simulated response with positive weight is 0"
    return(out)
  }
  if (zinb && !any(y[pos] == 0)) {
    out$message <- "no simulated response with positive weight is 0"
    return(out)
  }
  n_warn <- 0L
  muffle <- function(wcnd) {
    n_warn <<- n_warn + 1L
    invokeRestart("muffleWarning")
  }
  X <- object$x[rows, , drop = FALSE]
  fit <- tryCatch(withCallingHandlers({
    if (zinb) {
      ctrl <- object$control
      ctrl$trace <- FALSE
      ctrl$nb_control$trace <- FALSE
      ctrl$start <- start
      if (!is.null(maxit)) ctrl$maxit <- as.integer(maxit)
      zinb_em(X, object$z[rows, , drop = FALSE], y,
              object$offset$count[rows], object$offset$zero[rows], w, ctrl)
    } else {
      ctrl <- object$control
      ctrl$trace <- FALSE
      nbreg_fit(X, y, object$offset[rows], w, ctrl, start = start,
                compute_se = FALSE)
    }
  }, warning = muffle), nbzinbreg_error = function(e) e)
  out$warnings <- n_warn
  if (inherits(fit, "condition")) {
    out$message <- conditionMessage(fit)
    return(out)
  }
  if (zinb) {
    out$coefficients <- list(count = fit$beta, zero = fit$gamma)
    out$diverging <- setdiff(fit$gradient$excluded, "Log(theta)")
    out$status <- if (!isTRUE(fit$converged)) {
      "not_converged"
    } else if (isTRUE(fit$zero_collapsed)) {
      "zero_boundary"
    } else if (length(out$diverging)) {
      "diverging"
    } else if (!identical(fit$theta_bound, "none")) {
      "theta_bound"
    } else {
      "ok"
    }
  } else {
    out$coefficients <- fit$coefficients
    out$diverging <- union(fit$diverging, structural_divergence(X, y, w))
    out$status <- if (!isTRUE(fit$converged)) {
      "not_converged"
    } else if (length(out$diverging)) {
      "diverging"
    } else if (!identical(fit$theta_bound, "none")) {
      "theta_bound"
    } else {
      "ok"
    }
  }
  out$theta <- fit$theta
  out$loglik <- fit$loglik
  out$engine <- fit
  out
}

# ---------------------------------------------------------------------------
# Poisson regression (for the overdispersion test)
# ---------------------------------------------------------------------------

# The Poisson MLE by IRLS, written like the NB engine's Poisson start but run
# to convergence: one weighted least-squares step per round, with
# backtracking so the log-likelihood never falls. Converged when the relative
# change in the log-likelihood is below 1e-10 and every coefficient's
# relative change is below 1e-5. `ll_converged` records whether the
# log-likelihood alone settled (it does even when a coefficient runs off to
# minus infinity for a group of zeros, where the coefficients cannot).
poisson_fit <- function(X, y, offset, w, maxit = 100L, max_halving = 30L) {
  beta <- poisson_initial_beta(X, y, offset, w)
  if (is.null(beta)) {
    abort_nbz("The Poisson fit failed: its first least-squares step failed.",
              "nbzinbreg_error_fit")
  }
  names(beta) <- colnames(X)
  eta <- lin_pred(X, beta, offset)
  ll <- pois_loglik(y, eta, w)
  if (!is.finite(ll)) {
    abort_nbz(paste(
      "The Poisson fit failed: its first step gave fitted means that are",
      "infinite, or zero where counts are positive."), "nbzinbreg_error_fit")
  }
  objective <- function(b) pois_loglik(y, lin_pred(X, b, offset), w)
  converged <- FALSE
  ll_converged <- FALSE
  it <- 0L
  for (it in seq_len(maxit)) {
    proposal <- wls_solve(X, w * exp(eta), (eta - offset) + mu_ratio(y, eta))
    if (is.null(proposal)) break
    step <- backtrack(beta, proposal, objective, ll, max_halving)
    if (!step$ok) break
    rel_ll <- rel_change(step$ll, ll)
    rel_beta <- max(abs(step$par - beta) / (abs(beta) + 0.1))
    beta <- step$par
    eta <- lin_pred(X, beta, offset)
    ll <- step$ll
    ll_converged <- rel_ll < 1e-10
    if (ll_converged && rel_beta < 1e-5) {
      converged <- TRUE
      break
    }
  }
  list(coefficients = beta, loglik = ll, linear.predictors = eta,
       converged = converged, ll_converged = ll_converged, iter = it)
}

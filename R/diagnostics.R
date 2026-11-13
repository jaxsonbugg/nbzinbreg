#' Observed and Expected Count Frequencies
#'
#' @description
#' Compares how often each count value (0, 1, 2, ...) occurs in the data with
#' how often the fitted model expects it to occur. The expected frequency of
#' a value \eqn{k} is the sum over observations of the fitted probability
#' \eqn{P(Y_i = k)}. The row for \eqn{k = 0} shows whether the model predicts
#' the right number of zeros.
#'
#' @details
#' Frequencies are weighted by the prior weights, and observations with
#' weight 0 are left out. So with frequency weights, the table is that of
#' the expanded data. For the negative binomial (NB) model,
#' \eqn{P(Y_i = k)} is the NB probability with mean \eqn{\mu_i} and the
#' fitted `theta`. For the zero-inflated model it is
#' \eqn{\pi_i I(k = 0) + (1 - \pi_i) P_{NB}(k)}.
#'
#' Counts above `max_count` are grouped in a last row labelled
#' `"<max_count + 1>+"` (for example `"11+"` when `max_count = 10`). Its
#' expected frequency is \eqn{\sum_i w_i (1 - \pi_i) P_{NB}(Y_i > m)}, with
#' \eqn{m} = `max_count` and \eqn{\pi_i = 0} for NB. Both columns therefore
#' add up to the total weight.
#'
#' The table needs one probability per observation and count value, so a
#' computation of more than 1e8 probabilities (observations times
#' `max_count + 2`) is refused before it starts.
#'
#' @section Errors and warnings:
#' * `nbzinbreg_error_input`: `object` of the wrong class, `max_count` not
#'   `NULL` or a whole number of at least 0, a table of more than 1e8
#'   probabilities, or unused arguments.
#'
#' @param object An `nbreg` or `zinbreg` fit.
#' @param max_count Largest count to tabulate separately; larger counts are
#'   grouped in the last row. By default the largest observed count (with
#'   positive weight).
#' @param ... Not used (an error if supplied).
#' @return A data frame with columns `count` (character: `"0"`, `"1"`, ...,
#'   and the grouped last row), `observed` and `expected`. It has
#'   `max_count + 2` rows, and attributes `max_count` and `total` (the total
#'   weight).
#' @seealso [rootogram()], [predict.zinbreg()]
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200))
#' d$y <- ifelse(runif(200) < 0.3, 0, rnbinom(200, size = 2, mu = exp(1 + d$x)))
#' expected_counts(nbreg(y ~ x, data = d), max_count = 8)
#' expected_counts(zinbreg(y ~ x | 1, data = d), max_count = 8)
#' @export
expected_counts <- function(object, max_count = NULL, ...) {
  check_fit_class(object, c("nbreg", "zinbreg"))
  check_no_dots(...length(), ...names(), "expected_counts",
                c("object", "max_count"))
  y <- unname(object$y)
  w <- object$prior.weights
  pos <- w > 0
  if (is.null(max_count)) {
    max_count <- max(y[pos])
  } else {
    check_whole(max_count, "max_count", min = 0)
  }
  max_count <- as.integer(round(max_count))
  n_pos <- sum(pos)
  cells <- n_pos * (as.numeric(max_count) + 2)
  if (cells > max_count_cells) {
    abort_nbz(sprintf(paste(
      "expected_counts() would compute %s probabilities (%d observations",
      "times %s count values), more than the limit of 1e8. Use a smaller",
      "`max_count`."), format(cells, big.mark = ",", scientific = FALSE),
      n_pos, format(as.numeric(max_count) + 2, big.mark = ",",
                    scientific = FALSE)), "nbzinbreg_error_input")
  }
  p <- count_model_parts(object)
  mu <- p$mu[pos]
  pi <- p$pi[pos]
  ww <- w[pos]
  yy <- y[pos]
  k <- 0:max_count
  # one count value at a time, so memory stays O(n)
  expected <- vapply(k, function(j) {
    sum(ww * (1 - pi) * stats::dnbinom(j, size = object$theta, mu = mu))
  }, numeric(1))
  expected[1L] <- expected[1L] + sum(ww * pi)
  tail_expected <- sum(ww * (1 - pi) * stats::pnbinom(
    max_count, size = object$theta, mu = mu, lower.tail = FALSE))
  observed <- vapply(k, function(j) sum(ww[yy == j]), numeric(1))
  out <- data.frame(
    count = c(as.character(k), paste0(max_count + 1L, "+")),
    observed = c(observed, sum(ww[yy > max_count])),
    expected = c(expected, tail_expected),
    stringsAsFactors = FALSE)
  attr(out, "max_count") <- max_count
  attr(out, "total") <- sum(ww)
  out
}

# Largest number of probabilities expected_counts() computes.
max_count_cells <- 1e8

# The NB means and the structural-zero probabilities of a fit (pi = 0 for
# nbreg), for every row.
count_model_parts <- function(object) {
  if (inherits(object, "zinbreg")) {
    b <- object$coefficients$count
    g <- object$coefficients$zero
    b[is.na(b)] <- 0
    g[is.na(g)] <- 0
    list(mu = exp(lin_pred(object$x, b, object$offset$count)),
         pi = stats::plogis(lin_pred(object$z, g, object$offset$zero)))
  } else {
    b <- object$coefficients
    b[is.na(b)] <- 0
    mu <- exp(lin_pred(object$x, b, object$offset))
    list(mu = mu, pi = numeric(length(mu)))
  }
}

#' Rootogram of Observed and Expected Count Frequencies
#'
#' @description
#' Draws a hanging rootogram (Kleiber and Zeileis, 2016). The square roots of
#' the expected frequencies from [expected_counts()] form a curve, and a bar
#' for each observed frequency hangs from it. A bar that stops above the
#' zero line marks a count the model over-predicts. A bar that crosses below
#' it marks a count the model under-predicts.
#'
#' @details
#' The counts 0 to `max_count` are drawn. The grouped last row of
#' [expected_counts()] is left out of the plot but kept in the returned
#' table. The plot uses base graphics.
#'
#' @section Errors and warnings:
#' As for [expected_counts()]. Arguments in `...` go to [graphics::plot()],
#' which reports its own errors.
#'
#' @param object An `nbreg` or `zinbreg` fit.
#' @param max_count Largest count to show; see [expected_counts()].
#' @param ... Graphical arguments for [graphics::plot()], such as `main`,
#'   `xlab`, `ylab`, `xlim` or `ylim`. They replace the defaults.
#' @return The table from [expected_counts()], invisibly.
#' @references Kleiber, C. and Zeileis, A. (2016). Visualizing count data
#'   regressions using rootograms. *The American Statistician*, 70(3),
#'   296-303.
#' @seealso [expected_counts()]
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200))
#' d$y <- ifelse(runif(200) < 0.3, 0, rnbinom(200, size = 2, mu = exp(1 + d$x)))
#' rootogram(nbreg(y ~ x, data = d), max_count = 10)
#' rootogram(zinbreg(y ~ x | 1, data = d), max_count = 10)
#' @export
rootogram <- function(object, max_count = NULL, ...) {
  tab <- expected_counts(object, max_count = max_count)
  shown <- seq_len(nrow(tab) - 1L)
  k <- shown - 1L
  top <- sqrt(tab$expected[shown])
  bottom <- top - sqrt(tab$observed[shown])
  defaults <- list(
    xlab = "Count", ylab = "sqrt(Frequency)",
    main = sprintf("Hanging rootogram (%s)", class(object)[1L]),
    xlim = c(-0.5, max(k) + 0.5), ylim = range(0, bottom, top))
  dots <- list(...)
  args <- c(dots, defaults[setdiff(names(defaults), names(dots))])
  do.call(graphics::plot, c(list(x = NA_real_, y = NA_real_, type = "n"),
                            args))
  graphics::rect(k - 0.4, bottom, k + 0.4, top, col = "grey80",
                 border = "grey40")
  graphics::abline(h = 0, lty = 2)
  graphics::lines(k, top, type = "b", pch = 19, col = "firebrick")
  invisible(tab)
}

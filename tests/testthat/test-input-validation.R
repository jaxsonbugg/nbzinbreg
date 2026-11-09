# Input validation: every documented error is triggered here, and each test
# checks both the condition class and a key phrase of the message.

expect_input_error <- function(expr, pattern) {
  expect_error(expr, class = "nbzinbreg_error_input", regexp = pattern)
}

# ---- nbreg_control() --------------------------------------------------------

test_that("nbreg_control() returns validated settings", {
  ctrl <- nbreg_control()
  expect_s3_class(ctrl, "nbreg_control")
  expect_identical(ctrl$maxit, 100L)
  expect_identical(nbreg_control(maxit = 20)$maxit, 20L)
  expect_null(ctrl$theta_start)
  expect_equal(nbreg_control(theta_start = 2)$theta_start, 2)
})

test_that("nbreg_control() rejects invalid tolerances", {
  for (arg in c("tol", "irls_tol", "theta_tol")) {
    expect_input_error(do.call(nbreg_control, setNames(list(0), arg)),
                       paste0("`", arg, "` must be a number strictly between 0 and 1"))
    expect_input_error(do.call(nbreg_control, setNames(list(1), arg)),
                       "strictly between 0 and 1")
    expect_input_error(do.call(nbreg_control, setNames(list(-1e-8), arg)),
                       "strictly between 0 and 1")
    expect_input_error(do.call(nbreg_control, setNames(list("a"), arg)),
                       paste0("`", arg, "` must be a single finite number"))
    expect_input_error(do.call(nbreg_control, setNames(list(c(1e-8, 1e-9)), arg)),
                       "single finite number")
    expect_input_error(do.call(nbreg_control, setNames(list(NA_real_), arg)),
                       "single finite number")
    expect_input_error(do.call(nbreg_control, setNames(list(Inf), arg)),
                       "single finite number")
  }
})

test_that("nbreg_control() rejects invalid iteration limits", {
  for (arg in c("maxit", "irls_maxit", "theta_maxit", "max_halving")) {
    expect_input_error(do.call(nbreg_control, setNames(list(0), arg)),
                       paste0("`", arg, "` must be a whole number of at least 1"))
    expect_input_error(do.call(nbreg_control, setNames(list(2.5), arg)),
                       "whole number")
    expect_input_error(do.call(nbreg_control, setNames(list(-3), arg)),
                       "whole number")
    expect_input_error(do.call(nbreg_control, setNames(list(NULL), arg)),
                       "single finite number")
  }
})

test_that("nbreg_control() rejects iteration limits above the integer maximum", {
  big <- .Machine$integer.max
  for (arg in c("maxit", "irls_maxit", "theta_maxit", "max_halving")) {
    expect_input_error(do.call(nbreg_control, setNames(list(big + 1), arg)),
                       paste0("`", arg, "` must be at most 2147483647"))
    expect_input_error(do.call(nbreg_control, setNames(list(3e9), arg)),
                       "must be at most 2147483647")
    ctrl <- expect_silent(do.call(nbreg_control, setNames(list(big), arg)))
    expect_identical(ctrl[[arg]], big)
  }
})

test_that("nbreg_control() rejects an invalid trace flag", {
  expect_input_error(nbreg_control(trace = NA), "`trace` must be TRUE or FALSE")
  expect_input_error(nbreg_control(trace = "yes"), "`trace` must be TRUE or FALSE")
  expect_input_error(nbreg_control(trace = c(TRUE, FALSE)), "TRUE or FALSE")
  expect_input_error(nbreg_control(trace = 1), "TRUE or FALSE")
})

test_that("nbreg_control() rejects invalid theta bounds and thresholds", {
  expect_input_error(nbreg_control(theta_min = 0), "`theta_min` must be a positive number")
  expect_input_error(nbreg_control(theta_max = -1), "`theta_max` must be a positive number")
  expect_input_error(nbreg_control(theta_min = 10, theta_max = 5),
                     "`theta_min` .* must be smaller than `theta_max`")
  expect_input_error(nbreg_control(theta_min = 5, theta_max = 5),
                     "must be smaller than")
  expect_input_error(nbreg_control(theta_warn = 1e9),
                     "`theta_warn` .* must lie above `theta_min`")
  expect_input_error(nbreg_control(theta_warn = 1e-9), "`theta_warn`")
  expect_input_error(nbreg_control(theta_warn = 0), "`theta_warn` must be a positive number")
  expect_input_error(nbreg_control(theta_start = 0), "`theta_start` must be a positive number")
  expect_input_error(nbreg_control(theta_start = 1e9),
                     "`theta_start` .* must lie between `theta_min`")
  expect_input_error(nbreg_control(theta_start = "1"), "single finite number")
})

# ---- zinbreg() --------------------------------------------------------------

test_that("zinbreg() rejects bad input with informative errors", {
  skip("planned: zinbreg() implementation phase")
  # Planned cases: as for nbreg(), plus two-part formula problems and data
  # without any zeros.
})

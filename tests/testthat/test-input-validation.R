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

# ---- zinbreg_control() ------------------------------------------------------

test_that("zinbreg_control() returns validated settings", {
  ctrl <- zinbreg_control()
  expect_s3_class(ctrl, "zinbreg_control")
  expect_identical(names(ctrl), names(formals(zinbreg_control)))
  expect_identical(ctrl$maxit, 5000L)
  expect_identical(ctrl$zero_maxit, 50L)
  expect_identical(ctrl$hessian, "analytic")
  expect_identical(zinbreg_control(hessian = "numeric")$hessian, "numeric")
  expect_null(ctrl$start)
  expect_s3_class(ctrl$nb_control, "nbreg_control")
  expect_identical(ctrl$nb_control$tol, 1e-10)
  # a plain list for nb_control is checked and converted
  ctrl2 <- zinbreg_control(nb_control = list(irls_maxit = 7))
  expect_s3_class(ctrl2$nb_control, "nbreg_control")
  expect_identical(ctrl2$nb_control$irls_maxit, 7L)
  st <- zinbreg_control(start = list(count = c(a = 1L, b = 2L), zero = -1,
                                     theta = 3))$start
  expect_identical(st$count, c(a = 1, b = 2))
  expect_identical(st$theta, 3)
})

test_that("zinbreg_control() rejects invalid tolerances and limits", {
  for (arg in c("tol", "monotone_tol", "grad_tol", "zero_tol")) {
    expect_input_error(do.call(zinbreg_control, setNames(list(0), arg)),
                       paste0("`", arg, "` must be a number strictly between 0 and 1"))
    expect_input_error(do.call(zinbreg_control, setNames(list(1), arg)),
                       "strictly between 0 and 1")
    expect_input_error(do.call(zinbreg_control, setNames(list(NA_real_), arg)),
                       "single finite number")
    expect_input_error(do.call(zinbreg_control, setNames(list("a"), arg)),
                       "single finite number")
  }
  for (arg in c("maxit", "zero_maxit")) {
    expect_input_error(do.call(zinbreg_control, setNames(list(0), arg)),
                       paste0("`", arg, "` must be a whole number of at least 1"))
    expect_input_error(do.call(zinbreg_control, setNames(list(1.5), arg)),
                       "whole number")
    expect_input_error(do.call(zinbreg_control,
                               setNames(list(.Machine$integer.max + 1), arg)),
                       "must be at most 2147483647")
  }
})

test_that("zinbreg_control() rejects invalid flags and options", {
  expect_input_error(zinbreg_control(trace = NA), "`trace` must be TRUE or FALSE")
  expect_input_error(zinbreg_control(hessian_check = "yes"),
                     "`hessian_check` must be TRUE or FALSE")
  expect_input_error(zinbreg_control(hessian = "numerical"),
                     "`hessian` must be one of \"analytic\" or \"numeric\", not \"numerical\"")
  expect_input_error(zinbreg_control(hessian = "num"), "`hessian` must be one of")
  expect_input_error(zinbreg_control(hessian = NA_character_), "`hessian` must be one of")
  expect_input_error(zinbreg_control(hessian = 1), "`hessian` must be one of")
})

test_that("zinbreg_control() checks nb_control", {
  expect_input_error(zinbreg_control(nb_control = "fast"),
                     "`nb_control` must be created by `nbreg_control\\(\\)`")
  expect_input_error(zinbreg_control(nb_control = list(tolerance = 1)),
                     "Unknown setting.*in `nb_control`.*`tolerance`")
  expect_input_error(zinbreg_control(nb_control = list(tol = 2)),
                     "`tol` must be a number strictly between 0 and 1")
  bad <- nbreg_control()
  bad$irls_maxit <- -1
  expect_input_error(zinbreg_control(nb_control = bad),
                     "`irls_maxit` must be a whole number")
  extra <- nbreg_control()
  extra$speed <- 1
  expect_input_error(zinbreg_control(nb_control = extra),
                     "`nb_control` is a modified nbreg_control object")
})

test_that("zinbreg_control() checks start", {
  expect_input_error(zinbreg_control(start = c(1, 2)),
                     "`start` must be NULL or a list")
  expect_input_error(zinbreg_control(start = data.frame(count = 1)),
                     "`start` must be NULL or a list")
  expect_input_error(zinbreg_control(start = list(1, 2, 3)),
                     "Every element of `start` must be named")
  expect_input_error(zinbreg_control(start = list(count = 1, zero = 1)),
                     "`start` is missing `theta`")
  expect_input_error(
    zinbreg_control(start = list(count = 1, zero = 1, theta = 1, extra = 2)),
    "`start` has unknown element.*`extra`")
  expect_input_error(
    zinbreg_control(start = list(count = 1, count = 2, zero = 1, theta = 1)),
    "`start` gives `count` more than once")
  expect_input_error(zinbreg_control(start = list(count = "a", zero = 1, theta = 1)),
                     "`start\\$count` must be a non-empty numeric vector")
  expect_input_error(zinbreg_control(start = list(count = 1, zero = c(1, NA),
                                                  theta = 1)),
                     "`start\\$zero` must be a non-empty numeric vector of finite")
  expect_input_error(zinbreg_control(start = list(count = numeric(), zero = 1,
                                                  theta = 1)),
                     "`start\\$count` must be a non-empty")
  expect_input_error(zinbreg_control(start = list(count = 1, zero = 1, theta = 0)),
                     "`start\\$theta` must be a positive number")
  expect_input_error(zinbreg_control(start = list(count = 1, zero = 1,
                                                  theta = c(1, 2))),
                     "`start\\$theta` must be a single finite number")
  expect_input_error(zinbreg_control(start = list(count = 1, zero = 1, theta = 1e9)),
                     "`start\\$theta` .* must lie between `theta_min`")
})

test_that("zinbreg_control objects are revalidated by resolve_control()", {
  resolve <- function(x) {
    nbzinbreg:::resolve_control(x, zinbreg_control, "zinbreg_control")
  }
  ctrl <- zinbreg_control()
  expect_identical(resolve(ctrl), ctrl)
  bad <- ctrl
  bad$grad_tol <- 5
  expect_input_error(resolve(bad), "`grad_tol` must be a number strictly between")
  bad_nb <- ctrl
  bad_nb$nb_control$tol <- -1
  expect_input_error(resolve(bad_nb), "`tol` must be a number strictly between")
  extra <- ctrl
  extra$speed <- 1
  expect_input_error(resolve(extra),
                     "`control` is a modified zinbreg_control object with unknown field")
  expect_input_error(resolve(list(tolerance = 1)),
                     "Unknown setting.*`tolerance`.*`\\?zinbreg_control`")
  expect_input_error(resolve(list(tol = 1e-9, tol = 1e-10)), "more than once")
  expect_input_error(resolve("x"),
                     "`control` must be created by `zinbreg_control\\(\\)`")
  expect_identical(resolve(list(maxit = 10))$maxit, 10L)
})

# ---- nbreg(): arguments -----------------------------------------------------

d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4), t = 1:8)

test_that("nbreg() checks the formula", {
  expect_input_error(nbreg(), "`formula` is missing")
  expect_input_error(nbreg("y ~ x", d), "not a character string.*as.formula")
  expect_input_error(nbreg(1, d), "`formula` must be a formula")
  expect_input_error(nbreg(no_such_formula, d),
                     "`formula` could not be evaluated")
  expect_input_error(nbreg(~ x, d), "`formula` has no response")
  expect_input_error(nbreg(y ~ x | x, d), "Two-part formulas.*zinbreg")
})

test_that("nbreg() checks the data argument", {
  expect_input_error(nbreg(y ~ x, as.matrix(d)), "not a matrix.*as.data.frame")
  expect_input_error(nbreg(y ~ x, 5), "`data` must be a data frame")
  expect_input_error(nbreg(y ~ x, no_such_data), "`data` could not be evaluated")
  expect_input_error(nbreg(y ~ x, d[0, ]), "`data` has no rows")
})

test_that("nbreg() rejects unused arguments with specific hints", {
  expect_input_error(nbreg(y ~ x, d, family = "nb"),
                     "Unused argument.*`family`.*no `family` or `link`")
  expect_input_error(nbreg(y ~ x, d, link = "log"), "no `family` or `link`")
  expect_input_error(nbreg(y ~ x, d, init.theta = 2), "theta_start")
  expect_input_error(nbreg(y ~ x, d, tol = 1e-9),
                     "`tol`.*control = nbreg_control")
  expect_input_error(
    nbreg(y ~ x, d, NULL, nbreg_control(), NULL, , , TRUE, "extra"),
    "unnamed extra argument")
})

test_that("nbreg() checks model and control", {
  expect_input_error(nbreg(y ~ x, d, model = "yes"), "`model` must be TRUE or FALSE")
  expect_input_error(nbreg(y ~ x, d, model = NA), "`model` must be TRUE or FALSE")
  expect_input_error(nbreg(y ~ x, d, control = "fast"),
                     "`control` must be created by `nbreg_control")
  expect_input_error(nbreg(y ~ x, d, control = list(tolerance = 1)),
                     "Unknown setting.*`tolerance`")
  expect_input_error(nbreg(y ~ x, d, control = list(1e-9)), "must be named")
  expect_input_error(nbreg(y ~ x, d, control = list(tol = 1e-9, tol = 1e-10)),
                     "gives setting\\(s\\) `tol` more than once")
  expect_input_error(
    nbreg(y ~ x, d, control = list(tol = 1e-9, maxit = 5, tol = 1e-10,
                                   maxit = 6)),
    "`tol`, `maxit` more than once")
  expect_input_error(nbreg(y ~ x, d, control = list(tol = 2)),
                     "`tol` must be a number strictly between 0 and 1")
})

# ---- nbreg(): building the model frame --------------------------------------

test_that("nbreg() explains model-frame problems", {
  expect_input_error(nbreg(y ~ zz, d), "Variable `zz` was not found")
  expect_input_error(nbreg(y ~ x, d, offset = 1:3),
                     "`offset` must have one value per row")
  expect_input_error(nbreg(y ~ x, d, weights = 1:3),
                     "`weights` must have one value per row")
  dn <- d
  dn$x[3] <- NA
  expect_input_error(nbreg(y ~ x, dn, na.action = stats::na.fail),
                     "na.action = na.fail")
  expect_input_error(nbreg(y ~ x, d, subset = x > 100),
                     "No observations are left")
})

# ---- nbreg(): response, weights, offset -------------------------------------

test_that("nbreg() checks the response", {
  dd <- d
  dd$f <- factor(dd$y)
  dd$ch <- as.character(dd$y)
  dd$lg <- dd$y > 2
  expect_input_error(nbreg(f ~ x, dd), "must be numeric counts, not a factor")
  expect_input_error(nbreg(ch ~ x, dd), "must be numeric counts, not a character")
  expect_input_error(nbreg(lg ~ x, dd), "must be numeric counts, not a logical")
  expect_input_error(nbreg(cbind(y, y) ~ x, d), "single column.*2 columns")
  expect_input_error(nbreg(I(y - 2) ~ x, d),
                     "non-negative counts. Found 3 negative value.*-1, -2 and -1")
  expect_input_error(nbreg(I(y + 0.25) ~ x, d),
                     "whole numbers .counts.. Found 8 non-integer value.*1.25")
  expect_input_error(nbreg(I(c(Inf, y[-1])) ~ x, d), "1 infinite value")
  expect_input_error(nbreg(I(c(NA, y[-1])) ~ x, d, na.action = stats::na.pass),
                     "1 missing value")
  expect_input_error(nbreg(y ~ x, data.frame(x = 1:8, y = 0)),
                     "Every response value .* is zero")
})

test_that("nbreg() checks the weights", {
  expect_input_error(nbreg(y ~ x, d, weights = rep("a", 8)),
                     "`weights` must be numeric")
  expect_input_error(nbreg(y ~ x, d, weights = c(-1, rep(1, 7))),
                     "`weights` must be non-negative. Found 1 negative")
  expect_input_error(nbreg(y ~ x, d, weights = c(Inf, rep(1, 7))),
                     "`weights` contains 1 infinite")
  expect_input_error(nbreg(y ~ x, d, weights = c(NA, rep(1, 7)),
                           na.action = stats::na.pass),
                     "`weights` contains 1 missing")
  expect_input_error(nbreg(y ~ x, d, weights = rep(0, 8)), "All `weights` are zero")
})

test_that("nbreg() checks the offset", {
  expect_input_error(nbreg(y ~ x, d, offset = rep("a", 8)),
                     "offset must be a numeric vector")
  expect_input_error(nbreg(y ~ x, d, offset = log(c(0, 1:7))),
                     "offset contains 1 missing or infinite")
  expect_input_error(nbreg(y ~ x + offset(log(t - 1)), d),
                     "offset contains 1 missing or infinite")
})

# ---- nbreg(): model matrix --------------------------------------------------

test_that("nbreg() checks the model matrix", {
  expect_input_error(nbreg(y ~ 0, d), "no coefficients")
  expect_input_error(nbreg(y ~ I(c(Inf, x[-1])), d),
                     "contain infinite values")
  expect_error(nbreg(y ~ x + I(2 * x), d), class = "nbzinbreg_error_rank",
               regexp = "`I\\(2 \\* x\\)` is a linear combination")
  expect_input_error(nbreg(y ~ x, d[1:2, ]), "Not enough data")
  expect_input_error(nbreg(y ~ x, d, weights = c(1, 1, rep(0, 6))),
                     "Not enough data: 2 observation")
})

test_that("nbreg() explains model-matrix failures and names the factor", {
  dg <- d
  dg$g <- factor(rep(c("a", "b"), 4))
  dg$h <- factor(rep(c("p", "q"), each = 4))
  dg$ch <- "k"
  # a factor reduced to one level by `subset`
  expect_input_error(nbreg(y ~ x + g, dg, subset = g == "a"),
                     "Factor\\(s\\) `g` has only one level.*subset")
  # a constant character covariate
  expect_input_error(nbreg(y ~ x + ch, dg), "`ch` has only one level")
  # two factors reduced to one level each
  expect_input_error(nbreg(y ~ g + h, dg, subset = g == "a" & h == "p"),
                     "`g`, `h` have only one level")
  # the original R message is kept for reference. R's own wording differs
  # between versions ("contrasts can be applied only to factors with 2 or
  # more levels", and in R-devel "cannot set contrasts for variable 'ch'"),
  # so only the fact that it is kept, and that it concerns contrasts, is checked
  expect_input_error(nbreg(y ~ x + ch, dg), "Original error: .*contrasts")
  expect_match("Original error: contrasts can be applied only to factors with 2 or more levels",
               "Original error: .*contrasts")
  expect_match("Original error: cannot set contrasts for variable 'ch':",
               "Original error: .*contrasts")
})

test_that("nbreg() revalidates nbreg_control objects", {
  ctrl <- nbreg_control()
  bad_value <- ctrl
  bad_value$tol <- -1
  expect_input_error(nbreg(y ~ x, d, control = bad_value),
                     "`tol` must be a number strictly between 0 and 1")
  bad_whole <- ctrl
  bad_whole$maxit <- 2.5
  expect_input_error(nbreg(y ~ x, d, control = bad_whole),
                     "`maxit` must be a whole number")
  bad_bounds <- ctrl
  bad_bounds$theta_min <- 1e9
  expect_input_error(nbreg(y ~ x, d, control = bad_bounds),
                     "`theta_min` .* must be smaller than `theta_max`")
  extra <- ctrl
  extra$speed <- "fast"
  expect_input_error(nbreg(y ~ x, d, control = extra),
                     "modified nbreg_control object with unknown field.*`speed`")
  missing_field <- ctrl
  missing_field[["irls_tol"]] <- NULL
  expect_input_error(nbreg(y ~ x, d, control = missing_field),
                     "missing field.*`irls_tol`")
  both <- missing_field
  both$speed <- 1
  expect_input_error(nbreg(y ~ x, d, control = both),
                     "missing field.*`irls_tol`.* and unknown field.*`speed`")
  expect_input_error(nbreg(y ~ x, d,
                           control = structure("x", class = "nbreg_control")),
                     "not a list")
  expect_input_error(nbreg(y ~ x, d,
                           control = structure(list(1, 2),
                                               class = "nbreg_control")),
                     "malformed nbreg_control")
  # a valid edit is accepted, checked and used
  good <- ctrl
  good$tol <- 1e-10
  expect_identical(nbreg(y ~ x, d, control = good)$control$tol, 1e-10)
})

test_that("every nbreg() error carries the parent class", {
  expect_error(nbreg(), class = "nbzinbreg_error")
  expect_error(nbreg(y ~ x + I(2 * x), d), class = "nbzinbreg_error")
})

# ---- zinbreg() --------------------------------------------------------------

dz <- data.frame(
  x = rep(1:10, 3) / 10,
  z = rep(c(0.3, -0.1, 0.8, 0.2, -0.5), 6),
  y = c(0, 1, 0, 3, 0, 2, 0, 0, 5, 1, 0, 4, 0, 1, 0, 0, 2, 7, 0, 1,
        0, 3, 0, 0, 1, 6, 0, 2, 0, 1),
  t = rep(1:3, 10)
)

test_that("zinbreg() checks the formula", {
  expect_input_error(zinbreg(), "`formula` is missing")
  expect_input_error(zinbreg("y ~ x", dz), "not a character string")
  expect_input_error(zinbreg(1, dz), "`formula` must be a formula")
  expect_input_error(zinbreg(no_such_formula, dz),
                     "`formula` could not be evaluated")
  expect_input_error(zinbreg(~ x | z, dz), "`formula` has no response")
  expect_input_error(zinbreg(y ~ x | z | t, dz),
                     "has 3 parts on the right of `~`.*at most\\s+two")
  expect_input_error(zinbreg(y | t ~ x | z, dz),
                     "exactly one response on the left of `~`, not 2 parts")
  expect_input_error(zinbreg(y ~ x | 0, dz),
                     "zero part has no coefficients.*`\\| 1`.*`nbreg\\(\\)`")
  expect_input_error(zinbreg(y ~ x | -1, dz), "zero part has no coefficients")
  expect_input_error(zinbreg(y ~ 0 | z, dz),
                     "count part has no coefficients")
})

test_that("zinbreg() checks data, dots, model and control", {
  expect_input_error(zinbreg(y ~ x | z, as.matrix(dz)), "not a matrix")
  expect_input_error(zinbreg(y ~ x | z, 5), "`data` must be a data frame")
  expect_input_error(zinbreg(y ~ x | z, dz[0, ]), "`data` has no rows")
  expect_input_error(zinbreg(y ~ x | z, dz, dist = "negbin"),
                     "Unused argument.*`dist`.*no `dist`, `link` or `family`")
  expect_input_error(zinbreg(y ~ x | z, dz, link = "logit"),
                     "zero-inflated negative binomial")
  expect_input_error(zinbreg(y ~ x | z, dz, init.theta = 1),
                     "zinbreg_control\\(start =")
  expect_input_error(zinbreg(y ~ x | z, dz, EM = TRUE),
                     "`EM`.*control = zinbreg_control")
  expect_input_error(zinbreg(y ~ x | z, dz, model = NA),
                     "`model` must be TRUE or FALSE")
  expect_input_error(zinbreg(y ~ x | z, dz, control = "fast"),
                     "`control` must be created by `zinbreg_control\\(\\)`")
  # the other model's control object
  expect_input_error(zinbreg(y ~ x | z, dz, control = nbreg_control()),
                     paste0("created by `nbreg_control\\(\\)`, but ",
                            "`zinbreg_control\\(\\)` is needed.*",
                            "nb_control = nbreg_control"))
  expect_input_error(nbreg(y ~ x, d, control = zinbreg_control()),
                     paste0("created by `zinbreg_control\\(\\)`, but ",
                            "`nbreg_control\\(\\)` is needed here.$"))
  expect_input_error(zinbreg(y ~ x | z, dz, control = list(tol = -1)),
                     "`tol` must be a number strictly between 0 and 1")
  expect_input_error(zinbreg(y ~ x | z, dz, control = list(speed = 1)),
                     "Unknown setting.*`speed`")
  bad <- zinbreg_control()
  bad$nb_control$theta_max <- -1
  expect_input_error(zinbreg(y ~ x | z, dz, control = bad),
                     "`theta_max` must be a positive number")
})

test_that("zinbreg() explains model-frame problems in either part", {
  expect_input_error(zinbreg(y ~ x | nope, dz), "Variable `nope` was not found")
  expect_input_error(zinbreg(y ~ nope | z, dz), "Variable `nope` was not found")
  expect_input_error(zinbreg(y ~ x | z, dz, offset = 1:3),
                     "`offset` must have one value per row")
  expect_input_error(zinbreg(y ~ x | z, dz, weights = 1:3),
                     "`weights` must have one value per row")
  dn <- dz
  dn$z[4] <- NA
  expect_input_error(zinbreg(y ~ x | z, dn, na.action = stats::na.fail),
                     "na.action = na.fail")
  expect_input_error(zinbreg(y ~ x | z, dz, subset = x > 100),
                     "No observations are left")
})

test_that("zinbreg() checks the response, weights and offsets", {
  expect_input_error(zinbreg(I(y - 1) ~ x | z, dz), "non-negative counts")
  expect_input_error(zinbreg(I(y + 0.5) ~ x | z, dz), "whole numbers")
  dd <- dz
  dd$f <- factor(dd$y)
  expect_input_error(zinbreg(f ~ x | z, dd), "not a factor")
  expect_input_error(zinbreg(y ~ x | z, transform(dz, y = 0)),
                     "Every response value .* is zero")
  expect_input_error(zinbreg(y ~ x | z, dz, weights = c(-1, rep(1, 29))),
                     "`weights` must be non-negative")
  expect_input_error(zinbreg(y ~ x | z, dz, weights = rep(0, 30)),
                     "All `weights` are zero")
  expect_input_error(zinbreg(y ~ x | z, dz, offset = rep("a", 30)),
                     "offset.*numeric")
  expect_input_error(zinbreg(y ~ x | z, dz, offset = log(t - 1)),
                     "offset contains 10 missing or infinite")
  expect_input_error(zinbreg(y ~ x + offset(log(t - 1)) | z, dz),
                     "^The offset contains 10 missing or infinite")
  expect_input_error(zinbreg(y ~ x | z + offset(log(t - 1)), dz),
                     "^The zero-part offset contains 10 missing or infinite")
})

test_that("zinbreg() needs zeros among the observations with positive weight", {
  no_zeros <- transform(dz, y = y + 1)
  expect_input_error(zinbreg(y ~ x | z, no_zeros),
                     "no zeros among the observations with positive weight.*`nbreg\\(\\)`")
  # zeros that carry weight 0 do not count
  w <- ifelse(dz$y == 0, 0, 1)
  expect_input_error(zinbreg(y ~ x | z, dz, weights = w),
                     sprintf(paste0("positive\\s+weight \\(%d zero\\(s\\) have ",
                                    "weight 0 and do not count\\)"),
                             sum(dz$y == 0)))
  # one zero with positive weight is enough to try
  w[which(dz$y == 0)[1]] <- 1
  expect_error(suppressWarnings(zinbreg(y ~ x | z, dz, weights = w)), NA)
})

test_that("zinbreg() checks both model matrices", {
  dd <- dz
  dd$z2 <- 2 * dd$z
  expect_error(zinbreg(y ~ x | z + z2, dd), class = "nbzinbreg_error_rank",
               regexp = "zero-part model matrix is rank-deficient: `z2`")
  dd$x2 <- dd$x * 3
  expect_error(zinbreg(y ~ x + x2 | z, dd), class = "nbzinbreg_error_rank",
               regexp = "count-part model matrix is rank-deficient: `x2`")
  expect_input_error(zinbreg(y ~ x | I(c(Inf, z[-1])), dz),
                     "of the zero part contain infinite values")
  dd$g <- factor(rep("a", 30))
  expect_input_error(zinbreg(y ~ x | z + g, dd), "`g` has only one level")
  expect_input_error(zinbreg(y ~ x | z, dz[1:4, ]), "Not enough data")
  expect_input_error(
    zinbreg(y ~ x | z, dz, weights = c(rep(1, 5), rep(0, 25))),
    "Not enough data: 5 observation\\(s\\) with positive weight for 5")
})

test_that("zinbreg() checks start against the model", {
  st <- function(count, zero) {
    zinbreg_control(start = list(count = count, zero = zero, theta = 1))
  }
  expect_input_error(zinbreg(y ~ x | z, dz, control = st(1, c(0, 0))),
                     "`start\\$count` has 1 value\\(s\\), but the count part has 2")
  expect_input_error(zinbreg(y ~ x | z, dz, control = st(c(0, 0), 1:3)),
                     "`start\\$zero` has 3 value\\(s\\), but the zero part has 2")
  expect_input_error(zinbreg(y ~ x | z, dz,
                             control = st(c(a = 0, b = 0), c(0, 0))),
                     "names of `start\\$count` .* do not match")
  expect_input_error(zinbreg(y ~ x | z, dz,
                             control = st(c(x = 0, x = 0), c(0, 0))),
                     "do not match")
  # an invalid starting point for the fit
  expect_input_error(zinbreg(y ~ x | z, dz, control = st(c(800, 0), c(0, 0))),
                     "starting values give an invalid log-likelihood")
})

test_that("every zinbreg() error carries the parent class", {
  expect_error(zinbreg(), class = "nbzinbreg_error")
  expect_error(zinbreg(y ~ x | z + I(2 * z), dz), class = "nbzinbreg_error")
})

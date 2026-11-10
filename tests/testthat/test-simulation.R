# The simulation study's functions (inst/benchmark/simulation_functions.R):
# the data generator, the per-replication records and classes, the profile
# classes, the summaries with Monte Carlo SEs, and saving.

sim_env <- function() {
  path <- system.file("benchmark", "simulation_functions.R",
                      package = "nbzinbreg")
  if (!nzchar(path)) testthat::skip("simulation functions not available")
  env <- new.env()
  sys.source(path, envir = env, keep.source = FALSE)
  env
}

# A fake zinbreg fit with the fields the record uses, from a real fit.
real_fit <- function(s, seed = 1) {
  des <- s$sim_design(500, 651)
  y <- s$sim_response(des, s$sim_truth(), seed)$y
  zinbreg(y ~ x1 + x2 | x1 + z, data = data.frame(y = y, des))
}

# ---- the data generator --------------------------------------------------------

test_that("the design and responses are reproducible and leave the RNG alone", {
  s <- sim_env()
  set.seed(99)
  state <- .Random.seed
  d1 <- s$sim_design(500, 651)
  expect_identical(.Random.seed, state)
  expect_identical(d1, s$sim_design(500, 651))
  expect_false(identical(d1, s$sim_design(500, 652)))
  y1 <- s$sim_response(d1, s$sim_truth(), 652)
  expect_identical(y1, s$sim_response(d1, s$sim_truth(), 652))
  expect_identical(.Random.seed, state)
  expect_identical(names(d1), c("x1", "x2", "z"))
  expect_true(all(d1$x2 %in% 0:1))
  # every structural zero is a zero
  expect_true(all(y1$y[y1$structural == 1] == 0))
})

test_that("large samples match the design's exact moments", {
  s <- sim_env()
  tr <- s$sim_truth()
  des <- s$sim_design(200000, 7)
  # covariate distributions
  expect_lt(abs(mean(des$x1)), 4 / sqrt(2e5))
  expect_lt(abs(stats::sd(des$x1) - 1), 0.01)
  expect_lt(abs(mean(des$x2) - 0.5), 4 * 0.5 / sqrt(2e5))
  m <- s$design_moments(des, tr)
  sim <- s$sim_response(des, tr, 8)
  n <- nrow(des)
  p0 <- m[["zero_share"]]
  expect_lt(abs(mean(sim$y == 0) - p0), 4 * sqrt(p0 * (1 - p0) / n))
  ps <- m[["structural_share"]]
  expect_lt(abs(mean(sim$structural) - ps), 4 * sqrt(ps * (1 - ps) / n))
  expect_lt(abs(mean(sim$y) - m[["mean_y"]]),
            4 * stats::sd(sim$y) / sqrt(n))
  # the formulas themselves, by hand at one row
  lin <- s$sim_linear(des[1, ], tr)
  mu <- exp(1 + 0.5 * des$x1[1] - 0.5 * des$x2[1])
  pi <- stats::plogis(-0.8 + 0.6 * des$x1[1] - 0.5 * des$z[1])
  expect_equal(c(lin$mu, lin$pi, lin$f0),
               c(mu, pi, (1.5 / (1.5 + mu))^1.5))
})

test_that("the true values are named as coef() of the fit names them", {
  s <- sim_env()
  tr <- s$sim_truth()
  expect_identical(names(tr$params),
                   c("count_(Intercept)", "count_x1", "count_x2",
                     "zero_(Intercept)", "zero_x1", "zero_z", "log(theta)",
                     "theta"))
  expect_equal(tr$params[["log(theta)"]], log(1.5))
  fit <- real_fit(s)
  expect_identical(names(coef(fit)), names(tr$params)[1:6])
})

# ---- records and classes ---------------------------------------------------------

test_that("a record holds the fit's estimates, SEs and Wald intervals", {
  s <- sim_env()
  fit <- real_fit(s)
  rec <- s$replication_record(3, fit, truth = s$sim_truth())
  expect_identical(rec$class, "usable")
  expect_equal(rec$estimate[1:6], coef(fit))
  expect_equal(rec$se[1:6], sqrt(diag(fit$vcov)), ignore_attr = TRUE)
  expect_equal(rec$estimate[["log(theta)"]], log(fit$theta))
  expect_equal(rec$se[["log(theta)"]], fit$SE.logtheta)
  expect_equal(rec$se[["theta"]], fit$SE.theta)
  z <- stats::qnorm(0.975)
  expect_equal(rec$lower, rec$estimate - z * rec$se)
  expect_equal(rec$upper, rec$estimate + z * rec$se)
  tv <- s$sim_truth()$params
  expect_identical(rec$covers, rec$lower <= tv & tv <= rec$upper)
  expect_true(all(rec$valid))
})

test_that("each replication class is detected with its reason", {
  s <- sim_env()
  fit <- real_fit(s)
  cls <- function(f, ...) {
    r <- s$replication_record(1, f, ...)
    c(r$class, r$reason)
  }
  expect_identical(cls(NULL, error = "boom"), c("error", "boom"))
  f <- fit
  f$converged <- FALSE
  expect_identical(cls(f)[1], "not converged")
  f <- fit
  f$gradient_check$ok <- FALSE
  expect_match(cls(f)[2], "gradient check passed = FALSE")
  f <- fit
  f$diverging <- "zero_z"
  expect_identical(cls(f), c("boundary", "diverging: zero_z"))
  f <- fit
  f$theta_bound <- "upper"
  expect_identical(cls(f), c("boundary", "theta on its upper bound"))
  f <- fit
  f$pi[] <- 1e-10
  expect_match(cls(f)[2], "zero part collapsed")
  f <- fit
  f$vcov_full[1, 2] <- f$vcov_full[1, 2] + 1
  expect_identical(cls(f)[1], "invalid covariance")
  expect_match(cls(f)[2], "not symmetric")
  f <- fit
  f$SE.logtheta <- NaN
  expect_identical(cls(f), c("invalid values", "not valid: log(theta)"))
  f <- fit
  f$SE.theta <- -1
  expect_identical(cls(f), c("invalid values", "not valid: theta"))
  # a boundary replication keeps its valid parameters flagged valid
  f <- fit
  f$diverging <- "zero_z"
  f$vcov["zero_z", "zero_z"] <- NA
  rec <- s$replication_record(1, f)
  expect_false(rec$valid[["zero_z"]])
  expect_true(all(rec$valid[c("count_(Intercept)", "count_x1", "count_x2")]))
})

test_that("profile intervals are classified, and unbounded ones are evaluated", {
  s <- sim_env()
  tr <- s$sim_truth()
  ci <- function(lo, hi, st) {
    m <- matrix(c(lo, hi), 1, dimnames = list("theta", c("2.5 %", "97.5 %")))
    structure(m, status = matrix(st, 1, dimnames = list("theta", NULL)),
              class = "nbzinbreg_confint")
  }
  p <- s$profile_record(ci(1, 2, c("finite", "finite")), truth = tr)
  expect_identical(c(p$class, p$covers), c("finite", "TRUE"))
  p <- s$profile_record(ci(1.6, 2, c("finite", "finite")), truth = tr)
  expect_false(p$covers)
  p <- s$profile_record(ci(1.2, Inf, c("finite", "unbounded")), truth = tr)
  expect_identical(p$class, "unbounded")
  expect_true(p$covers)
  p <- s$profile_record(ci(1.7, Inf, c("finite", "unbounded")), truth = tr)
  expect_false(p$covers)
  p <- s$profile_record(ci(1.2, NA, c("finite", "beyond theta_max")), truth = tr)
  expect_identical(p$class, "unresolved")
  expect_true(is.na(p$covers))
  expect_match(p$detail, "beyond theta_max")
  p <- s$profile_record(NULL, error = "no", truth = tr)
  expect_identical(c(p$class, p$detail), c("error", "no"))
  p <- s$profile_record(NULL, fit_usable = FALSE, truth = tr)
  expect_identical(p$class, "fit not usable")
})

# ---- summaries -------------------------------------------------------------------

# A synthetic records table with known values: estimates truth + noise,
# SEs, intervals and coverage flags as a real record would have them.
synthetic_table <- function(s, R = 400, seed = 3, extra_classes = character()) {
  tr <- s$sim_truth()
  pn <- names(tr$params)
  set.seed(seed)
  tab <- data.frame(replication = seq_len(R), class = "usable",
                    reason = NA_character_, stringsAsFactors = FALSE)
  for (p in pn) {
    e <- tr$params[[p]] + stats::rnorm(R, 0.01, 0.1)
    se <- stats::runif(R, 0.08, 0.12)
    lo <- e - 1.96 * se
    hi <- e + 1.96 * se
    tab[[paste0("estimate_", p)]] <- e
    tab[[paste0("se_", p)]] <- se
    tab[[paste0("lower_", p)]] <- lo
    tab[[paste0("upper_", p)]] <- hi
    tab[[paste0("covers_", p)]] <- lo <= tr$params[[p]] & tr$params[[p]] <= hi
    tab[[paste0("valid_", p)]] <- TRUE
  }
  if (length(extra_classes)) {
    tab$class[seq_along(extra_classes)] <- extra_classes
  }
  tab
}

test_that("the summaries equal hand computations", {
  s <- sim_env()
  tab <- synthetic_table(s)
  tr <- s$sim_truth()
  sm <- s$summarise_simulation(tab, tr, B = 300, seed = 5)
  p <- "zero_x1"
  e <- tab[[paste0("estimate_", p)]]
  se <- tab[[paste0("se_", p)]]
  cv <- tab[[paste0("covers_", p)]]
  t <- tr$params[[p]]
  row <- sm[sm$parameter == p, ]
  R <- nrow(tab)
  expect_equal(row$mean, mean(e))
  expect_equal(row$bias, mean(e) - t)
  expect_equal(row$bias_mcse, stats::sd(e) / sqrt(R))
  expect_equal(row$rel_bias, (mean(e) - t) / t)
  expect_equal(row$mc_sd, stats::sd(e))
  expect_equal(row$mc_sd_mcse_normal, stats::sd(e) / sqrt(2 * (R - 1)))
  expect_equal(row$mean_se, mean(se))
  expect_equal(row$mean_se_mcse, stats::sd(se) / sqrt(R))
  expect_equal(row$se_ratio, mean(se) / stats::sd(e))
  expect_equal(row$coverage, mean(cv))
  expect_equal(row$coverage_mcse, sqrt(mean(cv) * (1 - mean(cv)) / R))
  # RMSE is computed directly, not as sqrt(bias^2 + sd^2)
  expect_equal(row$rmse, sqrt(mean((e - t)^2)))
  expect_equal(row$rmse^2,
               (mean(e) - t)^2 + stats::var(e) * (R - 1) / R)
  expect_false(isTRUE(all.equal(row$rmse, sqrt((mean(e) - t)^2 + stats::var(e)))))
  # no non-usable replications: the bounds equal the coverage
  expect_equal(row$coverage_lower_bound, row$coverage)
  expect_equal(row$coverage_upper_bound, row$coverage)
  expect_identical(row$n_valid_converged, as.numeric(R))
})

test_that("bootstrap MC SEs are reproducible and agree with analytic ones", {
  s <- sim_env()
  tab <- synthetic_table(s, R = 1000, seed = 11)
  a <- s$summarise_simulation(tab, B = 1000, seed = 5)
  b <- s$summarise_simulation(tab, B = 1000, seed = 5)
  expect_identical(a, b)
  c2 <- s$summarise_simulation(tab, B = 1000, seed = 6)
  expect_false(identical(a$rmse_mcse, c2$rmse_mcse))
  # normal synthetic estimates: the bootstrap SE of the SD is close to the
  # normal-theory value, and every bootstrap SE is positive and finite
  expect_true(all(abs(a$mc_sd_mcse / a$mc_sd_mcse_normal - 1) < 0.25))
  for (col in c("rel_bias_mcse", "mc_sd_mcse", "se_ratio_mcse", "rmse_mcse")) {
    expect_true(all(is.finite(a[[col]]) & a[[col]] > 0), info = col)
  }
  # the relative bias's bootstrap SE is about bias_mcse / |truth|
  i <- which(a$parameter == "count_x1")
  expect_lt(abs(a$rel_bias_mcse[i] / (a$bias_mcse[i] / 0.5) - 1), 0.15)
})

test_that("non-usable replications give coverage bounds and are counted", {
  s <- sim_env()
  tab <- synthetic_table(s, R = 100,
                         extra_classes = c("boundary", "error", "not converged"))
  sm <- s$summarise_simulation(tab, B = 50, seed = 1)
  p <- "count_x1"
  cv <- tab[[paste0("covers_", p)]][tab$class == "usable"]
  row <- sm[sm$parameter == p, ]
  expect_identical(row$n_usable, 97L)
  expect_equal(row$coverage, mean(cv))
  expect_equal(row$coverage_lower_bound, sum(cv) / 100)
  expect_equal(row$coverage_upper_bound, (sum(cv) + 3) / 100)
  # valid counts are over the converged replications (boundary included)
  expect_identical(row$n_valid_converged, 98)
  cl <- s$class_summary(tab)
  expect_identical(as.integer(cl$counts[c("usable", "boundary", "error",
                                          "not converged")]),
                   c(97L, 1L, 1L, 1L))
  expect_identical(cl$non_usable$replication, 1:3)
})

test_that("profile coverage and its bounds equal hand computations", {
  s <- sim_env()
  tab <- data.frame(
    profile_class = c("finite", "finite", "unbounded", "unresolved", "error",
                      "fit not usable", NA),
    profile_covers = c(TRUE, FALSE, TRUE, NA, NA, NA, NA),
    profile_seconds = c(1, 2, 3, 4, 5, NA, NA))
  ps <- s$profile_summary(tab)
  expect_identical(ps$n, 6L)
  expect_identical(ps$evaluable, 3L)
  expect_equal(ps$coverage, 2 / 3)
  expect_equal(ps$coverage_mcse, sqrt((2 / 3) * (1 / 3) / 3))
  expect_equal(ps$lower_bound, 2 / 6)
  expect_equal(ps$upper_bound, 5 / 6)
  expect_identical(as.integer(ps$counts), c(2L, 1L, 1L, 1L, 1L))
})

# ---- running and saving ------------------------------------------------------------

test_that("a small real run is complete in structure and reproducible", {
  skip_on_cran()
  s <- sim_env()
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  if (!dir.exists(file.path(root, "R"))) skip("package source tree not available")
  r1 <- s$run_simulation(R = 4L, n = 300L, seed = 21L, profile_R = 2L, B = 50L,
                         source_path = root, progress = FALSE)
  r2 <- s$run_simulation(R = 4L, n = 300L, seed = 21L, profile_R = 2L, B = 50L,
                         source_path = root, progress = FALSE)
  cols <- grep("^(estimate|se|covers)_", names(r1$replicates), value = TRUE)
  expect_identical(r1$replicates[cols], r2$replicates[cols])
  expect_identical(nrow(r1$replicates), 4L)
  expect_identical(sum(!is.na(r1$replicates$profile_class)), 2L)
  expect_identical(r1$summary$parameter, names(s$sim_truth()$params))
  expect_identical(r1$settings$rng_kind, RNGkind())
  expect_named(r1$design_moments, c("zero_share", "structural_share_of_zeros",
                                    "structural_share", "mean_y"))
  # saving: every file, with the status line first
  out <- tempfile("sim_out")
  saved <- s$save_simulation(r1, out)
  expect_setequal(list.files(out), c("simulation_results.rds",
                                     "simulation_summary.csv",
                                     "simulation_replicates.csv",
                                     "simulation_status.txt"))
  first <- readLines(file.path(out, "simulation_status.txt"))[1]
  expect_true(first %in% c("COMPLETE", "FAILED"))
  expect_identical(first == "COMPLETE", saved$status == "complete")
  back <- readRDS(file.path(out, "simulation_results.rds"))
  expect_identical(back$replicates, r1$replicates)
})

test_that("a run with unverified provenance or missing replications fails", {
  s <- sim_env()
  ok <- list(provenance = list(status = "verified", reasons = character()),
             replicates = data.frame(x = 1:3), settings = list(R = 3L))
  expect_identical(s$simulation_status(ok)$status, "complete")
  bad <- ok
  bad$provenance <- list(status = "unverified", reasons = "uncommitted changes")
  st <- s$simulation_status(bad)
  expect_identical(st$status, "failed")
  expect_match(st$reasons, "uncommitted changes")
  short <- ok
  short$settings$R <- 5L
  expect_match(s$simulation_status(short)$reasons, "only 3 of 5")
})

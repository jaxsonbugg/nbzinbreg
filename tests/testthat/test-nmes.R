# The NMES1988 analysis functions (inst/benchmark/nmes_functions.R): the
# descriptive summary, the model comparison, the effect tables, and the
# checkpointed bootstrap LR test with its exact resumption and its two
# separate verdicts (run status and inference).

nmes_env <- function() {
  path <- system.file("benchmark", "nmes_functions.R", package = "nbzinbreg")
  if (!nzchar(path)) testthat::skip("NMES functions not available")
  env <- new.env()
  sys.source(path, envir = env, keep.source = FALSE)
  env
}

# The NMES1988 fits take a few seconds; fit them once for the whole file.
nmes_cache <- new.env()
nmes_cached_fits <- function(s) {
  testthat::skip_if_not_installed("AER")
  if (is.null(nmes_cache$fits)) nmes_cache$fits <- s$nmes_fits()
  nmes_cache$fits
}

# A small zero-inflated data set (as in test-vuong.R) on which a bootstrap
# replicate takes a fraction of a second.
boot_fixture <- function() {
  set.seed(12)
  n <- 150
  x <- stats::rnorm(n)
  y <- stats::rnbinom(n, size = 2, mu = exp(1 + 0.5 * x))
  y[stats::runif(n) < stats::plogis(-1 + 0.5 * x)] <- 0
  d <- data.frame(y = y, x = x)
  list(nb = nbreg(y ~ x, d), zinb = zinbreg(y ~ x | x, d))
}

fake_prov <- function(status = "verified", fingerprint = "fp") {
  list(status = status,
       reasons = if (status != "verified") "test: unverified" else character(),
       fingerprint_installed = fingerprint, fingerprint_source = fingerprint,
       helper_md5 = c(benchmark_functions.R = "a", nmes_functions.R = "b"),
       git_commit = "test", git_clean = TRUE,
       r_version = R.version.string, platform = "test")
}

# ---- descriptives and the two Poisson zero counts -----------------------------

test_that("the unconditional Poisson zeros are n exp(-mean) and labelled as such", {
  s <- nmes_env()
  testthat::skip_if_not_installed("AER")
  d <- s$nmes_descriptives()
  y <- s$nmes_data()$visits
  expect_identical(d$n, length(y))
  expect_identical(d$observed_zeros, sum(y == 0))
  expect_equal(d$poisson_unconditional_zeros, length(y) * exp(-mean(y)))
  expect_equal(d$variance_to_mean, stats::var(y) / mean(y))
  expect_match(d$poisson_unconditional_label, "ignoring covariates")
  # the frequency table covers every observation, and the Poisson column
  # every probability
  expect_identical(as.integer(sum(d$freq$observed)), length(y))
  expect_equal(sum(d$freq$poisson_unconditional), length(y))
  expect_identical(d$freq$visits[nrow(d$freq)], "30+")
})

test_that("the comparison uses the fitted models' expected zeros", {
  s <- nmes_env()
  f <- nmes_cached_fits(s)
  cmp <- s$nmes_comparison(f)
  expect_identical(cmp$model, c("poisson", "nb", "zinb", "zinb_health"))
  expect_identical(anyDuplicated(cmp$label), 0L)
  # the Poisson regression: the sum of each person's exp(-mu_i)
  expect_equal(cmp$expected_zeros[1], sum(exp(-stats::fitted(f$poisson))))
  # ... which is not the unconditional reference
  d <- s$nmes_descriptives()
  expect_gt(abs(cmp$expected_zeros[1] - d$poisson_unconditional_zeros), 1)
  # NB and ZINB agree with the package's expected_counts()
  expect_equal(cmp$expected_zeros[2], expected_counts(f$nb)$expected[1])
  expect_equal(cmp$expected_zeros[3], expected_counts(f$zinb)$expected[1])
  # and with the ZINB zero probability written out by hand
  p0 <- f$zinb$pi + (1 - f$zinb$pi) *
    stats::dnbinom(0, size = f$zinb$theta, mu = f$zinb$mu)
  expect_equal(cmp$expected_zeros[3], sum(p0))
  expect_equal(cmp$AIC, unname(vapply(f, stats::AIC, 1)))
  expect_equal(cmp$BIC, unname(vapply(f, stats::BIC, 1)))
  expect_identical(cmp$df, c(8L, 9L, 15L, 17L))
  expect_true(all(cmp$observed_zeros == d$observed_zeros))
  expect_true(all(attr(f, "checks")$converged))
})

test_that("the count tables for the rootograms agree with the comparison", {
  s <- nmes_env()
  f <- nmes_cached_fits(s)
  cmp <- s$nmes_comparison(f)
  for (m in c("poisson", "nb", "zinb")) {
    tab <- s$nmes_count_table(f[[m]], 30L)
    expect_identical(nrow(tab), 32L, info = m)
    expect_equal(sum(tab$observed), length(f$nb$y), info = m)
    # the expected frequencies cover all the probability, and the first is
    # the comparison table's expected zeros
    expect_equal(sum(tab$expected), length(f$nb$y), tolerance = 1e-6,
                 info = m)
    expect_equal(tab$expected[1], cmp$expected_zeros[cmp$model == m],
                 info = m)
  }
  # the Poisson table by hand at one count
  expect_equal(s$nmes_count_table(f$poisson)$expected[4],
               sum(stats::dpois(3, stats::fitted(f$poisson))))
  # the panel draws without error (to a null device)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(s$nmes_rootogram_panel(s$nmes_count_table(f$nb), "NB"))
})

test_that("effect tables are exp() of the estimates and Wald limits", {
  s <- nmes_env()
  f <- nmes_cached_fits(s)
  irr <- s$nmes_effects(f$zinb, "count")
  or <- s$nmes_effects(f$zinb, "zero")
  expect_identical(irr$term, c("hospital", "healthpoor", "healthexcellent",
                               "chronic", "gendermale", "school",
                               "insuranceyes"))
  expect_identical(or$term, c("hospital", "chronic", "gendermale", "school",
                              "insuranceyes"))
  ci <- unclass(confint(f$zinb, method = "wald"))
  cf <- coef(f$zinb)
  expect_false("(Intercept)" %in% c(irr$term, or$term))
  expect_equal(or$ratio, unname(exp(cf[paste0("zero_", or$term)])))
  expect_equal(or$lower, unname(exp(ci[paste0("zero_", or$term), 1])))
  expect_equal(irr$upper, unname(exp(ci[paste0("count_", irr$term), 2])))
  sc <- summary(f$zinb)$coefficients$zero
  expect_equal(or$p_value, unname(sc[or$term, "Pr(>|z|)"]))
  nb <- s$nmes_effects(f$nb, "count")
  expect_equal(nb$ratio, unname(exp(coef(f$nb)[nb$term])))
  expect_error(s$nmes_effects(f$nb, "zero"), "no zero part")
  th <- s$nmes_theta(f$zinb)
  expect_equal(th[["theta"]], f$zinb$theta)
  expect_lt(th[["profile_lower"]], th[["theta"]])
  expect_gt(th[["profile_upper"]], th[["theta"]])
})

# ---- the checkpointed bootstrap ------------------------------------------------

test_that("chunks and a resume reproduce one uninterrupted vuong_test() call", {
  skip_on_cran()
  s <- nmes_env()
  bf <- boot_fixture()
  single <- suppressWarnings(
    vuong_test(bf$zinb, bf$nb, B = 6, seed = 651))$bootstrap
  td <- tempfile("nmes")
  od <- file.path(td, "boot")
  on.exit(unlink(td, recursive = TRUE))
  # run two of the three chunks, then stop as if interrupted
  r1 <- s$run_nmes_bootstrap(B = 6, chunk = 2, out_dir = od, fits = bf,
                             prov = fake_prov(), max_chunks = 2,
                             progress = FALSE)
  expect_false(r1$complete)
  expect_identical(r1$chunks_done, 2L)
  ck_file <- s$nmes_checkpoint_file(od)
  before <- readRDS(ck_file)
  # a half-written checkpoint and a disturbed random stream must not matter
  writeBin(as.raw(1:10), paste0(ck_file, ".tmp"))
  set.seed(99)
  invisible(stats::runif(5))
  r2 <- s$run_nmes_bootstrap(B = 6, chunk = 2, out_dir = od, fits = bf,
                             prov = fake_prov(), progress = FALSE)
  after <- readRDS(ck_file)
  expect_true(r2$complete)
  expect_false(file.exists(paste0(ck_file, ".tmp")))
  expect_identical(nrow(r2$resumes), 1L)
  expect_identical(r2$resumes$at_chunk, 3L)
  # the finished chunks were kept exactly as they were
  expect_identical(after$chunks[1:2], before$chunks)
  # the same data sets in the same order as one call with the same seed
  expect_identical(r2$replicates$lr_star, single$replicates)
  expect_identical(as.vector(r2$status), as.vector(single$status))
  expect_equal(r2$observed$statistic, single$statistic)
  # the pooled p-values use the package's formulas
  p <- r2$pool
  for (nm in c("p.value", "p.bounds", "used", "failures", "exceed",
               "conclusive", "reasons", "reject")) {
    expect_identical(p[[nm]], single[[nm]], info = nm)
  }
})

test_that("a checkpoint is only resumed by the identical run", {
  skip_on_cran()
  s <- nmes_env()
  bf <- boot_fixture()
  td <- tempfile("nmes")
  od <- file.path(td, "boot")
  on.exit(unlink(td, recursive = TRUE))
  s$run_nmes_bootstrap(B = 6, chunk = 2, out_dir = od, fits = bf,
                       prov = fake_prov(), max_chunks = 1, progress = FALSE)
  run <- function(...) {
    args <- utils::modifyList(list(B = 6, chunk = 2, out_dir = od, fits = bf,
                                   prov = fake_prov(), progress = FALSE),
                              list(...))
    do.call(s$run_nmes_bootstrap, args)
  }
  expect_error(run(B = 8), "different B")
  expect_error(run(chunk = 3), "different chunk")
  expect_error(run(seed = 652), "different seed")
  expect_error(run(maxit = 20000), "different maxit")
  expect_error(run(prov = fake_prov(fingerprint = "other")),
               "different fingerprint_installed")
  bf2 <- bf
  bf2$nb <- nbreg(y ~ 1, data.frame(y = bf$nb$y))
  expect_error(run(fits = bf2), "different fit_hash")
  expect_error(run(resume = FALSE), "already exists")
  # nothing was added by the refused attempts
  expect_length(readRDS(s$nmes_checkpoint_file(od))$chunks, 1L)
})

# Build results from a doctored checkpoint: `lr` are the replicate values
# (NA = failed), split into chunks of `chunk`.
doctored <- function(s, lr, obs = 10, chunk = 100L, B = length(lr),
                     prov = fake_prov(), drop_last = FALSE) {
  idx <- split(seq_along(lr), ceiling(seq_along(lr) / chunk))
  chunks <- lapply(seq_along(idx), function(i) {
    r <- lr[idx[[i]]]
    list(replicates = r, observed = obs, loglik = c(nb = -1, zinb = -1),
         status = table(factor(ifelse(is.na(r), "failure", "interior"),
                               levels = c("interior", "boundary",
                                          "boundary_accepted", "failure"))),
         monotone = 0L, first_failure = if (anyNA(r)) "x" else NA_character_,
         warnings = character(), seconds = 1)
  })
  if (drop_last) chunks <- chunks[-length(chunks)]
  ck <- list(header = list(B = as.integer(B), chunk = as.integer(chunk),
                           seed = 651L, alpha = 0.05, maxit = 5000L,
                           rng_kind = RNGkind()),
             provenance = prov, chunks = chunks,
             resumes = data.frame(time = character(), at_chunk = integer()),
             started = "now")
  s$nmes_boot_results(ck)
}

test_that("run status and inference are separate verdicts", {
  s <- nmes_env()
  # complete and conclusive: 200 usable replicates, none at or above LR
  r <- doctored(s, rep(1, 200))
  st <- s$nmes_bootstrap_status(r)
  expect_identical(st$run_status, "COMPLETE")
  expect_identical(st$inference, "CONCLUSIVE")
  expect_identical(r$pool$p.value, 1 / 201)
  expect_true(r$pool$reject)
  # complete but inconclusive: more than 10% failed
  r <- doctored(s, c(rep(1, 170), rep(NA, 30)))
  st <- s$nmes_bootstrap_status(r)
  expect_identical(st$run_status, "COMPLETE")
  expect_identical(st$inference, "INCONCLUSIVE")
  expect_match(st$inference_reasons, "more than 10%", all = FALSE)
  expect_true(is.na(r$pool$reject))
  # complete but inconclusive: failures could change the decision
  r <- doctored(s, c(rep(1, 185), rep(20, 5), rep(NA, 10)))
  st <- s$nmes_bootstrap_status(r)
  expect_identical(st$inference, "INCONCLUSIVE")
  expect_match(st$inference_reasons, "could change the decision", all = FALSE)
  expect_equal(unname(r$pool$p.bounds), c(6 / 201, 16 / 201))
  # complete but inconclusive: too few usable replicates
  st <- s$nmes_bootstrap_status(doctored(s, rep(1, 50), chunk = 25L))
  expect_identical(st$inference, "INCONCLUSIVE")
  expect_match(st$inference_reasons, "at least 99", all = FALSE)
  # missing chunks: the run failed, and the inference is not evaluated
  r <- doctored(s, rep(1, 200), drop_last = TRUE)
  expect_false(r$complete)
  expect_null(r$pool)
  st <- s$nmes_bootstrap_status(r)
  expect_identical(st$run_status, "FAILED")
  expect_identical(st$inference, "NOT EVALUATED")
  # unverified provenance: failed
  st <- s$nmes_bootstrap_status(doctored(s, rep(1, 200),
                                         prov = fake_prov("unverified")))
  expect_identical(st$run_status, "FAILED")
})

test_that("the maxit rule uses failure shares only", {
  s <- nmes_env()
  pilot <- function(n_fail, maxit, lr = 1) {
    r <- doctored(s, c(rep(lr, 50 - n_fail), rep(NA, n_fail)), chunk = 25L)
    r$settings$maxit <- as.integer(maxit)
    r
  }
  d <- s$nmes_maxit_decision(pilot(2, 5000))
  expect_identical(d$maxit, 5000L)
  expect_false(d$needs_high_pilot)
  d <- s$nmes_maxit_decision(pilot(5, 5000))
  expect_true(d$needs_high_pilot)
  expect_identical(s$nmes_maxit_decision(pilot(5, 5000),
                                         pilot(1, 20000))$maxit, 20000L)
  expect_identical(s$nmes_maxit_decision(pilot(5, 5000),
                                         pilot(5, 20000))$maxit, 5000L)
  # the replicate values (and so any p-value) do not enter the rule
  expect_identical(
    s$nmes_maxit_decision(pilot(5, 5000, lr = 0), pilot(1, 20000, lr = 99)),
    s$nmes_maxit_decision(pilot(5, 5000, lr = 99), pilot(1, 20000, lr = 0)))
})

test_that("saving writes the two status lines first", {
  s <- nmes_env()
  td <- tempfile("nmes")
  on.exit(unlink(td, recursive = TRUE))
  s$save_nmes_bootstrap(doctored(s, c(rep(1, 170), rep(NA, 30))), td)
  st <- readLines(file.path(td, "nmes_bootstrap_status.txt"))
  expect_identical(st[1], "COMPLETE")
  expect_match(st[2], "^INCONCLUSIVE: ")
  expect_true(file.exists(file.path(td, "nmes_bootstrap_results.rds")))
  csv <- utils::read.csv(file.path(td, "nmes_bootstrap_replicates.csv"))
  expect_identical(nrow(csv), 200L)
  s$save_nmes_bootstrap(doctored(s, rep(1, 200), drop_last = TRUE), td)
  st <- readLines(file.path(td, "nmes_bootstrap_status.txt"))
  expect_identical(st[1], "FAILED")
  expect_identical(st[3], "NOT EVALUATED")
})

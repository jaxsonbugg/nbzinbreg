# Investigation of the simulation study's results of record (not part of the
# package; run after inst/benchmark/simulation.R, with the package installed).
#
# The results of record flagged biases above 3 Monte Carlo SEs for the zero
# part and log(theta). This script checks (1) whether pscl::zeroinfl() shows
# the same on the same simulated data, recording both packages' log-likelihoods
# (it also finds data sets where the two packages end at different maxima), and
# (2) whether the biases shrink at n = 2000, as finite-sample bias of maximum
# likelihood should. Findings: dev/progress_log.md, simulation phase.
#
# Run from the repository root: Rscript dev/simulation_bias_investigation.R

suppressPackageStartupMessages({library(nbzinbreg); library(pscl)})
source(system.file("benchmark", "simulation_functions.R", package = "nbzinbreg"))
saved <- readRDS("inst/benchmark/simulation_results.rds")
tr <- sim_truth(); pn <- names(tr$params)[1:7]
des <- saved$design
out <- file.path(tempdir(), "")

# Part 1: both packages on data sets 1-300, with log-likelihoods
K <- 300
rows <- lapply(seq_len(K), function(r) {
  d <- data.frame(y = sim_response(des, tr, 651 + r)$y, des)
  f <- zinbreg(y ~ x1 + x2 | x1 + z, data = d)
  g <- zeroinfl(y ~ x1 + x2 | x1 + z, data = d, dist = "negbin",
                control = zeroinfl.control(reltol = 1e-12))
  e1 <- c(coef(f), log(f$theta)); e2 <- c(coef(g), log(g$theta))
  data.frame(r = r, maxdiff = max(abs(e1 - e2)), ll_ours = f$loglik, ll_zi = g$loglik,
             conv_ours = f$converged, grad_ok = isTRUE(f$gradient_check$ok),
             zi_code = g$optim$convergence, theta_ours = f$theta, theta_zi = g$theta,
             g0_ours = coef(f)[["zero_(Intercept)"]], g0_zi = coef(g)[["zero_(Intercept)"]])
})
p1 <- do.call(rbind, rows)
saveRDS(p1, paste0(out, "invest_part1.rds"))
dis <- p1[p1$maxdiff > 1e-3, ]
cat("data sets with max |diff| > 1e-3:", nrow(dis), "of", K, "\n")
dis$ll_diff <- dis$ll_ours - dis$ll_zi
print(dis, digits = 6, row.names = FALSE)
cat("agreeing data sets: max diff", signif(max(p1$maxdiff[p1$maxdiff <= 1e-3]), 2), "\n")

# Part 2: n = 2000, R = 300
des2 <- sim_design(2000, 651)
R2 <- 300
est2 <- t(vapply(seq_len(R2), function(r) {
  y <- sim_response(des2, tr, 651 + r)$y
  f <- zinbreg(y ~ x1 + x2 | x1 + z, data = data.frame(y = y, des2))
  c(coef(f), `log(theta)` = log(f$theta))
}, numeric(7)))
saveRDS(est2, paste0(out, "invest_part2.rds"))
b500 <- saved$summary$bias[1:7]; se500 <- saved$summary$bias_mcse[1:7]
b2000 <- unname(colMeans(est2) - tr$params[pn]); se2000 <- unname(apply(est2, 2, sd) / sqrt(R2))
tab <- data.frame(parameter = pn, bias_n500 = b500, mcse_n500 = se500,
                  bias_n2000 = b2000, mcse_n2000 = se2000, ratio = b500 / b2000,
                  sd_n500 = saved$summary$mc_sd[1:7], sd_n2000 = unname(apply(est2, 2, sd)))
print(format(tab, digits = 3), row.names = FALSE)

# Benchmark: nbzinbreg against the reference implementations.
#
# Status: header only. This script is written once nbreg() and zinbreg() are
# implemented.
#
# Plan:
# - Data sets: pscl::bioChemists and AER::NMES1988.
# - For each data set, fit NB with nbreg() and MASS::glm.nb(), and ZINB with
#   zinbreg() and pscl::zeroinfl(dist = "negbin"). The reference fits use
#   tight tolerances: glm.control(epsilon = 1e-10, maxit = 100) and
#   zeroinfl.control(reltol = 1e-12).
# - Report the two sets of estimates side by side with their maximum absolute
#   difference, for coefficients, standard errors, theta (NB) or log(theta)
#   (ZINB), and the log-likelihood. Required agreement: 1e-5 / 1e-5 / 1e-4 /
#   1e-4 for NB, and 1e-3 for every ZINB quantity.
# - Time each fit (for example with bench::mark()).
# - Save the results next to this script so the report can read them.
#
# The reference functions are used here on purpose: this script is not part
# of the fitting code, which never calls them.

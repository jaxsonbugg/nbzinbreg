# Simulation study: does zinbreg() recover known parameters?
#
# Status: header only. This script is written once zinbreg() is implemented.
#
# Plan:
# - Set a seed so the results are reproducible.
# - Simulate R = 100 data sets of size n = 500 from a ZINB model with one or
#   two covariates in each part and known (beta, gamma, theta).
# - Fit each data set with zinbreg().
# - Report, for each parameter: the bias, the mean of the estimated standard
#   errors compared with the Monte Carlo standard deviation of the estimates,
#   and the empirical coverage of 95% Wald intervals.
# - Save the results next to this script so the report can read them.

# Transitional tests for documented stubs.
#
# Every function in `stubbed` is still a stub: its body calls the internal
# not_implemented(), which signals an error of class
# "nbzinbreg_not_implemented".
#
# When a function is implemented:
#   1. remove its name from `stubbed` (and its entry from `stub_calls`);
#   2. add its real tests in its own test file.
# Nothing else in this file needs to change. The second test fails if the
# list and the code disagree, so a forgotten update cannot slip through.

stubbed <- c(
  # required fitting interface
  "nbreg", "zinbreg", "zinbreg_control",
  # required methods, nbreg
  "print.nbreg", "summary.nbreg", "print.summary.nbreg", "coef.nbreg",
  "vcov.nbreg", "logLik.nbreg", "fitted.nbreg", "residuals.nbreg",
  "nobs.nbreg", "predict.nbreg",
  # required methods, zinbreg
  "print.zinbreg", "summary.zinbreg", "print.summary.zinbreg",
  "coef.zinbreg", "vcov.zinbreg", "logLik.zinbreg", "fitted.zinbreg",
  "residuals.zinbreg", "nobs.zinbreg", "predict.zinbreg",
  # provisional extras
  "confint.nbreg", "confint.zinbreg", "anova.nbreg", "overdispersion_test",
  "vuong_test", "simulate.nbreg", "simulate.zinbreg", "boot_se",
  "expected_counts", "rootogram",
  "tidy.nbreg", "glance.nbreg", "tidy.zinbreg", "glance.zinbreg"
)

# A minimal call for each stub. Methods are called through their generic on a
# placeholder object, which also checks that the S3 method is registered.
nb <- structure(list(), class = "nbreg")
zinb <- structure(list(), class = "zinbreg")
stub_calls <- list(
  nbreg = function() nbreg(),
  zinbreg = function() zinbreg(),
  zinbreg_control = function() zinbreg_control(),
  print.nbreg = function() print(nb),
  summary.nbreg = function() summary(nb),
  print.summary.nbreg = function() print(structure(list(), class = "summary.nbreg")),
  coef.nbreg = function() coef(nb),
  vcov.nbreg = function() vcov(nb),
  logLik.nbreg = function() logLik(nb),
  fitted.nbreg = function() fitted(nb),
  residuals.nbreg = function() residuals(nb),
  nobs.nbreg = function() nobs(nb),
  predict.nbreg = function() predict(nb),
  print.zinbreg = function() print(zinb),
  summary.zinbreg = function() summary(zinb),
  print.summary.zinbreg = function() print(structure(list(), class = "summary.zinbreg")),
  coef.zinbreg = function() coef(zinb),
  vcov.zinbreg = function() vcov(zinb),
  logLik.zinbreg = function() logLik(zinb),
  fitted.zinbreg = function() fitted(zinb),
  residuals.zinbreg = function() residuals(zinb),
  nobs.zinbreg = function() nobs(zinb),
  predict.zinbreg = function() predict(zinb),
  confint.nbreg = function() confint(nb),
  confint.zinbreg = function() confint(zinb),
  anova.nbreg = function() anova(nb),
  overdispersion_test = function() overdispersion_test(nb),
  vuong_test = function() vuong_test(zinb, nb),
  simulate.nbreg = function() simulate(nb),
  simulate.zinbreg = function() simulate(zinb),
  boot_se = function() boot_se(nb),
  expected_counts = function() expected_counts(nb),
  rootogram = function() rootogram(nb),
  tidy.nbreg = function() tidy(nb),
  glance.nbreg = function() glance(nb),
  tidy.zinbreg = function() tidy(zinb),
  glance.zinbreg = function() glance(zinb)
)

test_that("every listed stub signals the not-implemented condition", {
  expect_setequal(names(stub_calls), stubbed)
  for (nm in stubbed) {
    expect_error(stub_calls[[nm]](), class = "nbzinbreg_not_implemented",
                 info = nm)
  }
})

test_that("the stub list matches the functions that are still stubs", {
  skip_if_not_installed("codetools")
  fns <- package_functions()
  fns$not_implemented <- NULL
  still_stub <- vapply(fns, function(f) {
    "not_implemented" %in% codetools::findGlobals(f, merge = FALSE)$functions
  }, logical(1))
  # Fails if a stub is missing from the list, or if an implemented function
  # was left in it.
  expect_setequal(names(fns)[still_stub], stubbed)
})

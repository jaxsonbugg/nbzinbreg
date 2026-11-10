# Guard for documented stubs.
#
# Every function in `stubbed` is still a stub: its body calls the internal
# not_implemented(), which signals an error of class
# "nbzinbreg_not_implemented". Every planned function is now implemented, so
# the list is empty; the second test fails if a stub is (re)introduced
# without being listed here.
#
# When a function is implemented:
#   1. remove its name from `stubbed` (and its entry from `stub_calls`);
#   2. add its real tests in its own test file.

stubbed <- character()

# A minimal call for each stub. Methods are called through their generic on a
# placeholder object, which also checks that the S3 method is registered.
nb <- structure(list(), class = "nbreg")
zinb <- structure(list(), class = "zinbreg")
stub_calls <- stats::setNames(list(), character())

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

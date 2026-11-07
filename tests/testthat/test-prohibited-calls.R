# The detector itself lives in helper-prohibited.R. These tests first check
# that the detector catches every prohibited pattern and ignores harmless
# ones, then run it over every function the package defines.

test_that("detector flags direct, qualified and indirect prohibited calls", {
  skip_if_not_installed("codetools")

  flagged <- list(
    direct_call      = function(x, f) optim(x, f),
    passed_as_value  = function(xs) lapply(xs, optimize),
    assigned         = function() { g <- uniroot; g },
    qualified        = function(d) stats::glm(y ~ x, data = d),
    triple_colon     = function(x, y) stats:::lm.fit(x, y),
    bad_package      = function(d) pscl::hurdle(y ~ x, data = d),
    other_fitter     = function(d) MASS::rlm(y ~ x, data = d),
    do_call_string   = function(args) do.call("glm", args),
    match_fun        = function() match.fun("optim"),
    get_string       = function() get("nlminb"),
    exported_value   = function() getExportedValue("MASS", "glm.nb"),
    parse_text       = function(f) eval(parse(text = "optim(1, f)")),
    str2lang_text    = function() str2lang("glm(y ~ x)"),
    in_default_arg   = function(fit = optim) fit,
    nested_function  = function() function(x) nlm(x, 1)
  )
  for (nm in names(flagged)) {
    expect_true(length(detect_prohibited(flagged[[nm]])) > 0,
                info = paste("should be flagged:", nm))
  }
})

test_that("detector ignores comments, local names and ordinary strings", {
  skip_if_not_installed("codetools")

  clean <- list(
    comment_mention = function(x) {
      # same parameterization as MASS::glm.nb(); we never call optim() or glm()
      x + 1
    },
    local_variable  = function() { lm <- 1; lm + 1 },
    class_check     = function(x) inherits(x, "glm"),
    argument_tag    = function() list(glm = 1, optim = 2),
    message_text    = function() stop("this package does not use optim() or glm()"),
    allowed_tools   = function(X, w, z) qr.coef(qr(sqrt(w) * X), sqrt(w) * z),
    allowed_derivs  = function(y, th) digamma(y + th) - digamma(th)
  )
  for (nm in names(clean)) {
    expect_identical(detect_prohibited(clean[[nm]]), character(),
                     info = paste("should not be flagged:", nm))
  }
})

test_that("no function defined by nbzinbreg uses a prohibited function", {
  skip_if_not_installed("codetools")

  fns <- package_functions()
  # Guard against a scoping mistake that would silently check nothing.
  expect_gt(length(fns), 0)

  hits <- Filter(length, lapply(fns, detect_prohibited))
  if (length(hits)) {
    report <- vapply(names(hits), function(nm) {
      paste0(nm, ": ", paste(hits[[nm]], collapse = ", "))
    }, character(1))
    fail(paste(c("Prohibited calls found:", report), collapse = "\n  "))
  } else {
    succeed()
  }
})

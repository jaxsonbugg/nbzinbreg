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
    nested_function  = function() function(x) nlm(x, 1),
    # names built at run time
    get_pasted       = function() get(paste0("op", "tim")),
    get_named_arg    = function() get(x = paste("nl", "minb", sep = "")),
    match_fun_built  = function(s) match.fun(sprintf("%s", s)),
    do_call_built    = function(a) do.call(paste0("gl", "m"), a),
    exported_built   = function() getExportedValue("stats", paste0("op", "tim")),
    # functions taken out of a namespace by name
    namespace_dollar = function() asNamespace("stats")$optim,
    namespace_brackets = function() asNamespace("stats")[["optim"]],
    get_namespace    = function(d) getNamespace("stats")$glm(y ~ x, data = d),
    environment_get  = function() as.environment("package:stats")[["nlminb"]]
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
    allowed_derivs  = function(y, th) digamma(y + th) - digamma(th),
    # passing a function held in a variable, or a qualified one, is fine
    do_call_symbol  = function(f, a) do.call(f, a),
    do_call_qualified = function(a) do.call(graphics::plot, a),
    do_call_literal = function(a) do.call("rbind", a),
    # a field that happens to share a prohibited name is not a lookup
    list_field      = function(fit) fit$glm,
    namespace_ok    = function() asNamespace("stats")$qlogis
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
  # Internal engines and S3 methods are checked; functions re-exported or
  # imported from other packages (generics::tidy, generics::glance,
  # Formula::Formula, numDeriv::hessian) are not, because they are not
  # the package's own code.
  expect_true(all(c("nbreg_fit", "zinb_em", "logit_irls", "nbreg",
                    "zinbreg", "predict.zinbreg", "tidy.zinbreg") %in%
                    names(fns)))
  expect_false(any(c("tidy", "glance", "Formula", "hessian") %in% names(fns)))

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

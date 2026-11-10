# Detector for prohibited function calls in the package's own code.
#
# The project rule: the fitting code must implement the estimation itself, so
# no package function may call a general-purpose optimizer or a ready-made
# GLM/NB/ZINB fitter. Tests, benchmarks and the report may use those reference
# functions; the package code may not.
#
# How it works, in plain terms:
# 1. Only functions that the package itself defines are checked. A function
#    counts as package-defined when its enclosing top-level environment is the
#    nbzinbreg namespace. Re-exported or imported functions (for example
#    generics::tidy) live in other namespaces, so they are skipped.
# 2. Each function's parsed body is examined, never its source text. Parsing
#    throws comments away, so a comment such as "matches glm.nb()" can never
#    trigger the detector, and neither can roxygen documentation.
# 3. Three kinds of use are flagged:
#    - direct use of a prohibited name, either called (`optim(...)`) or passed
#      as a value (`lapply(x, optim)`, `f <- optim`). codetools::findGlobals()
#      lists the global names a function uses and leaves out local variables,
#      so a local variable that happens to be called `lm` is not flagged;
#    - qualified use, `pkg::fn` or `pkg:::fn`, where either the function name
#      or the package is prohibited;
#    - indirect use by name, where a prohibited name is given as a string to a
#      function that looks functions up or builds calls (`do.call("glm", ...)`,
#      `match.fun("optim")`, `getExportedValue("MASS", "glm.nb")`, ...), or
#      inside code text passed to `parse(text = )`, `str2lang()` or
#      `str2expression()`.
#    Ordinary strings, such as error messages or `inherits(x, "glm")`, are not
#    flagged, because they do not call anything.
# 4. Two ways of hiding the name are flagged as well:
#    - a name built at run time and handed to a lookup function, for example
#      `get(paste0("op", "tim"))` or `do.call(sprintf("%s", f), args)`. The
#      detector cannot know what such a string will be, so any computed
#      (non-literal) function name given to `get()`, `get0()`, `mget()`,
#      `match.fun()`, `do.call()`, `getExportedValue()` or
#      `getFromNamespace()` is flagged. A plain variable
#      (`do.call(f, args)`) or a literal is not computed and stays allowed;
#    - a function taken out of a namespace by name, for example
#      `asNamespace("stats")$optim`, `asNamespace("stats")[["optim"]]` or
#      `getNamespace("stats")$glm`, where the object is a namespace or
#      environment lookup and the selected name is prohibited.
#    What remains possible (for example building a call with `as.call()` from
#    pieces) is beyond a static check; the rule in CLAUDE.md still applies,
#    and R/ is also checked by reading it.

prohibited_functions <- c(
  "optim", "optimize", "optimise", "nlm", "nlminb", "uniroot",
  "glm", "glm.fit", "lm", "lm.fit", "lm.wfit",
  "glm.nb", "theta.ml", "zeroinfl",
  "mle", "mle2", "optimx", "maxLik"
)

# Any qualified use of these packages inside R/ is flagged. MASS and pscl are
# only reference packages for tests, so the package code has no reason to
# call them at all.
prohibited_packages <- c(
  "countreg", "glmmTMB", "gamlss", "VGAM", "bbmle", "optimx", "maxLik",
  "stats4", "MASS", "pscl"
)

# Functions that look up another function, or build a call, from a string.
lookup_functions <- c(
  "do.call", "match.fun", "get", "get0", "mget", "getExportedValue",
  "getFromNamespace", "call", "as.name", "as.symbol"
)

# Functions that turn a string of R code into code.
parse_functions <- c("parse", "str2lang", "str2expression")

# For the lookup functions that take a function name: which argument holds
# it (by name, and by position when unnamed).
name_argument <- list(
  get = c("x", 1L), get0 = c("x", 1L), mget = c("x", 1L),
  match.fun = c("FUN", 1L), do.call = c("what", 1L),
  getExportedValue = c("name", 2L), getFromNamespace = c("x", 1L))

# Calls that return a namespace or an environment to look a name up in.
namespace_functions <- c("asNamespace", "getNamespace", "loadNamespace",
                         "as.environment", "environment", "baseenv",
                         "globalenv", "topenv", "parent.env")

# The argument of a call that holds the function name, or NULL.
name_arg_of <- function(x, head_name) {
  spec <- name_argument[[head_name]]
  args <- as.list(x)[-1]
  nms <- names(args)
  if (!is.null(nms) && spec[1] %in% nms) return(args[[spec[1]]])
  unnamed <- if (is.null(nms)) args else args[!nzchar(nms)]
  pos <- as.integer(spec[2])
  if (length(unnamed) >= pos) unnamed[[pos]] else NULL
}

# TRUE if `a` is a computed name: a call other than `::` / `:::` (a
# qualified function is a function, not a name to look up).
is_computed_name <- function(a) {
  is.call(a) && !(is.symbol(a[[1]]) &&
                    as.character(a[[1]]) %in% c("::", ":::", "function"))
}

is_empty_arg <- function(x) identical(x, quote(expr = ))

# Walk an expression and collect qualified and string-based prohibited uses.
# `from_text = TRUE` is used for code recovered from strings, where
# findGlobals() cannot run, so direct calls are flagged by name there.
walk_prohibited <- function(expr, from_text = FALSE) {
  hits <- character()

  walk <- function(x) {
    if (is.pairlist(x)) {
      # formal arguments of an inner function definition, with defaults
      for (i in seq_along(x)) if (!is_empty_arg(x[[i]])) walk(x[[i]])
      return(invisible())
    }
    if (!is.call(x)) return(invisible())

    head <- x[[1]]
    head_name <- if (is.symbol(head)) as.character(head) else ""

    if (head_name %in% c("::", ":::")) {
      pkg <- as.character(x[[2]])
      fn <- as.character(x[[3]])
      if (pkg %in% prohibited_packages || fn %in% prohibited_functions) {
        hits <<- c(hits, paste0(pkg, head_name, fn))
      }
      return(invisible())
    }

    if (from_text && head_name %in% prohibited_functions) {
      hits <<- c(hits, paste0(head_name, "() in parsed text"))
    }

    args <- as.list(x)[-1]
    strings <- unlist(Filter(is.character, args), use.names = FALSE)

    if (head_name %in% lookup_functions && length(strings)) {
      bad <- strings[strings %in% c(prohibited_functions, prohibited_packages)]
      if (length(bad)) {
        hits <<- c(hits, paste0(head_name, '("', bad, '")'))
      }
    }

    # a function name computed at run time and handed to a lookup function
    if (head_name %in% names(name_argument)) {
      a <- name_arg_of(x, head_name)
      if (!is.null(a) && is_computed_name(a)) {
        hits <<- c(hits, paste0(head_name, "(<computed name: ",
                                paste(deparse(a), collapse = " "), ">)"))
      }
    }

    # a prohibited function taken out of a namespace or environment by name
    if (head_name %in% c("$", "[[", "@") && length(x) >= 3L) {
      obj <- x[[2]]
      sel <- x[[3]]
      sel_name <- if (is.symbol(sel) || is.character(sel)) as.character(sel) else ""
      obj_head <- if (is.call(obj) && is.symbol(obj[[1]])) as.character(obj[[1]]) else ""
      if (obj_head %in% namespace_functions &&
          sel_name %in% prohibited_functions) {
        hits <<- c(hits, paste0(obj_head, "(...)", head_name, sel_name))
      }
    }

    if (head_name %in% parse_functions && length(strings)) {
      for (s in strings) {
        parsed <- tryCatch(parse(text = s, keep.source = FALSE),
                           error = function(e) NULL)
        for (p in parsed) {
          hits <<- c(hits, walk_prohibited(p, from_text = TRUE))
        }
      }
    }

    if (is.call(head)) walk(head)
    for (i in seq_along(args)) {
      if (!is_empty_arg(args[[i]])) walk(args[[i]])
    }
    invisible()
  }

  walk(expr)
  unique(hits)
}

# All prohibited uses in one function, as a character vector (empty if clean).
detect_prohibited <- function(f) {
  globals <- codetools::findGlobals(f, merge = FALSE)
  direct <- intersect(c(globals$functions, globals$variables),
                      prohibited_functions)
  formals_hits <- unlist(lapply(formals(f), function(d) {
    if (is_empty_arg(d)) character() else walk_prohibited(d)
  }), use.names = FALSE)
  unique(c(direct, walk_prohibited(body(f)), formals_hits))
}

# The functions defined by the package itself, keyed by name.
package_functions <- function(pkg = "nbzinbreg") {
  ns <- asNamespace(pkg)
  objs <- mget(ls(ns, all.names = TRUE), envir = ns)
  keep <- vapply(objs, function(f) {
    is.function(f) && !is.primitive(f) &&
      identical(topenv(environment(f)), ns)
  }, logical(1))
  objs[keep]
}

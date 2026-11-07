# Internal helpers shared across the package.

# Signal that an exported function is still a documented stub.
#
# Every planned function exists with its final (or provisional) signature and
# documentation, but until it is implemented its body calls this helper. The
# error has its own class, `nbzinbreg_not_implemented`, so tests can tell
# "not written yet" apart from any other error. When a function is
# implemented, its call to not_implemented() disappears and its name is
# removed from the `stubbed` list in tests/testthat/test-stubs.R.
not_implemented <- function(what) {
  cnd <- structure(
    class = c("nbzinbreg_not_implemented", "error", "condition"),
    list(
      message = sprintf(
        "`%s()` is not yet implemented in this development version of nbzinbreg.",
        what
      ),
      call = NULL
    )
  )
  stop(cnd)
}

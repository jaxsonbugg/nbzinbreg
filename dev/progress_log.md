# Progress log

Every project change is logged here with the date and the reason. Newest entries at the bottom.

## 2026-10-28: Repository scaffolding

**What:** Created the repository folder `C:\dev\nbzinbreg` with the layout described in `dev/requirements.md`, a `.gitignore` (secrets, R/Quarto artifacts, private files), `CLAUDE.md` (rules for Claude), a stub `README.md`, `LICENSE.md` (MIT), and the three notes in `dev/`. The course handout PDF is kept in `docs-private/`, which is git-ignored and never pushed.

**Why:** A clean, organized starting point so that the real work (the R package, tests, benchmarks, report) has a clear place to go from the first commit. The repository is public because GitHub Pages and `install_github()` need it to be.

**R environment:** R 4.3.1 was already installed at `C:\R-4.3.1` but not on the PATH, so `C:\R-4.3.1\bin` was added to the user PATH. Packages already present: usethis, roxygen2, testthat, MASS, pscl, AER, numDeriv, Formula, knitr, rmarkdown. Newly installed: pkgdown, bench, quarto, available, and devtools 2.4.5 (the newest build available as a Windows binary for R 4.3; the newest source version needs a newer `fs` that would require Rtools to compile). Some existing library DLLs (`rlang`, `cli`, `ps`) were locked by running R processes, so upgrading them in place failed (the old copies in `C:\R-4.3.1\library` were left intact and still load). Copies were installed into the per-user library `C:\Users\Buggb\AppData\Local\R\win-library\4.3` instead; they are the same versions, because the newest Windows binaries for R 4.3 are `rlang` 1.1.5, `cli` 3.6.4, `ps` 1.9.0, and `fs` 1.6.5.

**Package-name availability check for `nbzinbreg`:** the name is valid (letters only, starts with a letter) and no conflicts were found.

| Check | Result |
|---|---|
| CRAN package page (`cran.r-project.org/package=nbzinbreg`) | 404, not present |
| CRAN package index page (`/web/packages/nbzinbreg/index.html`) | 404, not present |
| CRAN archive (`/src/contrib/Archive/nbzinbreg/`) | 404, never published/archived |
| `available.packages()` from CRAN (25,181 packages), case-insensitive match | no match; no similar names (`nbzinb`, `zinbreg`, `nbreg`) |
| CRAN pages for near-names `nbzinb`, `zinbreg`, `nbreg`, `NBZINBREG`, `NBzinbreg` | all 404 |
| Bioconductor (`bioconductor.org/packages/nbzinbreg`) | redirects to the generic "removed packages" page, the same as for a made-up name, so not present |
| GitHub: `gh search repos nbzinbreg`, `github.com/jaxsonbugg/nbzinbreg`, `github.com/nbzinbreg` | no results / 404 |
| `available::available("nbzinbreg")` | valid name; available on CRAN, Bioconductor, and GitHub |

**GitHub:** the public repository `jaxsonbugg/nbzinbreg` was created and the first commit pushed to `main`. The first push was rejected because the GitHub CLI login token lacks the `workflow` scope, which GitHub requires for any file under `.github/workflows/`, even a placeholder. The placeholder was dropped from the commit. Before the CI workflow file (`R-CMD-check.yaml`) can be pushed, the token needs that scope, which is done with `gh auth refresh -h github.com -s workflow` (it opens a browser sign-in).

**Next:** generate the R package skeleton (DESCRIPTION, NAMESPACE, `.Rbuildignore`, license stub, CI workflow) with usethis, then write the README as the working status document.

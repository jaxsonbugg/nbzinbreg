# AI-use log

A running record of how Claude Code is used in this project. It is the raw material for the report's AI-use statement, which must say what Claude was asked to do, what was verified independently, and give at least one concrete example where Claude's output had to be corrected or rejected. Add an entry whenever Claude writes a non-trivial piece of code or the owner catches an error.

## Entry template

- **Date:**
- **Asked Claude to:**
- **What Claude produced:**
- **How it was checked** (tests, comparison with the reference function, reading the code against the equations):
- **Problems found and how they were corrected** (leave blank if none):

## Entries

### 2026-10-28: Repository scaffolding

- **Asked Claude to:** read the project description, propose a repository layout, create the folder structure, a `.gitignore`, `CLAUDE.md`, and starter notes, and verify the package name is free on CRAN.
- **What Claude produced:** the folder skeleton, `.gitignore`, `CLAUDE.md`, stub `README.md`, `LICENSE.md`, and the notes in `dev/`. It also ran the name-availability checks (results in `dev/progress_log.md`).
- **How it was checked:** the owner reviewed and approved the plan before any files were created.
- **Problems found and how they were corrected:** the first attempt to install a Python PDF library in a scratch folder was declined by the owner in favor of installing Poppler, which was then used to read the handout. A first installation attempt for `devtools` failed because of locked library files, which led to using a compatible binary version and the per-user library.

### 2026-11-07: Package skeleton

- **Asked Claude to:** plan, then generate the package skeleton with documented stubs, a prohibited-call detector, transitional stub tests, and the design document, without implementing any algorithm.
- **What Claude produced:** the files listed in `dev/progress_log.md` for this date.
- **How it was checked:** the owner reviewed the plan twice and required nine corrections before approving it. `devtools::test()`, `R CMD check`, `pkgdown::check_pkgdown()`, a local install, and `install_github()` in a fresh session were all run. The detector was tested on a deliberately bad function.
- **Problems found and how they were corrected:** in the owner's plan review, Claude's first plan worked around a missing GitHub permission with a partial push of all but the last commit, could have flagged functions imported from other packages as prohibited calls, treated every PDF as forbidden (but `report/report.pdf` is a required deliverable), and did not separate "skeleton verified" from "statistics verified"; all were fixed before any file was changed. During the work, one of Claude's detector self-tests referred to a package (`VGAM`) that is not a declared dependency, which `R CMD check` flagged. It was replaced with an equivalent test. Claude also first passed a Git Bash style path to R in the `install_github()` check, which R could not use; the check was rerun with a Windows path.

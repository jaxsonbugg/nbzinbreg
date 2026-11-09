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

### 2026-11-08: `nbreg()` implementation

- **Asked Claude to:** implement `nbreg()`, `nbreg_control()` and the required `nbreg` methods with exhaustive input checking, following the reviewed plan.
- **What Claude produced:** `R/nbreg-fit.R`, `R/nbreg.R`, `R/methods-nbreg.R`, the control validation and condition helpers, and their tests.
- **How it was checked:**
  - Every reference case was compared with `MASS::glm.nb()` at the required tolerances.
  - The θ derivatives were compared with `numDeriv`.
  - Every input error was triggered once by a probe script before tests were written.
  - The full `R CMD check` was run, plus a plain `R CMD check` without devtools, and CI on five platforms.
- **Problems found and how they were corrected:**
  - **In the owner's plan reviews (before any code was written), Claude's first plan:**
    - clipped η and floored μ, which would have silently changed the likelihood;
    - used the number of positive-weight observations instead of `sum(weights)` in the θ start;
    - treated a small fitted mean alone as evidence of divergence;
    - used an unjustified fixed step when the θ curvature was not negative.
  - **The owner's code review caught three problems after commits were pushed:**
    - `check_whole()` accepted values above the integer maximum;
    - `model.matrix()` errors were not classified;
    - edited `nbreg_control` objects were trusted without revalidation.
  - **Claude's own tests and probes caught:**
    - the near-integer response check rejecting −1e-10 as negative;
    - an unclassed offset error from `model.offset()`;
    - two wrong expectations in Claude's own tests (a print-format assumption, and an unexpected warning on a tiny data set).

### 2026-11-09: `zinbreg()` implementation

- **Asked Claude to:** plan, then implement, `zinbreg()`, `zinbreg_control()` and all required `zinbreg` methods with exhaustive input checking, following the reviewed plan. The work was pushed commit by commit for review.
- **What Claude produced:**
  - the control function and the shared control resolver;
  - the EM engine (`R/zinbreg-fit.R`);
  - the analytic Hessian and gradient check (`R/zinbreg-hessian.R`);
  - `zinbreg()` and the methods;
  - 640 new test expectations (498 before this phase, 1138 after).
- **How it was checked:**
  - every reference case was compared with `pscl::zeroinfl(reltol = 1e-12)` at the required tolerances, with both Hessian methods;
  - the analytic gradient and Hessian were compared with `numDeriv`;
  - every divergence and fallback scenario was first run as a probe script, before tests were written;
  - each regression test was confirmed to fail on the old code;
  - full `R CMD check`, plain `R CMD check`, and CI on five platforms after every pushed commit.
- **Problems found and how they were corrected:**
  - **In the owner's plan review (before any code was written),** six gaps in Claude's first ZINB plan:
    - non-converged M-steps could be treated as converged;
    - the logistic step was unstable when pi (1 − pi) underflows, and had no fallback start;
    - divergence detection was reused from the NB part instead of being specific to each part;
    - there was no gradient check at termination;
    - `type = "prob"` could allocate an enormous matrix;
    - a higher log-likelihood than the reference could have excused coefficient differences.
  - **The owner's code reviews caught, after pushed commits:**
    - the count M-step's rank-drop flag was not passed on;
    - a count part with every coefficient held fixed would have called the NB engine with a zero-column matrix;
    - coefficients the initial NB fit flagged were frozen, which could constrain a valid ZINB solution;
    - an unavailable gradient check could still count as convergence;
    - `print()` lacked the EM iteration count;
    - prediction standard errors could show 0 when no covariance matrix existed;
    - Pearson residuals were NaN when the fitted variance underflowed.

    Each was fixed in its own commit with regression tests.
  - **Claude's own probes and tests caught:**
    - the first EM-level drift check flagged the bioChemists zero intercept during ordinary slow convergence;
    - `tau − pi` cancelled to 0 when both rounded to 1;
    - a bad starting value left the count M-step stuck (log-likelihood −678.2 instead of −340.9), which led to the M-step restart;
    - a stray `"count_"` name came from `paste0()` on an empty vector;
    - the Hessian was symmetric only to 4e-15;
    - the monotonicity warning used a wrong column name;
    - several wrong expectations in Claude's own tests (an index, a count of zeros, and a tolerance tighter than EM's stopping rule).
  - **A judgment the owner had to make:** with weights, `zeroinfl()` disagreed with `zinbreg()` by 0.03. Claude did not relax the tolerance. It investigated and showed that zeroinfl's weighted fit stops at a point that is not a maximum. The owner agreed to use `zeroinfl()` on repeated rows as the reference for integer frequency weights.
  - **A plan value changed with evidence:** `grad_tol` became 1e-8 instead of 1e-6, because 1e-6 could allow errors above the 1e-3 tolerance on bioChemists. This is recorded in `dev/progress_log.md`.

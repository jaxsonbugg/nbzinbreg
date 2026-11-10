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

### 2026-11-09: The extras

- **Asked Claude to:** plan, then implement, the extras: `confint()` with profile intervals for theta, `anova()`, `overdispersion_test()`, `vuong_test()` with a bootstrap LR test, `simulate()`, `boot_se()`, `expected_counts()`, `rootogram()`, `tidy()` and `glance()`. The errors were to be comprehensive, and the work was done in eight commits, each pushed and reviewed before the next.
- **What Claude produced:**
  - `R/extras-utils.R`, `R/profile.R`, `R/anova.R`, `R/vuong.R`, and the extras in `R/inference.R`, `R/simulate.R`, `R/diagnostics.R` and `R/tidiers.R`;
  - fixed-theta EM and zero-part separation detection in the ZINB engine;
  - one test file per extra (about 650 new expectations).
- **How it was checked:**
  - every extra was compared with a reference: `confint.default(glm.nb)`, `confint(zeroinfl)`, `anova()` of `glm.nb` fits, zeroinfl log-likelihoods, `pscl::odTest()`, `pscl::vuong()` (also on repeated rows for weights), and `colSums(pscl::predprob())`;
  - profile limits were checked against profile log-likelihoods maximized independently (`glm(negative.binomial(theta))`, and `optim()` in tests only);
  - each regression test was confirmed to fail on the previous code;
  - the full suite with the slow tests, `devtools::check()`, plain `R CMD check`, and CI on five platforms after every commit.
- **Problems found and how they were corrected:**
  - **In the owner's plan review (before any code was written),** eight gaps in Claude's first plan:
    - the bootstrap LR test did not handle replicates at the pi = 0 boundary;
    - failed replicates had no effect on the p-value;
    - frequency weights were drawn once per row instead of independently;
    - profile-interval endpoints at a bound were not classified;
    - the overdispersion test did not separate the Poisson limit from a truncating bound;
    - `anova()` gave chi-squared p-values for irregular fits;
    - the weighted Vuong statistic was undefined;
    - several smaller points.
  - **The owner's code reviews caught, after pushed commits:**
    - `check_parm()` with `Inf`;
    - no allocation limit for bootstrap data;
    - weights that round to 0;
    - an integer overflow in the profile retry limit;
    - an unbounded endpoint certified outside the region;
    - a p-value for a non-converged NB fit;
    - unconverged boundary replicates counted as failures;
    - separation concluded from extreme probabilities alone;
    - `boot_se()` running on a non-converged fit;
    - an unconverged higher likelihood ignored in LR replicates.

    Each was fixed in its own commit with regression tests.
  - **Claude's own probes and tests caught:**
    - separation of the zero part by a continuous covariate, which the engine did not detect;
    - `boot_refit()` keeping a lower boundary result over a higher interior one;
    - an interior local maximum in the supplied ZINB fit (equal to zeroinfl's) on pure NB data, below a separated supremum;
    - a default-start refit crawling for thousands of EM iterations, which made the bootstrap 2.2 times slower;
    - several wrong expectations in Claude's own tests (a tolerance tighter than the stopping rule, an assumption that a boundary refit reaches the global maximum, `NULL` dropped from a list in a mock, and the plot region's 4% margin).

### 2026-11-10: The benchmark

- **Asked Claude to:** plan, then implement, the benchmark of handout Section 5 (both models, both data sets, side-by-side tables, the required tolerances, timing), and draft the report's Benchmarking section. The work was done in four commits, each pushed and reviewed.
- **What Claude produced:**
  - `inst/benchmark/benchmark_functions.R`: validity checks, comparisons, timing, the code-fingerprint provenance and saving;
  - `inst/benchmark/benchmark.R`;
  - `tests/testthat/test-benchmark.R`;
  - the results of record;
  - the report's Benchmarking section and the supplementary article.
- **How it was checked:**
  - every failure path was triggered in tests with doctored fits and tables;
  - a dry run showed that a failed run saves its results and exits with status 1;
  - the fingerprint was compared between source, `load_all()` and a real byte-compiled install;
  - the results of record were produced from a clean commit, with verified provenance;
  - the PDF was rendered and inspected page by page;
  - the test suite, both package checks and CI ran on every commit.
- **Problems found and how they were corrected:**
  - **In the owner's plan review,** five changes to Claude's first plan:
    - Claude proposed recording `git rev-parse HEAD`, which identifies the checkout, not the installed code; the owner required a verifiable check, which became the code fingerprint;
    - the first plan discarded all results when a tolerance failed;
    - it applied the theta tolerance to SE(theta), which the handout does not require;
    - it lacked explicit validity checks;
    - the report had no check that its saved timings came from the installed code.
  - **Claude's own work caught:**
    - a locale-dependent sort that made two identical fingerprints differ;
    - an unclosed file connection;
    - PDF tables wrapping negative numbers;
    - Quarto's cache files appearing in `vignettes/`;
    - an inaccurate sentence about when the Hessian is used.
  - **Claude also corrected its own claim** that the old pkgdown entry `anova.nbreg` was invalid: it resolves as an alias, so that commit fixed the index's content, not a failing check.

### 2026-11-10: The simulation study

- **Asked Claude to:** plan, then implement, the simulation study of handout Section 5.1. It was done in three commits, pushed as they came, with the report's Simulation study section.
- **What Claude produced:**
  - `inst/benchmark/simulation_functions.R`: the generator, the per-replication records and classes, the profile classes, the summaries with analytic and bootstrap Monte Carlo SEs, the status and saving;
  - `inst/benchmark/simulation.R`;
  - `tests/testthat/test-simulation.R`;
  - the results of record;
  - `dev/simulation_bias_investigation.R`;
  - the report section.
- **How it was checked:**
  - every summary and Monte Carlo SE by hand on synthetic records;
  - every replication and profile class from doctored inputs;
  - the generator against exact moments on a 200,000-row sample;
  - a pilot run before the full run;
  - the results of record from a clean commit with verified provenance;
  - the flagged biases investigated against `zeroinfl()` and at n = 2000 before committing;
  - the PDF inspected page by page;
  - the suite, both package checks and CI.
- **Problems found and how they were corrected:**
  - **In the owner's plan review,** four changes to Claude's first plan:
    - "usable" was defined too weakly;
    - RMSE was computed as sqrt(bias² + SD²) instead of directly;
    - Monte Carlo SEs were promised for every summary without a method for the nonlinear ones;
    - unbounded profile intervals were not evaluated.
  - **Claude's work caught:**
    - its own rough description of the design (45% zeros, against an exact 53.3%);
    - a crash in its investigation script;
    - a table too wide for the PDF;
    - provenance wording reused from the benchmark.
  - **A finding reported, not acted on:** the investigation showed that `zinbreg()` can converge to a local maximum (2 of 300 data sets). Claude did not change the engine, as agreed with the owner, and recorded it for the owner's decision.

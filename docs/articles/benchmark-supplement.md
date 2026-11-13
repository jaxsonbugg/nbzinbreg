# Benchmark: Supplementary Checks

The report’s Benchmarking section compares the required quantities of
the project description: the coefficients, their standard errors, theta
(NB) or log(theta) (ZINB), and the log-likelihood, each with a pass/fail
tolerance. This article adds comparisons that are **not graded
requirements**, so they have no tolerances. They are shown with two
measures of difference:

- the **absolute difference** \|ours − reference\|, the measure the
  requirements use;
- the **relative difference** \|ours − reference\| / \|reference\|,
  which puts quantities of different sizes on a common footing (a
  relative difference of 1e-5 means agreement to about five significant
  digits).

The data, formulas and settings are those of the report: our functions
at their default settings,
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) and
[`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html) with
tight tolerances.

## Standard Errors of theta and log(theta)

[`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md)
reports SE(theta) and
[`zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg.md)
reports SE(log theta), as the references do. Their agreement with the
references:

| Comparison | Quantity | Abs. difference | Rel. difference |
|:---|:---|---:|---:|
| bioChemists NB | SE(theta) | 1.13e-06 | 4.17e-06 |
| bioChemists ZINB (analytic Hessian) | SE(log theta) | 1.21e-06 | 8.91e-06 |
| NMES1988 NB | SE(theta) | 4.20e-07 | 1.25e-05 |
| NMES1988 ZINB (analytic Hessian) | SE(log theta) | 7.25e-08 | 2.07e-06 |

Table 1: Standard error of theta (NB) and of log(theta) (ZINB, analytic
Hessian), against the references.

## ZINB with the Numeric Hessian

[`zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg.md)
computes standard errors from the analytic Hessian (matrix of second
derivatives) of the log-likelihood by default; with
`zinbreg_control(hessian = "numeric")` it differentiates the
log-likelihood numerically instead
([`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html)).
The option only changes how the matrix used for the standard errors is
computed: the EM iterations and the gradient check that stops them are
the same, so the estimates are identical and this checks the standard
errors by a second, independent route.

| Comparison | Quantity | Max. abs. difference | Max. rel. difference |
|:---|:---|---:|---:|
| bioChemists ZINB (numeric Hessian) | coefficients | 4.75e-05 | 2.77e-04 |
| bioChemists ZINB (numeric Hessian) | standard errors | 2.82e-05 | 7.08e-05 |
| bioChemists ZINB (numeric Hessian) | log(theta) | 4.01e-06 | 4.11e-06 |
| bioChemists ZINB (numeric Hessian) | log-likelihood | 4.23e-09 | 2.73e-12 |
| bioChemists ZINB (numeric Hessian) | SE(log theta) | 1.21e-06 | 8.91e-06 |
| NMES1988 ZINB (numeric Hessian) | coefficients | 1.36e-05 | 1.48e-04 |
| NMES1988 ZINB (numeric Hessian) | standard errors | 6.38e-06 | 2.26e-05 |
| NMES1988 ZINB (numeric Hessian) | log(theta) | 2.33e-06 | 5.91e-06 |
| NMES1988 ZINB (numeric Hessian) | log-likelihood | 2.63e-09 | 2.17e-13 |
| NMES1988 ZINB (numeric Hessian) | SE(log theta) | 7.25e-08 | 2.07e-06 |
| bioChemists ZINB (analytic against numeric Hessian) | standard errors (ours: analytic against numeric Hessian) | 2.47e-07 | 2.72e-06 |
| NMES1988 ZINB (analytic against numeric Hessian) | standard errors (ours: analytic against numeric Hessian) | 5.93e-08 | 2.21e-07 |

Table 2: ZINB fits with the numeric Hessian against zeroinfl(), and our
analytic against our numeric standard errors.

Our two Hessian methods give standard errors that agree with each other
to 2.47e-07, far more closely than either agrees with zeroinfl() (at
least 6.38e-06). The remaining differences from zeroinfl() therefore
come from the small differences between the two packages’ estimates, not
from how the standard errors are computed.

## Validity Checks

Every comparison in the report and here is made only after these checks
pass for both fits of a pair. A failed check would be reported here with
its reason, and the comparison would be marked as not assessed.

**bioChemists NB**

| Check | Passed | Detail |
|:---|:--:|:---|
| our fit completed | yes | yes |
| reference fit completed | yes | yes |
| our fit converged | yes | converged = TRUE, 4 iterations |
| reference fit converged | yes | glm.nb converged = TRUE, theta warning: none |
| finite values (ours) | yes | all finite |
| finite values (reference) | yes | all finite |
| same observations | yes | 915 observations, identical responses |
| same parameterization | yes | 6 coefficients with identical names in the same order |
| valid covariance (ours) | yes | finite, symmetric, positive definite; condition number 4.46e+03 |
| valid covariance (reference) | yes | finite, symmetric, positive definite; condition number 4.46e+03 |

**bioChemists ZINB (analytic Hessian)**

| Check | Passed | Detail |
|:---|:--:|:---|
| our fit completed | yes | yes |
| reference fit completed | yes | yes |
| our fit converged | yes | converged = TRUE, gradient check passed = TRUE (Newton decrement 8.47e-09), diverging coefficients: none, 123 EM iterations |
| reference fit converged | yes | zeroinfl optim convergence code 0, converged = TRUE |
| finite values (ours) | yes | all finite |
| finite values (reference) | yes | all finite |
| same observations | yes | 915 observations, identical responses |
| same parameterization | yes | 12 coefficients with identical names in the same order |
| valid covariance (ours) | yes | finite, symmetric, positive definite; condition number 5.35e+05 |
| valid covariance (reference) | yes | finite, symmetric, positive definite; condition number 5.35e+05 |

**bioChemists ZINB (numeric Hessian)**

| Check | Passed | Detail |
|:---|:--:|:---|
| our fit completed | yes | yes |
| reference fit completed | yes | yes |
| our fit converged | yes | converged = TRUE, gradient check passed = TRUE (Newton decrement 8.47e-09), diverging coefficients: none, 123 EM iterations |
| reference fit converged | yes | zeroinfl optim convergence code 0, converged = TRUE |
| finite values (ours) | yes | all finite |
| finite values (reference) | yes | all finite |
| same observations | yes | 915 observations, identical responses |
| same parameterization | yes | 12 coefficients with identical names in the same order |
| valid covariance (ours) | yes | finite, symmetric, positive definite; condition number 5.35e+05 |
| valid covariance (reference) | yes | finite, symmetric, positive definite; condition number 5.35e+05 |

**NMES1988 NB**

| Check | Passed | Detail |
|:---|:--:|:---|
| our fit completed | yes | yes |
| reference fit completed | yes | yes |
| our fit converged | yes | converged = TRUE, 3 iterations |
| reference fit converged | yes | glm.nb converged = TRUE, theta warning: none |
| finite values (ours) | yes | all finite |
| finite values (reference) | yes | all finite |
| same observations | yes | 4406 observations, identical responses |
| same parameterization | yes | 8 coefficients with identical names in the same order |
| valid covariance (ours) | yes | finite, symmetric, positive definite; condition number 2.08e+03 |
| valid covariance (reference) | yes | finite, symmetric, positive definite; condition number 2.08e+03 |

**NMES1988 ZINB (analytic Hessian)**

| Check | Passed | Detail |
|:---|:--:|:---|
| our fit completed | yes | yes |
| reference fit completed | yes | yes |
| our fit converged | yes | converged = TRUE, gradient check passed = TRUE (Newton decrement 7.57e-09), diverging coefficients: none, 57 EM iterations |
| reference fit converged | yes | zeroinfl optim convergence code 0, converged = TRUE |
| finite values (ours) | yes | all finite |
| finite values (reference) | yes | all finite |
| same observations | yes | 4406 observations, identical responses |
| same parameterization | yes | 14 coefficients with identical names in the same order |
| valid covariance (ours) | yes | finite, symmetric, positive definite; condition number 1.07e+05 |
| valid covariance (reference) | yes | finite, symmetric, positive definite; condition number 1.07e+05 |

**NMES1988 ZINB (numeric Hessian)**

| Check | Passed | Detail |
|:---|:--:|:---|
| our fit completed | yes | yes |
| reference fit completed | yes | yes |
| our fit converged | yes | converged = TRUE, gradient check passed = TRUE (Newton decrement 7.57e-09), diverging coefficients: none, 57 EM iterations |
| reference fit converged | yes | zeroinfl optim convergence code 0, converged = TRUE |
| finite values (ours) | yes | all finite |
| finite values (reference) | yes | all finite |
| same observations | yes | 4406 observations, identical responses |
| same parameterization | yes | 14 coefficients with identical names in the same order |
| valid covariance (ours) | yes | finite, symmetric, positive definite; condition number 1.07e+05 |
| valid covariance (reference) | yes | finite, symmetric, positive definite; condition number 1.07e+05 |

## Timing with the Numeric Hessian

The numeric Hessian needs many more evaluations of the log-likelihood,
so it is slower. Timings are from the saved benchmark results (median of
10 runs, with the range).

| Data | Fit | Median (s) | Range (s) |
|:---|:---|---:|---:|
| bioChemists | zinbreg() (default settings, analytic Hessian) | 1.12 | 1.09 to 1.17 |
| bioChemists | zinbreg() (default settings, numeric Hessian) | 1.47 | 1.41 to 1.56 |
| NMES1988 | zinbreg() (default settings, analytic Hessian) | 2.35 | 2.29 to 2.4 |
| NMES1988 | zinbreg() (default settings, numeric Hessian) | 4.49 | 4.39 to 4.81 |

Table 3: zinbreg() with the analytic and the numeric Hessian:
computation time in seconds.

The installed package’s code fingerprint
(f6f29a27eff95d07ab6302d5a428954e) matches the code that produced the
timings (commit 65f2720b0fc178074ae477c7b23d55ba54768b71).


<!-- README.md is generated from README.Rmd. Please edit that file -->

# nbzinbreg

<!-- badges: start -->

[![R-CMD-check](https://github.com/jaxsonbugg/nbzinbreg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jaxsonbugg/nbzinbreg/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Negative binomial (NB) and zero-inflated negative binomial (ZINB)
regression in R, with the estimation written from scratch: iteratively
reweighted least squares (IRLS) and Newton-Raphson for NB, and the EM
algorithm with IRLS steps for ZINB. No general-purpose optimizer is
used. NB results are validated against `MASS::glm.nb()`; ZINB results
will be validated against `pscl::zeroinfl()`.

<!-- status: start -->

<!-- Temporary section: delete everything from "status: start" to
"status: end" once the package is finished and no work remains. -->

## Status

Negative binomial regression is implemented and validated; the
zero-inflated model is not implemented yet.

| Part | Status |
|----|----|
| Package skeleton (installs, documented, CI) | ✅ Done |
| `nbreg()`: negative binomial regression | ✅ Done |
| `zinbreg()`: zero-inflated NB regression | ⬜ Not started |
| Methods (`print`, `summary`, `predict`, …) | 🚧 In progress (done for `nbreg`) |
| Extras (confidence intervals, model tests, diagnostics) | ⬜ Not started |
| Benchmarks against `glm.nb()` and `zeroinfl()` | ⬜ Not started |
| Simulation study | ⬜ Not started |
| NMES1988 data analysis | ⬜ Not started |
| Report and website | ⬜ Not started |

Details: [progress
log](https://github.com/jaxsonbugg/nbzinbreg/blob/main/dev/progress_log.md)
· [design
notes](https://github.com/jaxsonbugg/nbzinbreg/blob/main/dev/design.md)
· [requirements
checklist](https://github.com/jaxsonbugg/nbzinbreg/blob/main/dev/requirements.md)

<!-- status: end -->

## Installation

``` r
# install.packages("devtools")
devtools::install_github("jaxsonbugg/nbzinbreg")
```

## Usage

### Negative binomial regression

``` r
library(nbzinbreg)
data("bioChemists", package = "pscl")

fit_nb <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
summary(fit_nb)
#> 
#> Call:
#> nbreg(formula = art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.2677 -0.8754 -0.2613  0.4983  6.6573 
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.256144   0.137348   1.865 0.062191 .  
#> femWomen    -0.216418   0.072636  -2.979 0.002887 ** 
#> marMarried   0.150489   0.082097   1.833 0.066791 .  
#> kid5        -0.176415   0.052813  -3.340 0.000837 ***
#> phd          0.015271   0.035873   0.426 0.670326    
#> ment         0.029082   0.003214   9.048  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 2.264 (SE 0.2712)
#> Log-likelihood: -1560.96 on 7 df
#> AIC: 3135.92;  BIC: 3169.65
#> Observations: 915;  residual df: 909
#> Outer passes: 4;  converged: yes
```

### Zero-inflated negative binomial regression

*Planned interface: this does not run yet.*

``` r
library(nbzinbreg)
data("bioChemists", package = "pscl")

# Same covariates in both parts; use y ~ x | z to choose the zero-part
# covariates separately.
fit_zinb <- zinbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
summary(fit_zinb)
```

## Check results

`R CMD check`: 0 errors, 0 warnings, 0 notes on all CI platforms (macOS,
Windows and Ubuntu; R release, devel and the previous release). A local
Windows run may show one NOTE, “unable to verify current time”, which
comes from the machine’s clock check, not from the package.

## License

MIT © Jaxson Bugg

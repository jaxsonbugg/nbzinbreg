
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nbzinbreg

<!-- badges: start -->

[![R-CMD-check](https://github.com/jaxsonbugg/nbzinbreg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jaxsonbugg/nbzinbreg/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Negative binomial (NB) and zero-inflated negative binomial (ZINB)
regression in R, with the estimation written from scratch: iteratively
reweighted least squares (IRLS) and Newton-Raphson for NB, and the EM
algorithm with IRLS steps for ZINB. No general-purpose optimizer is
used. NB results are validated against `MASS::glm.nb()`, and ZINB
results against `pscl::zeroinfl()`.

<!-- status: start -->

<!-- Temporary section: delete everything from "status: start" to
"status: end" once the package is finished and no work remains. -->

## Status

Negative binomial and zero-inflated negative binomial regression are
implemented and validated, with all required methods and the extras
(confidence intervals, model tests, bootstrap and diagnostics). The
benchmarks, simulation study, data analysis and report come next.

| Part                                                    | Status         |
|---------------------------------------------------------|----------------|
| Package skeleton (installs, documented, CI)             | ✅ Done        |
| `nbreg()`: negative binomial regression                 | ✅ Done        |
| `zinbreg()`: zero-inflated NB regression                | ✅ Done        |
| Methods (`print`, `summary`, `predict`, …)              | ✅ Done        |
| Extras (confidence intervals, model tests, diagnostics) | ✅ Done        |
| Benchmarks against `glm.nb()` and `zeroinfl()`          | ⬜ Not started |
| Simulation study                                        | ⬜ Not started |
| NMES1988 data analysis                                  | ⬜ Not started |
| Report and website                                      | ⬜ Not started |

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

``` r
library(nbzinbreg)
data("bioChemists", package = "pscl")

# Same covariates in both parts; use y ~ x | z to choose the zero-part
# covariates separately.
fit_zinb <- zinbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
summary(fit_zinb)
#> 
#> Call:
#> zinbreg(formula = art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.2942 -0.7601 -0.2909  0.4448  6.4155 
#> 
#> Count model coefficients (negative binomial with log link):
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.4167473  0.1435962   2.902  0.00371 ** 
#> femWomen    -0.1955061  0.0755925  -2.586  0.00970 ** 
#> marMarried   0.0975830  0.0844519   1.155  0.24789    
#> kid5        -0.1517319  0.0542061  -2.799  0.00512 ** 
#> phd         -0.0006999  0.0362696  -0.019  0.98460    
#> ment         0.0247861  0.0034924   7.097 1.27e-12 ***
#> Log(theta)   0.9763605  0.1354682   7.207 5.71e-13 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept)  -0.1917     1.3228  -0.145  0.88475   
#> femWomen      0.6359     0.8489   0.749  0.45378   
#> marMarried   -1.4994     0.9386  -1.597  0.11017   
#> kid5          0.6284     0.4428   1.419  0.15582   
#> phd          -0.0377     0.3080  -0.122  0.90257   
#> ment         -0.8823     0.3162  -2.790  0.00527 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 2.655 (SE 0.3596)
#> Log-likelihood: -1549.99 on 13 df
#> AIC: 3125.98;  BIC: 3188.63
#> Observations: 915;  residual df: 902
#> EM iterations: 123;  converged: yes
#> Gradient check: Newton decrement 8.47e-09 (grad_tol 1e-08)
#> Starting values: count part from negative binomial fit; zero part from logistic fit of I(y == 0)
```

### Inference and diagnostics

``` r
# Profile-likelihood interval for theta, Wald intervals for the coefficients
confint(fit_nb)
#>                2.5 %   97.5 %
#> (Intercept) -0.01305  0.52534
#> femWomen    -0.35878 -0.07405
#> marMarried  -0.01042  0.31140
#> kid5        -0.27993 -0.07290
#> phd         -0.05504  0.08558
#> ment         0.02278  0.03538
#> theta        1.80734  2.89853

# Is there overdispersion? (Poisson against NB)
overdispersion_test(fit_nb)
#> 
#>  Likelihood-ratio test for overdispersion (Poisson against negative
#>  binomial; boundary-corrected p-value)
#> 
#> data:  art ~ fem + mar + kid5 + phd + ment
#> LR = 180.2, p-value < 2.2e-16
#> alternative hypothesis: true theta is less than Inf
#> sample estimates:
#>    theta 
#> 2.264388

# Observed against expected frequencies of the first few counts
head(expected_counts(fit_zinb), 4)
#>   count observed expected
#> 1     0      275 285.4333
#> 2     1      246 233.9321
#> 3     2      178 165.3441
#> 4     3       84 101.2073
```

Model comparison with `anova()` and `vuong_test()` (including a
bootstrap likelihood-ratio test of NB against ZINB), bootstrap standard
errors (`boot_se()`), `simulate()`, `rootogram()`, and `tidy()` /
`glance()` are also available; see the function reference.

## Check results

`R CMD check`: 0 errors, 0 warnings, 0 notes on all CI platforms (macOS,
Windows and Ubuntu; R release, devel and the previous release). A local
Windows run may show one NOTE, “unable to verify current time”, which
comes from the machine’s clock check, not from the package.

## License

MIT © Jaxson Bugg

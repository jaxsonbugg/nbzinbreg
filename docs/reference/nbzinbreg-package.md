# nbzinbreg: Negative Binomial and Zero-Inflated Negative Binomial Regression

Fits negative binomial (NB2) regression by iteratively reweighted least
squares (IRLS) for the regression coefficients, alternated with
Newton-Raphson for the dispersion parameter, and zero-inflated negative
binomial regression by an EM algorithm whose M-steps are IRLS fits. No
general-purpose optimizer is used. Results are checked against the
reference implementations in 'MASS' and 'pscl'.

## See also

Useful links:

- <https://github.com/jaxsonbugg/nbzinbreg>

- <https://jaxsonbugg.github.io/nbzinbreg/>

- Report bugs at <https://github.com/jaxsonbugg/nbzinbreg/issues>

## Author

**Maintainer**: Jaxson Bugg <buggjm@miamioh.edu> \[copyright holder\]

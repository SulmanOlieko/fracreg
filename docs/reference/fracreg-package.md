# Fractional Response Regressions

Provides comprehensive tools for the estimation and specification
analysis of fractional response models. It supports univariate one-part,
hurdle two-part, and double-inflated three-part models. The package also
incorporates estimators for panel data settings (CRE, GMM, QML) and
addresses unobserved heterogeneity and endogeneity via correlated random
effects and control function approaches. It supports various link
functions, calculates average and conditional partial effects
analytically across all model types, and provides robust specification
tests such as RESET, P test, and GGOFF tests.

## Details

|          |            |
|----------|------------|
| Package: | fracreg    |
| Type:    | Package    |
| Version: | 1.0.0      |
| Date:    | 2026-07-05 |
| License: | GPL-3      |

## Acknowledgements

This package builds upon, consolidates, and modernises the fractional
regression frameworks originally implemented in the `frm`, `frmhet`, and
`frmpd` R packages developed by Joaquim J.S. Ramalho. As those original
packages have been deprecated and removed from the active CRAN
repository, `fracreg` serves as an actively maintained successor,
ensuring these econometric tools remain available to the R community.

Furthermore, we acknowledge James Ji (@f1kidd) and A. John Woodill
(@johnwoodill), the authors of the `fmlogit` R package on GitHub, whose
foundational work on fractional multinomial logit models inspired the
implementation of `fracregmlogit`. We also extend our gratitude to Ariel
Rokem and Kendrick Kay, the authors of the `fracridge` package, whose
methodological contributions to fractional ridge regression are
incorporated into the `fracregridge` functionalities of this package.

## References

Ramalho, J. J. S. *frm: Fractional Regression Models*. R package.
Formerly available on CRAN, currently archived.

Ramalho, J. J. S. *frmhet: Fractional Regression Models under
Heterogeneity*. R package. Formerly available on CRAN, currently
archived.

Ramalho, J. J. S. *frmpd: Fractional Regression Models for Panel Data*.
R package. Formerly available on CRAN, currently archived.

Papke, L. E. and Wooldridge, J. M. (1996), "Econometric methods for
fractional response variables with an application to 401(k) plan
participation rates", *Journal of Applied Econometrics*, 11(6), 619-632.

Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
estimating and testing empirical strategies for fractional response
models", *Journal of Economic Surveys*, 25(1), 19-68.

Ramalho, E. A., Ramalho, J. J. S., and Murteira, J. M. R. (2014), "A
two-part fractional response model for the spatial distribution of
vineyards in Portugal", *Journal of Applied Econometrics*, 29(4),
607-630.

Fang, K., & Ma, S. (2013), "Three-part model for fractional response
variables with application to Chinese household health insurance
coverage", *Journal of Applied Statistics*, 40(5), 925-940.

Ji, J., and Woodill, A. J., *fmlogit: Fractional Multinomial Logit*. R
package repository. \<https://github.com/f1kidd/fmlogit\>.

Rokem, A., and Kay, K., *fracridge: Fractional Ridge Regression*.
Package repository. \<https://github.com/nrdg/fracridge\>.

## See also

Useful links:

- <https://sulmanolieko.github.io/fracreg/>

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

# Fractional Response Models - Partial Effects

`fracreg.pe` is used to compute average and/or conditional partial
effects in fractional response models.

## Usage

``` r
fracreg.pe(object, APE = TRUE, CPE = FALSE, at = NULL, which.x = NULL, 
           variance = TRUE, table = TRUE)
```

## Arguments

- object:

  an object containing the results of an `fracreg` command.

- APE:

  a logical value indicating whether average partial effects are to be
  computed.

- CPE:

  a logical value indicating whether conditional partial effects are to
  be computed.

- at:

  a numeric vector containing the covariates' values at which the
  conditional partial effects are to be computed or the strings `"mean"`
  (the default) or `"median"`, in which cases the covariates are
  evaluated at their mean or median values (or mode, in case of dummy
  variables), respectively.

- which.x:

  a vector containing the names of the covariates to which the partial
  effects are to be computed.

- variance:

  a logical value indicating whether the variance of the estimated
  partial effects should be calculated. Defaults to `TRUE` whenever
  `table = TRUE`.

- table:

  a logical value indicating whether a summary table with the results
  should be printed.

## Details

`fracreg.pe` calculates partial effects for fractional response models
estimated via `fracreg`. `fracreg.pe` may be used to compute average or
conditional partial effects for: (i) one-part fractional response
models; (ii) the binary components of two-part and three-part fractional
response models; (iii) the fractional components of two-part and
three-part fractional response models; and (iv) two-part and three-part
fractional response models overall. For calculating standard errors, it
is taken into account the option that was previously chosen for
estimating the model. See Ramalho, Ramalho and Murteira (2011) and Fang
and Ma (2013) for details on the computation of partial effects in the
fractional response framework.

## Value

`fracreg.pe` returns a list with the following element:

- PE.p:

  a named vector of partial effects.

If `variance = TRUE` or `table = TRUE`, the previous list also contains
the following element:

- PE.sd:

  a named vector of standard errors of the estimated partial effects.

When both average and conditional partial effects are requested, two
lists containing the previous elements are returned, indexed by the
prefixes `ape` and `cpe`.

## References

Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
estimating and testing empirical strategies for fractional response
models", *Journal of Economic Surveys*, 25(1), 19-68.

Fang, K., & Ma, S. (2013), "Three-part model for fractional response
variables with application to Chinese household health insurance
coverage", *Journal of Applied Statistics*, 40(5), 925-940.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreg`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.md),
for fitting fractional response models.  
[`fracreg.reset`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.reset.md)
and
[`fracreg.ggoff`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.ggoff.md),
for specification tests.  
[`fracreg.ptest`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.ptest.md),
for non-nested hypothesis tests.

## Examples

``` r
N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1

#Computing average partial effects for a logit fractional response model
res <- fracreg(y,X,linkfrac="logit",table=FALSE)
fracreg.pe(res)
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                             Average partial effects 
#> -------------------------------------------------------------------------------- 
#>                              Fractional logit model 
#> -------------------------------------------------------------------------------- 
#>    Estimate Std. Error z value Pr(>|z|)    
#> X1 0.178302   0.009269   19.24   <2e-16 ***
#> X2 0.186834   0.011511   16.23   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 17:25:25 
#> -------------------------------------------------------------------------------- 

#Computing average partial effects for a binary logit + fractional probit
#two-part model
res <- fracreg(y,X,linkbin="logit",linkfrac="probit",type="2P",inf=1,table=FALSE)
fracreg.pe(res)
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                             Average partial effects 
#> -------------------------------------------------------------------------------- 
#>                 Binary logit + Fractional probit two-part model 
#> -------------------------------------------------------------------------------- 
#>    Estimate Std. Error z value Pr(>|z|)    
#> X1  0.10841    0.01274   8.511   <2e-16 ***
#> X2  0.12638    0.01486   8.502   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 17:25:25 
#> -------------------------------------------------------------------------------- 

#Computing conditional partial effects for X2 in the logit component
#of a two-part fractional response model, with the covariates evaluated
#at their median values
res <- fracreg(y,X,linkfrac="logit",type="2Pfrac",inf=1,table=FALSE)
fracreg.pe(res,APE=FALSE,CPE=TRUE,at="median",which.x="X2")
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                           Conditional partial effects 
#> -------------------------------------------------------------------------------- 
#>                  Fractional logit component of a two-part model 
#> -------------------------------------------------------------------------------- 
#>    Estimate Std. Error z value Pr(>|z|)    
#> X2  0.18783    0.01698   11.06   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 17:25:25 
#> -------------------------------------------------------------------------------- 
#> 
#> Note: covariates evaluated at median (or mode, for dummies) values

#Computing average partial effects for a three-part double-inflated model
y3p <- y
y3p[1:20] <- 0
y3p[21:40] <- 1
res3p <- fracreg(y3p,X,linkbin=c("logit","probit"),linkfrac="logit",type="3P",table=FALSE)
fracreg.pe(res3p)
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                             Average partial effects 
#> -------------------------------------------------------------------------------- 
#>        Three-part model - binary logit , binary probit + fractional logit 
#> -------------------------------------------------------------------------------- 
#>    Estimate Std. Error z value Pr(>|z|)    
#> X1  0.17427    0.01288   13.54   <2e-16 ***
#> X2  0.16785    0.01545   10.86   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 17:25:25 
#> -------------------------------------------------------------------------------- 
```

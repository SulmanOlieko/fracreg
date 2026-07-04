# P Test for Fractional Regression Models

`fracreg.ptest` is used to perform the P test to evaluate the
specification of alternative, non-nested fractional regression models by
testing against each other.

## Usage

``` r
fracreg.ptest(object1, object2, version = "Wald", table = TRUE)
```

## Arguments

- object1:

  an object containing the results of an `fracreg` command.

- object2:

  an object containing the results of another `fracreg` command.

- version:

  a vector containing the test versions to use. Available options:
  `Wald` (the default) and `LM`. Both options may be chosen at the same
  time and are computed in a robust way.

- table:

  a logical value indicating whether a summary table with the test
  results should be printed.

## Details

`fracreg.ptest` applies the P test statistic proposed by Davidson and
MacKinnon (1981) to fractional regression models estimated via
`fracreg`. `fracreg.ptest` may be used to test against each other two
alternative specifications for the link function in: (i) one-part
fractional regression models; (ii) the binary components of two-part and
three-part fractional regression models; (iii) the fractional components
of two-part and three-part fractional regression models; and (iv)
two-part and three-part fractional regression models. In addition,
`fracreg.ptest` may be used to test one-part models against two-part or
three-part models and in cases where the link functions are the same but
the regressors are non-nested. See Ramalho, Ramalho and Murteira (2011)
for details on the application of the P test in the fractional
regression framework.

## Value

`fracreg.reset` returns a named vector with the test results.

## References

Davidson, R. and J.G. MacKinnon (1981), "Several tests for model
specification on the presence of alternative hypotheses",
*Econometrica*, 49(3), 781-793.

Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
estimating and testing empirical strategies for fractional regression
models", *Journal of Economic Surveys*, 25(1), 19-68.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreg`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.md),
for fitting fractional regression models.  
[`fracreg.reset`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.reset.md)
and
[`fracreg.ggoff`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.ggoff.md),
for specification tests.  
[`fracreg.pe`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.pe.md),
for computing partial effects.

## Examples

``` r
N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1

#Testing logit versus loglog specifications for standard fractional
#regression models using a LM version of the P test
res1 <- fracreg(y,X,linkfrac="logit",table=FALSE)
res2 <- fracreg(y,X,linkfrac="loglog",table=FALSE)
fracreg.ptest(res1,res2,"LM")
#> 
#> *** P test ***
#> 
#> H0:  Fractional logit model
#> H1:  Fractional loglog model
#> 
#>  Version Statistic p-value  
#>       LM     0.930   0.335  
#> 
#> H0:  Fractional loglog model
#> H1:  Fractional logit model
#> 
#>  Version Statistic p-value  
#>       LM     3.582   0.058 *
#> 

#Testing a logit one-part fractional regression model versus a binary logit +
#fractional probit two-part model using a Wald version of the P test
res1 <- fracreg(y,X,linkfrac="logit",table=FALSE)
res2 <- fracreg(y,X,linkbin="logit",linkfrac="probit",type="2P",inf=1,table=FALSE)
fracreg.ptest(res1,res2,"Wald")
#> 
#> *** P test ***
#> 
#> H0:  Fractional logit model
#> H1:  Binary logit + Fractional probit two-part model
#> 
#>  Version Statistic p-value    
#>     Wald     2.194   0.029 ** 
#> 
#> H0:  Binary logit + Fractional probit two-part model
#> H1:  Fractional logit model
#> 
#>  Version Statistic p-value    
#>     Wald    12.341   0.000 ***
#> 
```

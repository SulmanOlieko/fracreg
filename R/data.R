#' 401(k) Plan Participation Data
#'
#' A cross-sectional dataset on 401(k) plan participation rates and firm characteristics, 
#' widely used in empirical applications of fractional response models 
#' (e.g., Papke & Wooldridge, 1996). The data is derived from the 'wooldridge' package.
#'
#' @format A data frame with 1,534 observations and 10 variables:
#' \describe{
#'   \item{prate}{Participation rate: proportion of eligible employees participating in the 401(k) plan (0 to 1).}
#'   \item{mrate}{Match rate: the firm's contribution matching rate per dollar.}
#'   \item{totpart}{Total number of participants.}
#'   \item{totelg}{Total number of eligible employees.}
#'   \item{age}{Age of the 401(k) plan in years.}
#'   \item{totemp}{Total number of firm employees.}
#'   \item{sole}{Indicator variable: 1 if the 401(k) is the sole retirement plan offered.}
#'   \item{ltotemp}{Natural log of total employees.}
#'   \item{age_sq}{Square of plan age.}
#'   \item{mrate_sq}{Square of match rate.}
#' }
#' @source Papke, L. E., & Wooldridge, J. M. (1996). "Econometric Methods for Fractional Response Variables with an Application to 401(k) Plan Participation Rates." Journal of Applied Econometrics, 11(6), 619-632.
#' @examples
#' data("fracreg_k401k")
#' summary(fracreg_k401k$prate)
"fracreg_k401k"

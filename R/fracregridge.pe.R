#' @title Partial Effects for Fractional Ridge Regression
#'
#' @description
#' Because Fractional Ridge Regression fits a linear model without a link function, the partial effects are mathematically identical to the estimated ridge coefficients. This function serves as a wrapper to maintain API compatibility with the rest of the \code{fracreg} package, printing a brief notification and returning the standard coefficient tables.
#' @param object An object of class \code{fracregridge}.
#' @param APE logical. Ignored for ridge regression.
#' @param CPE logical. Ignored for ridge regression.
#' @param at numeric vector. Ignored for ridge regression.
#' @param variance logical. Ignored for ridge regression.
#' @param table logical. Ignored for ridge regression.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' An object of class \code{fracreg.pe} containing the standard coefficient tables.
#'
#' @author Sulman Olieko Owili <oliekosulman@gmail.com>
#'
#' @seealso
#' \code{\link{fracregridge}}, \code{\link{fracreg.pe}}
#'
#' @examples
#' # Generate random data
#' set.seed(123)
#' y <- rnorm(100)
#' X <- matrix(rnorm(1000), 100, 10)
#' colnames(X) <- paste0("X", 1:10)
#' 
#' # Fit Fractional Ridge Regression
#' mod <- fracregridge(y, X, fracs = c(0.3, 0.5))
#' 
#' # Compute Partial Effects (identical to coefficients)
#' pe <- fracregridge.pe(mod)
#' print(pe)
#' @export
fracregridge.pe <- function(object, APE=TRUE, CPE=FALSE, at=NULL, variance=TRUE, table=FALSE, ...) {
    if(missing(object)) stop("object is missing")
    if(is.null(object$class) || object$class != "fracregridge") stop("object is not the output of a fracregridge command")
    
    res <- list(
        call = match.call(),
        type = "fracregridge.pe",
        class = "fracregridge.pe",
        table.info = object$table.info
    )
    
    class(res) <- c("fracregridge.pe", "fracregridge")
    return(res)
}

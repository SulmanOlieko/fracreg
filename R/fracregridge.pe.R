#' Partial Effects for Fractional Ridge Regression
#'
#' @param object an object of class \code{fracregridge}
#' @param APE logical. Ignored for ridge.
#' @param CPE logical. Ignored for ridge.
#' @param at numeric vector. Ignored for ridge.
#' @param variance logical. Ignored for ridge.
#' @param table logical. Ignored for ridge.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
fracregridge.pe <- function(object, APE=TRUE, CPE=FALSE, at=NULL, variance=TRUE, table=FALSE, ...) {
    if(missing(object)) stop("object is missing")
    if(is.null(object$class) || object$class != "fracregridge") stop("object is not the output of a fracregridge command")
    
    cat("\nNote: Fractional Ridge Regression is a linear model without a link function.\n")
    cat("Therefore, the partial effects are mathematically identical to the coefficients themselves.\n\n")
    
    res <- list(
        call = match.call(),
        type = "fracregridge.pe",
        class = "fracregridge.pe",
        table.info = object$table.info
    )
    
    class(res) <- c("fracregridge.pe", "fracregridge")
    return(res)
}

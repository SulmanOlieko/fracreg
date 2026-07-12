#' @title Extract Model Coefficients for fracreg
#' @description Extracts the estimated coefficients from a fitted \code{fracreg} model.
#' @param object A fitted model object of class \code{fracreg}.
#' @param ... Further arguments passed to or from other methods.
#' @return A named vector of coefficients.
#' @exportS3Method coef fracreg
coef.fracreg <- function(object, ...) {
    if(!is.null(object$p)) {
        return(object$p)
    } else if (!is.null(object$table.info$p)) {
        return(object$table.info$p)
    } else if (object$type %in% c("2P", "3P")) {
        # For multi-part models, coefficients are in the individual sub-model components
        warning("For multi-part models, use coef() on the individual sub-models (e.g., object$resBIN, object$resFRAC)")
        return(NULL)
    }
    return(NULL)
}

#' @title Extract Fitted Values for fracreg
#' @description Extracts the fitted conditional mean values from a \code{fracreg} model.
#' @param object A fitted model object of class \code{fracreg}.
#' @param ... Further arguments passed to or from other methods.
#' @return A numeric vector of fitted values.
#' @exportS3Method fitted fracreg
fitted.fracreg <- function(object, ...) {
    if (object$type == "2P") {
        return(object$yhat2P)
    } else if (object$type == "3P") {
        return(object$yhat3P)
    }
    return(object$yhat)
}

#' @title Extract Model Residuals for fracreg
#' @description Extracts the response residuals from a fitted \code{fracreg} model.
#' @param object A fitted model object of class \code{fracreg}.
#' @param ... Further arguments passed to or from other methods.
#' @return A numeric vector of residuals.
#' @exportS3Method residuals fracreg
residuals.fracreg <- function(object, ...) {
    if (object$type %in% c("2P", "3P")) {
        y <- object$ybase
    } else {
        y <- object$table.info$y
    }
    return(y - fitted(object))
}

#' @title Extract the Number of Observations for fracreg
#' @description Extracts the number of observations used to estimate a \code{fracreg} model.
#' @param object A fitted model object of class \code{fracreg}.
#' @param ... Further arguments passed to or from other methods.
#' @return An integer denoting the number of observations.
#' @exportS3Method nobs fracreg
nobs.fracreg <- function(object, ...) {
    if (object$type %in% c("2P", "3P")) {
        return(length(object$ybase))
    }
    return(length(object$table.info$y))
}

#' @title Extract Log-Likelihood for fracreg
#' @description Extracts the log-pseudolikelihood or log-likelihood from a fitted \code{fracreg} model.
#' @param object A fitted model object of class \code{fracreg}.
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class \code{logLik}.
#' @exportS3Method logLik fracreg
logLik.fracreg <- function(object, ...) {
    if (object$type %in% c("2P", "3P")) {
        warning("For multi-part models, use logLik() on the individual sub-models (e.g., object$resBIN, object$resFRAC)")
        return(NA)
    }
    if(is.null(object$table.info$LL)) {
        return(NA)
    }
    val <- object$table.info$LL
    attr(val, "df") <- length(coef(object))
    attr(val, "nobs") <- nobs(object)
    class(val) <- "logLik"
    return(val)
}

#' @title Extract Covariance Matrix for fracreg
#' @description Extracts the estimated variance-covariance matrix of the parameters.
#' @param object A fitted model object of class \code{fracreg}.
#' @param ... Further arguments passed to or from other methods.
#' @return A matrix of the estimated covariances.
#' @exportS3Method vcov fracreg
vcov.fracreg <- function(object, ...) {
    if(!is.null(object$p.var)) {
        return(object$p.var)
    } else if (object$type %in% c("2P", "3P")) {
        warning("For multi-part models, use vcov() on the individual sub-models (e.g., object$resBIN, object$resFRAC)")
        return(NULL)
    }
    return(NULL)
}

#' @title Predict Method for fracreg
#' @description Predicts conditional mean values from a fitted \code{fracreg} model.
#' @param object A fitted model object of class \code{fracreg}.
#' @param newdata An optional data frame or matrix in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param ... Further arguments passed to or from other methods.
#' @return A numeric vector of predicted values.
#' @exportS3Method predict fracreg
predict.fracreg <- function(object, newdata = NULL, ...) {
    if (is.null(newdata)) {
        return(fitted(object))
    }
    
    if (object$type %in% c("2P", "3P")) {
        stop("predict() with newdata is currently not supported natively for multi-part aggregate models. You can predict on the individual sub-components.")
    }
    
    # We must construct the linear predictor
    x.names <- object$x.names
    
    # Simple check to map data frame to matrix
    if(is.data.frame(newdata)) {
        # Check if intercept is needed
        if("(Intercept)" %in% x.names && !"(Intercept)" %in% colnames(newdata)) {
            newdata <- cbind(1, as.matrix(newdata))
            colnames(newdata)[1] <- "(Intercept)"
        } else {
            newdata <- as.matrix(newdata)
        }
    } else if (is.matrix(newdata)) {
        if("(Intercept)" %in% x.names && !"(Intercept)" %in% colnames(newdata)) {
            newdata <- cbind(1, newdata)
            colnames(newdata)[1] <- "(Intercept)"
        }
    }
    
    # Ensure correct ordering and matching
    if(!all(x.names %in% colnames(newdata))) {
        stop("newdata does not contain all the necessary columns")
    }
    
    X <- newdata[, x.names, drop=FALSE]
    beta <- coef(object)
    
    xb <- as.vector(X %*% beta)
    
    link <- object$link
    pred <- fracreg.links(link)$linkinv(xb)
    return(pred)
}

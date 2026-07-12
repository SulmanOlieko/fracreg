#' @title Extract Model Coefficients for fracregpd
#' @description Extracts the estimated coefficients from a fitted \code{fracregpd} model.
#' @param object A fitted model object of class \code{fracregpd}.
#' @param ... Further arguments passed to or from other methods.
#' @return A named vector of coefficients.
#' @exportS3Method coef fracregpd
coef.fracregpd <- function(object, ...) {
    if(!is.null(object$p)) {
        return(object$p)
    } else if (!is.null(object$table.info$p)) {
        return(object$table.info$p)
    }
    return(NULL)
}

#' @title Extract Fitted Values for fracregpd
#' @description Extracts the fitted conditional mean values from a \code{fracregpd} model.
#' @param object A fitted model object of class \code{fracregpd}.
#' @param ... Further arguments passed to or from other methods.
#' @return A numeric vector of fitted values.
#' @exportS3Method fitted fracregpd
fitted.fracregpd <- function(object, ...) {
    link <- object$link
    yhat <- fracregpd.links(link)$linkinv(object$xbhat)
    return(yhat)
}

#' @title Extract Model Residuals for fracregpd
#' @description Extracts the response residuals from a fitted \code{fracregpd} model.
#' @param object A fitted model object of class \code{fracregpd}.
#' @param ... Further arguments passed to or from other methods.
#' @return A numeric vector of residuals.
#' @exportS3Method residuals fracregpd
residuals.fracregpd <- function(object, ...) {
    type <- object$type
    link <- object$link
    Hy <- object$Hy
    
    if (type == "QMLcre") {
        y <- Hy
    } else {
        if (link == "logit") {
            y <- Hy / (1 + Hy)
        } else if (link == "cloglog") {
            y <- 1 - exp(-Hy)
        } else {
            y <- Hy # Fallback
        }
    }
    
    return(y - fitted(object))
}

#' @title Extract the Number of Observations for fracregpd
#' @description Extracts the number of observations used to estimate a \code{fracregpd} model.
#' @param object A fitted model object of class \code{fracregpd}.
#' @param ... Further arguments passed to or from other methods.
#' @return An integer denoting the number of observations.
#' @exportS3Method nobs fracregpd
nobs.fracregpd <- function(object, ...) {
    return(length(object$Hy))
}

#' @title Extract Log-Likelihood for fracregpd
#' @description Extracts the log-pseudolikelihood or log-likelihood from a fitted \code{fracregpd} model.
#' @param object A fitted model object of class \code{fracregpd}.
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class \code{logLik}.
#' @exportS3Method logLik fracregpd
logLik.fracregpd <- function(object, ...) {
    if(is.null(object$table.info$LL)) {
        return(NA)
    }
    val <- object$table.info$LL
    attr(val, "df") <- length(coef(object))
    attr(val, "nobs") <- nobs(object)
    class(val) <- "logLik"
    return(val)
}

#' @title Extract Covariance Matrix for fracregpd
#' @description Extracts the estimated variance-covariance matrix of the parameters.
#' @param object A fitted model object of class \code{fracregpd}.
#' @param ... Further arguments passed to or from other methods.
#' @return A matrix of the estimated covariances.
#' @exportS3Method vcov fracregpd
vcov.fracregpd <- function(object, ...) {
    if(!is.null(object$p.var)) {
        return(object$p.var)
    }
    return(NULL)
}

#' @title Predict Method for fracregpd
#' @description Predicts conditional mean values from a fitted \code{fracregpd} model.
#' @param object A fitted model object of class \code{fracregpd}.
#' @param newdata An optional data frame or matrix in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param ... Further arguments passed to or from other methods.
#' @return A numeric vector of predicted values.
#' @exportS3Method predict fracregpd
predict.fracregpd <- function(object, newdata = NULL, ...) {
    if (is.null(newdata)) {
        return(fitted(object))
    }
    
    # We must construct the linear predictor
    x.names <- object$table.info$x.names
    
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
    pred <- fracregpd.links(link)$linkinv(xb)
    return(pred)
}

#' @title Extract Model Coefficients for fracregridge
#' @description Extracts the estimated coefficients from a fitted \code{fracregridge} model.
#' @param object A fitted model object of class \code{fracregridge}.
#' @param ... Further arguments passed to or from other methods.
#' @return A matrix or array of coefficients.
#' @exportS3Method coef fracregridge
coef.fracregridge <- function(object, ...) {
    return(object$coef)
}

#' @title Extract Fitted Values for fracregridge
#' @description Extracts the fitted conditional mean values from a \code{fracregridge} model.
#' @param object A fitted model object of class \code{fracregridge}.
#' @param ... Further arguments passed to or from other methods.
#' @return A matrix or array of fitted values.
#' @exportS3Method fitted fracregridge
fitted.fracregridge <- function(object, ...) {
    return(object$yhat)
}

#' @title Extract Model Residuals for fracregridge
#' @description Extracts the response residuals from a fitted \code{fracregridge} model.
#' @param object A fitted model object of class \code{fracregridge}.
#' @param ... Further arguments passed to or from other methods.
#' @return A matrix or array of residuals.
#' @exportS3Method residuals fracregridge
residuals.fracregridge <- function(object, ...) {
    y <- object$y
    yhat <- fitted(object)
    
    if (length(dim(yhat)) == 2) {
        # yhat is N x F
        res <- yhat
        y_vec <- as.numeric(y)
        for(f in 1:ncol(yhat)) {
            res[, f] <- y_vec - yhat[, f]
        }
        return(res)
    }
    
    if (!is.matrix(y)) {
        y <- as.matrix(y)
    }
    
    # y is now N x J
    # yhat is N x F x J
    res <- yhat
    for(i in 1:dim(yhat)[3]) {
        for(f in 1:dim(yhat)[2]) {
            res[, f, i] <- y[, i] - yhat[, f, i]
        }
    }
    return(res)
}

#' @title Extract the Number of Observations for fracregridge
#' @description Extracts the number of observations used to estimate a \code{fracregridge} model.
#' @param object A fitted model object of class \code{fracregridge}.
#' @param ... Further arguments passed to or from other methods.
#' @return An integer denoting the number of observations.
#' @exportS3Method nobs fracregridge
nobs.fracregridge <- function(object, ...) {
    return(object$stats.info[[1]]$n_obs)
}

#' @title Extract Log-Likelihood for fracregridge
#' @description Extracts the log-likelihood from a fitted \code{fracregridge} model.
#' @param object A fitted model object of class \code{fracregridge}.
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class \code{logLik}.
#' @exportS3Method logLik fracregridge
logLik.fracregridge <- function(object, ...) {
    warning("logLik() is not applicable for fracregridge models because they are estimated via an analytical SVD penalty solver (Ridge), not MLE.")
    return(NA)
}

#' @title Extract Covariance Matrix for fracregridge
#' @description Extracts the estimated variance-covariance matrices of the parameters.
#' @param object A fitted model object of class \code{fracregridge}.
#' @param ... Further arguments passed to or from other methods.
#' @return A matrix of the estimated covariances.
#' @exportS3Method vcov fracregridge
vcov.fracregridge <- function(object, ...) {
    warning("vcov() is not natively stored for fracregridge models across all fractions. See summary(object) for individual fraction standard errors.")
    return(NULL)
}

#' @title Predict Method for fracregridge
#' @description Predicts conditional mean values from a fitted \code{fracregridge} model.
#' @param object A fitted model object of class \code{fracregridge}.
#' @param newdata An optional data frame or matrix in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param ... Further arguments passed to or from other methods.
#' @return A matrix or array of predicted values.
#' @exportS3Method predict fracregridge
predict.fracregridge <- function(object, newdata = NULL, ...) {
    if (is.null(newdata)) {
        return(fitted(object))
    }
    
    # Simple check to map data frame to matrix
    x.names <- colnames(object$x)
    if(is.data.frame(newdata)) {
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
    
    if (is.matrix(beta)) {
        # Single equation
        return(X %*% beta)
    } else if (is.array(beta)) {
        # Multivariate
        pred <- array(NA, dim = c(nrow(X), dim(beta)[2], dim(beta)[3]))
        for(i in 1:dim(beta)[3]) {
            pred[,,i] <- X %*% beta[,,i]
        }
        return(pred)
    }
}

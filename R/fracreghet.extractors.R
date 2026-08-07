#' @title Extract Model Coefficients for fracreghet
#' @description Extracts the estimated coefficients from a fitted \code{fracreghet} model.
#' @param object A fitted model object of class \code{fracreghet}.
#' @param ... Further arguments passed to or from other methods.
#' @return A named vector of coefficients.
#' @exportS3Method coef fracreghet
coef.fracreghet <- function(object, ...) {
    if(!is.null(object$p)) {
        return(object$p)
    } else if (!is.null(object$table.info$p)) {
        return(object$table.info$p)
    }
    return(NULL)
}

#' @title Extract Fitted Values for fracreghet
#' @description Extracts the fitted conditional mean values from a \code{fracreghet} model.
#' @param object A fitted model object of class \code{fracreghet}.
#' @param ... Further arguments passed to or from other methods.
#' @return A numeric vector of fitted values.
#' @exportS3Method fitted fracreghet
fitted.fracreghet <- function(object, ...) {
    link <- object$link
    yhat <- fracreghet.links(link)$linkinv(object$xbhat)
    return(yhat)
}

#' @title Extract Model Residuals for fracreghet
#' @description Extracts the response residuals from a fitted \code{fracreghet} model.
#' @param object A fitted model object of class \code{fracreghet}.
#' @param ... Further arguments passed to or from other methods.
#' @return A numeric vector of residuals.
#' @exportS3Method residuals fracreghet
residuals.fracreghet <- function(object, ...) {
    type <- object$type
    link <- object$link
    Hy <- object$Hy
    adjust <- object$adjust
    
    if (type == "QMLxv") {
        y <- Hy
    } else if (type %in% c("LINx", "LINz", "LINxv")) {
        # Hy = linkfun(y_adj) => y_adj = linkinv(Hy) => y = linkinv(Hy) - adjust
        y_adj <- fracreghet.links(link)$linkinv(Hy)
        if (is.numeric(adjust)) {
            y <- y_adj - adjust
        } else {
            y <- y_adj
        }
    } else {
        # GMM estimators
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

#' @title Extract the Number of Observations for fracreghet
#' @description Extracts the number of observations used to estimate a \code{fracreghet} model.
#' @param object A fitted model object of class \code{fracreghet}.
#' @param ... Further arguments passed to or from other methods.
#' @return An integer denoting the number of observations.
#' @exportS3Method nobs fracreghet
nobs.fracreghet <- function(object, ...) {
    return(length(object$Hy))
}

#' @title Extract Log-Likelihood for fracreghet
#' @description Extracts the log-pseudolikelihood or log-likelihood from a fitted \code{fracreghet} model.
#' @param object A fitted model object of class \code{fracreghet}.
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class \code{logLik}.
#' @exportS3Method logLik fracreghet
logLik.fracreghet <- function(object, ...) {
    if (object$type != "QMLxv") {
        return(NA) # GMM/LIN do not have log-likelihood
    }
    yhat <- fitted(object)
    y <- object$Hy # Because for QMLxv, Hy is y
    eps <- 1e-16
    LL <- sum(ifelse(y > 0, y * log(pmax(yhat, eps)), 0) + ifelse(y < 1, (1-y) * log(pmax(1-yhat, eps)), 0))
    
    attr(LL, "df") <- length(coef(object))
    attr(LL, "nobs") <- nobs(object)
    class(LL) <- "logLik"
    return(LL)
}

#' @title Extract Covariance Matrix for fracreghet
#' @description Extracts the estimated variance-covariance matrix of the parameters.
#' @param object A fitted model object of class \code{fracreghet}.
#' @param ... Further arguments passed to or from other methods.
#' @return A matrix of the estimated covariances.
#' @exportS3Method vcov fracreghet
vcov.fracreghet <- function(object, ...) {
    if(!is.null(object$p.var)) {
        return(object$p.var)
    } else if (!is.null(object$table.info$p.var)) {
        return(object$table.info$p.var)
    }
    return(NULL)
}

#' @title Predict Method for fracreghet
#' @description Predicts conditional mean values from a fitted \code{fracreghet} model.
#' @param object A fitted model object of class \code{fracreghet}.
#' @param newdata An optional data frame or matrix in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param ... Further arguments passed to or from other methods.
#' @return A numeric vector of predicted values.
#' @exportS3Method predict fracreghet
predict.fracreghet <- function(object, newdata = NULL, ...) {
    if (is.null(newdata)) {
        return(fitted(object))
    }
    
    stop("predict() with newdata is currently not natively supported for fracreghet due to the secondary variance equation. Fitted values are available without newdata.")
}

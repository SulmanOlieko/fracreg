#' @title Extract Model Coefficients for fracregmlogit
#' @description Extracts the estimated coefficients from a fitted \code{fracregmlogit} model.
#' @param object A fitted model object of class \code{fracregmlogit}.
#' @param ... Further arguments passed to or from other methods.
#' @return A matrix of coefficients for each choice equation.
#' @exportS3Method coef fracregmlogit
coef.fracregmlogit <- function(object, ...) {
    return(object$coefficient)
}

#' Extract Fitted Values, Residuals, and Predictions
#' 
#' @name fitted.fracregmlogit
#' @aliases residuals.fracregmlogit
#' @aliases predict.fracregmlogit
#' Extract fitted dependent variable, residuals, or predictions from a fractional multinomial logit model. 
#' @param object A "fracregmlogit" object.
#' @param newdata A new X matrix to perform model prediction. If NULL, defaults to the original dataset. 
#' X can be a vector with length k, or a matrix with k columns, where k is the number of explanatory 
#' variables in the original model. 
#' @param newbeta A new augmented matrix of coefficients that can be used to predict outcome variables. 
#' Feeds into \code{object$coefficient}, which contains the baseline coefficient. Useful for constructing
#' confidence intervals via simulation or bootstrapping. 
#' @param ... Additional arguments.
#' @return An object of class \code{data.frame} containing numeric values where each column corresponds to one of the choice alternatives in the response variable matrix and each row corresponds to an observation. Specifically:
#' \itemize{
#'   \item \code{fitted}: Returns the estimated fitted fractional response values (choice shares or predicted conditional probabilities).
#'   \item \code{residuals}: Returns the response residuals (the actual observed shares minus the estimated fitted shares).
#'   \item \code{predict}: Returns the predicted choice shares or conditional probabilities computed from the specified model object and \code{newdata} or \code{newbeta}.
#' }
#' @seealso \code{\link{fracregmlogit}}
#' @examples
#' data("fracreg_spending")
#' df <- na.omit(fracreg_spending)
#' X = df[,2:5]
#' y = df[,6:11]
#' results1 = fracregmlogit(y, X)
#' 
#' # Extract fitted values
#' fit = fitted(results1)
#' 
#' # Extract residuals
#' res = residuals(results1)
#' 
#' # Predict using the first observation from the original dataset
#' pred = predict(results1, newdata = X[1,])
#' @rdname fitted.fracregmlogit
#' @exportS3Method fitted fracregmlogit
fitted.fracregmlogit <- function(object, ...) {
  j=length(object$estimates)+1; k=dim(object$estimates[[1]])[1]; N=dim(object$y)[1]
  betamat_aug = object$coefficient; X=object$X; y=object$y
  sum_expxb = rowSums(exp(X %*% t(betamat_aug))) # sum of the exp(x'b)s
  yhat = y
  for(i in 1:j){
    expxb = exp(X %*% betamat_aug[i,]) # individual exp(x'b)
    yhat[,i] = expxb / sum_expxb
  }
  return(as.data.frame(yhat))
}

#' @rdname fitted.fracregmlogit
#' @exportS3Method residuals fracregmlogit
residuals.fracregmlogit <- function(object, ...) {
  yhat = fitted(object)
  return(as.data.frame(object$y-yhat))
}

#' @rdname fitted.fracregmlogit
#' @exportS3Method predict fracregmlogit
predict.fracregmlogit <- function(object, newdata=NULL, newbeta = NULL, ...) {
  if(length(newdata)==0) return(fitted(object))
  if(length(newbeta)>0) object$coefficient = newbeta
  j=length(object$estimates)+1; k=dim(object$estimates[[1]])[1]; N=dim(object$y)[1]
  betamat_aug = object$coefficient;
  newdata = as.matrix(newdata)
  if(length(newdata) == dim(newdata)[1]) newdata = t(newdata) # vector
  if(k != dim(newdata)[2]+1) stop(paste("Dimension of newdata is wrong. Should be",k-1,"instead of",dim(newdata)[2]))
  X = cbind(newdata,1); N = dim(X)[1]
  yhat = matrix(ncol=j,nrow=N); colnames(yhat) = colnames(object$y)
  sum_expxb = rowSums(exp(X %*% t(betamat_aug))) # sum of the exp(x'b)s
  for(i in 1:j){
    expxb = exp(X %*% betamat_aug[i,]) # individual exp(x'b)
    yhat[,i] = expxb / sum_expxb
  }
  return(as.data.frame(yhat))
}

#' @title Extract the Number of Observations for fracregmlogit
#' @description Extracts the number of observations used to estimate a \code{fracregmlogit} model.
#' @param object A fitted model object of class \code{fracregmlogit}.
#' @param ... Further arguments passed to or from other methods.
#' @return An integer denoting the number of observations.
#' @exportS3Method nobs fracregmlogit
nobs.fracregmlogit <- function(object, ...) {
    return(unname(object$count["Obs"]))
}

#' @title Extract Log-Likelihood for fracregmlogit
#' @description Extracts the log-pseudolikelihood from a fitted \code{fracregmlogit} model.
#' @param object A fitted model object of class \code{fracregmlogit}.
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class \code{logLik}.
#' @exportS3Method logLik fracregmlogit
logLik.fracregmlogit <- function(object, ...) {
    if(is.null(object$likelihood)) {
        return(NA)
    }
    val <- object$likelihood
    # df is roughly number of choices * number of params
    k <- dim(object$coefficient)[2]
    j <- dim(object$coefficient)[1] - 1
    attr(val, "df") <- k * j
    attr(val, "nobs") <- nobs(object)
    class(val) <- "logLik"
    return(val)
}

#' @title Extract Covariance Matrix for fracregmlogit
#' @description Extracts the estimated variance-covariance matrices of the parameters.
#' @param object A fitted model object of class \code{fracregmlogit}.
#' @param ... Further arguments passed to or from other methods.
#' @return A list of covariance matrices for each choice equation.
#' @exportS3Method vcov fracregmlogit
vcov.fracregmlogit <- function(object, ...) {
    return(object$vcov)
}

#' Fractional Multinomial Logit Average Partial Effects
#' 
#' Calculate average partial effects (APE) of independent variables from a fractional multinomial logit model. 
#' 
#' @param object A "fracregmlogit" object.
#' @param effect Can be "marginal", for marginal effects; or "discrete", for discrete changes from
#' the min to the max. 
#' @param marg.type Type of marginal or discrete effects to be computed. Default to "atmean", the effect at 
#' the mean of all covariates. Also takes "aveacr", the averaged effects across all observations. See details. 
#' @param se Whether to calculate standard errors for those margins. See details. 
#' @param varlist A string vector which provides the names of variables to calculate 
#' the marginal effect for. If missing, all variables except the constant will be calculated. 
#' Use "constant" if you wish to compute the marginal effect of the constant. 
#' @param at Specify values of the X-matrix at which the partial effect will be retrieved. Expects a vector input
#' of length K-1. Only supported for \code{marg.type="atmean"}. See \code{predict.fracregmlogit(newdata)}. 
#' @param R Number of times to sample for the Krinsky-Robb standard error. Default to 1000. 
#' @details This module calculates the average partial effects (APEs) from a fractional multinomial logit model.
#' Partial effects are the counterpart of the marginal effects in a linear model setting. In linear models, 
#' usually the parameter estimate itself represents the marginal effect (if the variable in question is continuous). 
#' In logit models, however, the parameter estimate at hand is the effect on the log-ratio between the choice variable
#' and the baseline variable. This function is intended to extract APEs from the 
#' coefficient estimates computed from the fractional multinomial logit models.
#' 
#' This function allows for two types of partial effects: marginal effects, and discrete effects.
#' A marginal effect represents how a unit change in one continuous variable x may influence the choice variable y. 
#' The estimation of marginal effects is very straightforward. However, special care is needed when averaging 
#' the marginal effect across observations to acquire the APE. One approach is to use the estimate of the marginal effect while setting
#' other explanatory variables at the mean. We call this the marginal effect at the mean (MEM), which corresponds
#' to the option \code{marg.type="atmean"}. Another approach is to take the average of marginal effects for each
#' individual. We call this the average marginal effect (AME), which corresponds to the option \code{marg.type="aveacr"}. 
#' 
#' The discrete effect represents how a discrete change in one specific x, discrete or continuous, influences the choice variable y. 
#' This is more useful for categorical variables, as calculating the "marginal effect" makes little sense
#' for them. In this function, we calculate the discrete effect by changing the explanatory variable from 
#' its minimum to its maximum. For a binary variable, this is just the difference between 0 and 1. Similar 
#' to the marginal effect case, we also have the discrete effect at the mean (DEM), corresponding to \code{marg.type="atmean"}
#' and the average discrete effect (ADE), corresponding to \code{marg.type="aveacr"}.
#' 
#' Standard errors are provided for the effects by using the Krinsky-Robb (KR) method. Krinsky-Robb is a simulation-based
#' method that calculates the empirical value of a function given a known distribution of its variables. Here 
#' we provide Krinsky-Robb standard errors for MEM and DEM, and the user can specify how many times of 
#' simulation \code{R} the Krinsky-Robb algorithm should run. 
#' 
#' The user can also specify a subset of explanatory variables when calculating effects. This is done through
#' specifying string vectors containing the column names of the explanatory variables to \code{varlist}. As the
#' KR standard error can be computationally intensive, it is advised to calculate it only for the variables of interest. 
#' 
#' @return The function returns an object of class "fracregmlogit.pe". It contains the following components:
#' @return \code{effects} A matrix of calculated effects.
#' @return \code{se} A matrix of standard errors corresponding to the effects. Shows up if se=TRUE for the 
#' input parameter.
#' @return \code{ztable} A list of matrices containing effects, standard errors, z-stats and p-values.
#' @return \code{R} Number of simulation times for Krinsky-Robb standard error calculation. Null if se=FALSE.  
#' @return \code{expl} String message explaining the effects calculated.    
#' 
#' @seealso \code{\link{fracregmlogit}} for the model estimation, \code{\link{plot.fracregmlogit.pe}} for plotting effects.
#' 
#' @examples
#' data("fracreg_spending")
#' X = fracreg_spending[,2:5]
#' y = fracreg_spending[,6:11]
#' results1 = fracregmlogit(y, X)
#' 
#' # Calculate marginal effects at the mean (without standard errors for speed)
#' pe_marg = fracregmlogit.pe(results1, effect="marginal", se=FALSE)
#' 
#' # Calculate discrete effects for specific variables with standard errors
#' pe_disc = fracregmlogit.pe(results1, effect="discrete", 
#'                            varlist = colnames(results1$X)[c(1,3)], 
#'                            se=TRUE, R=50)
#' summary(pe_disc)
#' @export fracregmlogit.pe
fracregmlogit.pe<-function(object,effect=c("marginal","discrete"),
                          marg.type="atmean",se=FALSE,varlist = NULL,at=NULL,R=1000){
  effect = match.arg(effect)
  j=length(object$estimates)+1; K=dim(object$estimates[[1]])[1]; N=dim(object$y)[1]
  betamat = object$coefficient
  R = R # for Krinsky-Robb sampling
  # determine variables
  Xnames = colnames(object$X); ynames = colnames(object$y)
  if(length(varlist)==0){
    varlist=Xnames[-K]
    var_colNo = c(1:(K-1))
    k = length(var_colNo)
  }else{
    var_colNo = unlist(lapply(varlist, function(x) {which(Xnames == x)}))
    if(length(varlist) != length(var_colNo)) stop("Unrecognized varlist input. Please double check your spelling")
    k = length(var_colNo)
  }
  
  xmarg = matrix(ncol=k,nrow=j)
  se_mat = matrix(ncol=k,nrow=j)
  marg_list = list()
  
  # Pre-draw Multivariate Normal vectors for Krinsky-Robb sampling if required
  if (se == TRUE) {
    # Initialize a 3D array: R simulations x J choices x (K+1) parameters
    betamat_R = array(0, dim=c(R, j, K))
    for (i in 2:j) {
      # Use MASS::mvrnorm to capture covariances between variables within each equation
      betamat_R[, i, ] = MASS::mvrnorm(R, betamat[i, ], object$vcov[[i]])
    }
    # Note: betamat_R[, 1, ] remains perfectly 0 as required for the baseline choice.
  }

  if(effect == "marginal"){
    # calculate marginal effects
    yhat = predict(object); yhat = as.matrix(yhat)
    for(c in var_colNo){
      c1 = which(var_colNo == c)
      if(marg.type == "aveacr"){
        # this is the average marginal effect for all observations
        beta_bar = as.vector(yhat %*% betamat[,c])
        betak_long = matrix(rep(betamat[,c],N),nrow=N,byrow=TRUE)
        marg_mat =  yhat * (betak_long-beta_bar)
        xmarg[,c1] = colMeans(marg_mat)
        marg_list[[c1]] = marg_mat
      }
      if(marg.type == "atmean"){
        # this is the marginal effect at the mean
        if(is.null(at)) at = colMeans(object$X[,-K])
        yhat_mean = predict(object,newdata=at)
        beta_bar = sum(yhat_mean * betamat[,c])
        betak = betamat[,c]
        marg_vec = yhat_mean * (betak - beta_bar)
        xmarg[,c1] = as.numeric(marg_vec) 
      }
    }
    
    if (se == TRUE) {
      # Vectorized Multivariate Krinsky-Robb standard error calculation for marginal effects
      x_mean = c(colMeans(object$X[,-K]), 1)
      
      Xbeta_R = matrix(0, nrow=R, ncol=j)
      for (i in 2:j) {
        Xbeta_R[, i] = betamat_R[, i, ] %*% x_mean
      }
      exp_Xbeta_R = exp(Xbeta_R)
      yhat_R = exp_Xbeta_R / rowSums(exp_Xbeta_R) # R x j
      
      for (c in var_colNo) {
        c1 = which(var_colNo == c)
        
        # Calculate the beta_bar (weighted average parameter) for each of the R draws
        beta_bar_R = rowSums(yhat_R * betamat_R[, , c])
        
        # Calculate the marginal effect spread across the R draws for each choice j
        for (i in 1:j) {
          marg_matrix_i = yhat_R[, i] * (betamat_R[, i, c] - beta_bar_R)
          se_mat[i, c1] = sd(marg_matrix_i)
        }
      }
    }
  }
  
  if(effect=="discrete"){
    for(c in var_colNo){
      c1 = which(var_colNo == c)
      if(marg.type == "aveacr"){
        Xmin <- Xmax <- object$X[,-K]
        Xmin[,c] = min(object$X[,c])
        Xmax[,c] = max(object$X[,c])
        yhat_min = predict(object,newdata=Xmin)
        yhat_max = predict(object,newdata=Xmax)
        ydisc = yhat_max - yhat_min
        xmarg[,c1] = colMeans(ydisc)
        marg_list[[c1]] = ydisc
      }
      if(marg.type == "atmean"){
        if(is.null(at)) at = colMeans(object$X[,-K])
        Xmin <- Xmax <- at
        Xmin[c] = min(object$X[,c])
        Xmax[c] = max(object$X[,c])
        yhat_min = predict(object,newdata=Xmin)
        yhat_max = predict(object,newdata=Xmax)
        ydisc = yhat_max - yhat_min
        xmarg[,c1] = as.numeric(ydisc)
      }
    }
    
    if (se == TRUE) {
      # Vectorized Multivariate Krinsky-Robb standard error calculation for discrete effects
      for (c in var_colNo) {
        c1 = which(var_colNo == c)
        
        Xmin <- Xmax <- c(colMeans(object$X[,-K]), 1)
        Xmin[c] = min(object$X[,c])
        Xmax[c] = max(object$X[,c])
        
        Xbeta_min_R = matrix(0, nrow=R, ncol=j)
        Xbeta_max_R = matrix(0, nrow=R, ncol=j)
        for (i in 2:j) {
          Xbeta_min_R[, i] = betamat_R[, i, ] %*% Xmin
          Xbeta_max_R[, i] = betamat_R[, i, ] %*% Xmax
        }
        
        yhat_min_R = exp(Xbeta_min_R) / rowSums(exp(Xbeta_min_R))
        yhat_max_R = exp(Xbeta_max_R) / rowSums(exp(Xbeta_max_R))
        
        ydisc_R = yhat_max_R - yhat_min_R # R x j
        for (i in 1:j) {
          se_mat[i, c1] = sd(ydisc_R[, i])
        }
      }
    }
  }
  
  # generating hypothesis testing tables.
  listmat = list()
  if(se){
    for(i in 1:k){
      tabout = matrix(ncol=4,nrow=j)
      tabout[,1:2] = cbind(xmarg[,i],se_mat[,i])
      tabout[,3] = tabout[,1] / tabout[,2]
      tabout[,4] = 2*(1-pnorm(abs(tabout[,3])))
      colnames(tabout) = c("estimate","std","z","p-value")
      rownames(tabout) = ynames
      listmat[[i]] = tabout
    }
    names(listmat)=varlist
  }
  
  colnames(xmarg) <- colnames(se_mat) <- varlist
  rownames(xmarg) <- rownames(se_mat) <-colnames(object$y)
  outlist=list()
  outlist$effects = xmarg
  if(se==TRUE){outlist$se = se_mat; outlist$ztable = listmat}
  if(marg.type=="aveacr") {names(marg_list)=varlist; outlist$marg.list = marg_list}
  marg.type.out = ifelse(marg.type=="atmean","at the mean,","average across observations,")

  # please include this in the file
  outlist$R = ifelse(se,R,0)

  outlist$expl = paste(effect,"effect",marg.type.out,
                       ifelse(se==TRUE,"Krinsky-Robb standard error calculated","standard error not computed"))
  return(structure(outlist,class="fracregmlogit.pe"))
}

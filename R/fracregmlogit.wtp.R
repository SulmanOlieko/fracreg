#' "Willingness to Pay" for fracregmlogit models
#' 
#' Evaluate the "Willingness to Pay" given a set of arbitrary values for outcome variables. 
#' Usually used for policy evaluations where the total magnitude of marginal change matters.
#' 
#' @param object A \code{fracregmlogit.pe} object.
#' @param wtp.vec A numeric vector containing the arbitrary outcome values to be evaluated for each choice j.
#' @param varlist A string vector which provides the names of variables to calculate 
#' the wtp for. If missing, all variables in the object will be calculated.
#' @param indv.obs A logical value indicating whether to return individual observations.
#' @return A "fracregmlogit.wtp" object containing the estimates, standard error, z-stats, and p-value.
#' @details This function calculates the aggregate effect of a variable on the 
#' "willingness to pay" by linearly multiplying the average partial effect with ex-ante (arbitrary) 
#' willingness to pay numbers associated with each choice. 
#' 
#' Suppose there are three choices A, B, C, each with a willingness to pay (or cost, profit, budget),
#' of 100, 200, and 300. The discrete effects of variable X on A, B and C are 0.5, 0.5, and -1, with 
#' standard errors 0.2, 0.3 and 0.5. The aggregated discrete effect of X on the total willingness 
#' to pay (or cost), is thus 100*0.5 + 200*0.5 + 300*(-1) = -150. The standard error can also be
#' calculated to be 162.8, assuming that the standard error is independent. 
#' A simple z-test is provided to test whether the aggregate effect is different from zero. 
#' 
#' Note that if the input \code{fracregmlogit.pe} object has no standard error computation, then no standard error
#' will be computed for the willingness to pay.
#' 
#' @seealso \code{\link{fracregmlogit.pe}}, \code{\link{plot.fracregmlogit}}
#' 
#' @examples
#' \donttest{
#' data("fracreg_spending")
#' X = fracreg_spending[,2:5]
#' y = fracreg_spending[,6:11]
#' results1 = fracregmlogit(y, X)
#' pe = fracregmlogit.pe(results1)
#' # Assuming arbitrary WTP values for the 6 choices
#' wtp_est = wtp(pe, wtp.vec = c(1, 2, 3, 4, 5, 6), varlist = "houseval")
#' summary(wtp_est)
#' }
#' @export wtp
wtp = function(object, wtp.vec, varlist=NULL, indv.obs=FALSE){
  if(!inherits(object, "fracregmlogit.pe")) stop("Expect an fracregmlogit.pe object. Wrong object type given.")
  j=ncol(object$effects); k=nrow(object$effects)
  Xnames = rownames(object$effects); ynames = colnames(object$effects)
  
  if(length(varlist)==0){
    varlist=Xnames
  }else{
    if (!all(varlist %in% Xnames)) stop("Some variables in varlist not found in object.")
  }
  k = length(varlist)
  if(length(wtp.vec)!=j) stop("Wrong length of wtp.vec. Please check specification again.")
  
  # wtp calcs
  betamat = object$effects[varlist, , drop=FALSE]; semat = object$se[varlist, , drop=FALSE]
  wtp_mean = as.vector(betamat %*% wtp.vec)

  if(object$R>0){ # prevent a bug that does not output R in the fracregmlogit.pe module. 
  wtp_se = as.vector(sqrt((semat^2) %*% (wtp.vec^2)))
  # output tables
  tabout = matrix(ncol=4,nrow=k)
  tabout[,1] = wtp_mean
  tabout[,2] = wtp_se
  tabout[,3] = tabout[,1] / tabout[,2]
  tabout[,4] = 2*(1-pnorm(abs(tabout[,3])))
  colnames(tabout) = c("Coefficient", "Std. Error", "z value", "Pr(>|z|)")
  rownames(tabout) = varlist
  }else tabout = wtp_mean
  if(indv.obs){
    if(is.null(object$marg.list)) {
      stop("Individual observations (indv.obs=TRUE) require marg.type='aveacr' in fracregmlogit.pe.")
    }
    wtp_mat = matrix(ncol = k, nrow=nrow(object$marg.list[[1]]))
    for(c1 in 1:length(varlist)){
      v = varlist[c1]
      wtp_mat[,c1] = as.matrix(object$marg.list[[v]]) %*% wtp.vec
    }
    colnames(wtp_mat) = varlist
  }
  # output list
  outlist = list()
  outlist$wtp = tabout
  if(indv.obs) outlist$wtp.obs = wtp_mat
  return(structure(outlist,class="fracregmlogit.wtp"))
}
  

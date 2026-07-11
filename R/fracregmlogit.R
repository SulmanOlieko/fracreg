#' Estimate Fractional Multinomial Logit Models
#' 
#' Used to estimate fractional multinomial logit models using quasi-maximum
#' likelihood estimation following Papke and Wooldridge (1996).
#' @param y the dependent variable (N*J). Can be a matrix or a named data frame.
#'   The first column of the matrix is automatically treated as the baseline.
#' @param X independent variable (N*K). Can be a matrix or a named data frame.
#'   If there is no intercept term in the X, then an intercept term is
#'   automatically added.
#' @param beta0 Initial value for beta used in optimisation. Uses a 1*K(J-1)
#'   vector. Default to a vector of zeros.
#' @param MLEmethod Method of optimisation. Goes into
#'   \code{maxLik(method=MLEmethod))}. Choose from "NR","BFGS","CG","BHHH","SANN",or "NM".  
#'   Default to "CG", the conjugate gradients method. See Details. 
#' @param maxit Maximum number of iterations.
#' @param abstol Tolerance.
#' @param cluster A vector of clusters to be used for clustered standard error computation. 
#' Default to NULL, no cluster computed. 
#' @param reps Number of bootstrap replications to be computed for clustered standard errors.
#' @param ... additional parameters that go into \code{maxLik()}
#' @return The function returns an object of class "fracregmlogit". Use \code{fracregmlogit.pe}, \code{predict}, 
#'  \code{residuals}, \code{fitted} to extract various useful features of the value returned by 
#' \code{fracregmlogit}. 
#' @return An object of class "fracregmlogit" contains the following components: 
#' @return \code{estimates}   A list of matrices containing parameter estimates,
#'   standard errors, and hypothesis testing results.
#' @return \code{baseline}    The baseline choice
#' @return \code{likelihood}  The likelihood value
#' @return \code{conv_code}   Convergence diagnostics code. 
#' @return \code{convergence} Convergence messages. 
#' @return \code{count}       Provides dataset information
#' @return \code{y}           The dependent variable data frame.
#' @return \code{X}           The independent variable data frame. Augmented by factor dummy transformation
#' , constant term added. 
#' @return \code{rowNo}       A vector of row numbers from the original X and y that is used for estimation.
#' @return \code{coefficient} Matrix of estimated coefficients. Augmented with the baseline coefficient
#' (which is a vector of zeros). 
#' @return \code{vcov}        A list of matrices containing the robust variance covariance matrix for each choice
#' variable. 
#' @return \code{cluster}     The vector of clusters. 
#' @return \code{reps}        Number of bootstrap replications for clustered standard error
#' @details The fractional multinomial logit model is the expansion of the multinomial
#' logit to fractional responses. Unlike standard multinomial logit models,
#' which only consider 0-1 responses, the fractional multinomial logit model considers the
#' case where the response variable is fractions that sum up to one. Examples
#' of this type of data include percentages of budget spent in education, defence,
#' public health; fractions of a population that have middle school, high
#' school, college, or post-college education, etc.
#' 
#' This function follows Papke and Wooldridge (1996)'s paper, in which they
#' proposed a quasi-maximum likelihood estimator for fractional response data.
#' The likelihood function used here is a standard multinomial likelihood
#' function, see Buis (2008) and <http://maartenbuis.nl/software/likelihoodFmlogit.pdf> for
#' the likelihood used here. Robust standard errors are provided following Papke
#' and Wooldridge (1996), in which they proposed an asymptotically consistent
#' estimator of variance.
#' 
#' Maximisation is done by calling \code{\link[maxLik]{maxLik}}. maxLik is a wrapper function 
#' for different maximisation methods in R. These include most methods provided by \code{\link[maxLik]{maxLik}},  
#' but also other methods such as BHHH (Berndt-Hall-Hall-Hausman). 
#' 
#' MLE convergence can be a problem in R, especially if the dataset is large with many explanatory variables. 
#' It is recommended to call CG (Conjugate Gradients) or BHHH (Berndt-Hall-Hall-Hausman).
#' The conjugate gradients method is usually faster, but could lead to non-convergence under 
#' certain scenarios. BHHH is slower, but has better convergence properties.
#' 
#' @seealso \code{\link{fracregmlogit.pe}} for computing partial effects, \code{\link{plot.fracregmlogit.pe}} for plotting effects, \code{\link{fitted.fracregmlogit}} for residuals and predictions.
#' 
#' @examples
#' data("fracreg_spending")
#' X = fracreg_spending[,2:5]
#' y = fracreg_spending[,6:11]
#' 
#' # Fit the fractional multinomial logit model
#' results1 = fracregmlogit(y, X)
#' 
#' # View estimates
#' summary(results1)
#' 
#' # Compute marginal effects
#' pe = fracregmlogit.pe(results1, effect="marginal", marg.type="aveacr", se=TRUE, R=50)
#' summary(pe)
#' 
#' # Plot effects for 'houseval'
#' plot(pe, varlist="houseval")
#' 
#' @references Papke, L. E. and Wooldridge, J. M. (1996), Econometric methods
#'   for fractional response variables with an application to 401(k) plan
#'   participation rates. J. Appl. Econ., 11: 619-632.
#' @references Buis, M. L. (2008), fmlogit: Stata module fitting a fractional 
#'   multinomial logit model by quasi maximum likelihood. Statistical Software Components, 
#'   Boston College Department of Economics.
#' @references Mullahy, J. (2015), Multivariate fractional regression estimation of 
#'   econometric share models. Journal of Econometric Methods, 4(1): 71-100.
#' @references Murteira, J. M. R., and Ramalho, J. J. S. (2016), Regression analysis of 
#'   multivariate fractional data. Econometric Reviews, 35(4): 515-552.
#' @references Ji, J., and Woodill, A. J., fmlogit: Fractional Multinomial Logit. R package repository. <https://github.com/f1kidd/fmlogit>.
#' @importFrom maxLik maxLik
#' @export fracregmlogit




fracregmlogit=function(y, X, beta0 = NULL, MLEmethod = "CG", maxit = 5e+05, 
                          abstol = 1e-05,cluster=NULL,reps=1000, ...){
  start.time = proc.time()
  
  if(length(cluster)!=nrow(y) & !is.null(cluster)){
    warning("Length of the cluster does not match the data. Cluster is ignored.")
    cluster = NULL
  }
  Xclass = sapply(X, class)
  Xfac = which(Xclass %in% c("factor", "character"))
  if (length(Xfac) > 0) {
    Xfacnames = colnames(X)[Xfac]
    strformFac = paste(Xfacnames, collapse = "+")
    Xdum = model.matrix(as.formula(paste("~", strformFac, 
                                         sep = "")), data = X)[, -1]
    X = cbind(X, Xdum)
    X = X[, -Xfac]
  }
  Xnames = colnames(X)
  ynames = colnames(y)
  X = as.matrix(X)
  y = as.matrix(y)
  n = dim(X)[1]
  j = dim(y)[2]
  k = dim(X)[2]
  xy = cbind(X, y)
  xy = na.omit(xy)
  row.remain = setdiff(1:n, attr(xy, "na.action"))
  X = xy[, 1:k]
  y = xy[, (k + 1):(k + j)]
  n = dim(y)[1]
  remove(xy)
  # adding in the constant term
  if(k==1){
    # check if the input X is constant
    if(length(unique(X))==1){ # X is constant
      Xnames = "constant"
      X = as.matrix(as.numeric(X),nrow=1)
      colnames(X) = Xnames
      X = as.matrix(X)
      k=0
    }else{ # one single variable of input
      Xnames = "X1"
      X = as.matrix(X)
      k = dim(X)[2]
      X = cbind(X, rep(1, n))
      Xnames = c(Xnames, "constant")
      colnames(X) = Xnames
    }
  }else{ # normal cases
    X = X[, apply(X, 2, function(x) length(unique(x)) != 1)]
    Xnames = colnames(X)
    k = dim(X)[2]
    X = cbind(X, rep(1, n))
    Xnames = c(Xnames, "constant")
    colnames(X) = Xnames
  }
  
  
  testcols <- function(X) {
    m = crossprod(as.matrix(X))
    ee = eigen(m)
    evecs <- split(zapsmall(ee$vectors), col(ee$vectors))
    mapply(function(val, vec) {
      if (val != 0) 
        NULL
      else which(vec != 0)
    }, zapsmall(ee$values), evecs)
  }
  collinear = unique(unlist(testcols(X)))
  while (length(collinear) > 0) {
    if (qr(X)$rank == dim(X)[2]) 
      print("Model may suffer from multicollinearity problems.")
    break
    if ((k + 1) %in% collinear) 
      collinear = collinear[-length(collinear)]
    X = X[, -collinear[length(collinear)]]
    Xnames = colnames(X)
    k = k - 1
    collinear = unique(unlist(testcols(X)))
  }

  QMLE_Obs <- function(betas) {
    betas = matrix(betas, nrow = j - 1, byrow = T)
    betamat = rbind(rep(0, k + 1), betas)
    
    # Pre-calculate to save time
    Xbeta = X %*% t(betamat)
    # log_sum_exp vector of length n
    log_sum_exp = log(rowSums(exp(Xbeta)))
    
    # L is n x j matrix
    L = y * (Xbeta - log_sum_exp)
    
    # log likelihood for each observation
    llf = rowSums(L)
    return(llf)
  }
  
  grad_Obs <- function(betas) {
    betas = matrix(betas, nrow = j - 1, byrow = T)
    betamat = rbind(rep(0, k + 1), betas)
    
    Xbeta = X %*% t(betamat)
    # G is n x j matrix
    G = exp(Xbeta)
    sum_G = rowSums(G)
    G = G / sum_G
    
    # gradient for each observation n and parameter i (for i = 2..j)
    grad = matrix(0, nrow = n, ncol = (j - 1) * (k + 1))
    
    for (i in 2:j) {
      grad[, ((i-2)*(k+1) + 1):((i-1)*(k+1))] = X * (y[, i] - G[, i])
    }
    
    return(grad)
  }

  if (length(beta0) == 0){
    beta0 = rep(0, (k + 1) * (j - 1))
  }
  if (length(beta0) != (k + 1) * (j - 1)) {
    beta0 = rep(0, (k + 1) * (j - 1))
    warning("Wrong length of beta0 given. Use default setting instead.")
  }
  opt <- maxLik::maxLik(QMLE_Obs, grad = grad_Obs, start = beta0, method = MLEmethod, 
                control = list(iterlim = maxit, tol = abstol), ...)
  betamat = matrix(opt$estimate, ncol = k + 1, byrow = T)
  betamat_aug = rbind(rep(0, k + 1), betamat)
  colnames(betamat_aug) = Xnames
  rownames(betamat_aug) = ynames
  sigmat = matrix(nrow = j - 1, ncol = k + 1)
  vcov = list()
  
  if(is.null(cluster)==F){
    cluster = cluster[row.remain]
    clusters <- names(table(cluster))
    for (i in 1:j) {
      sterrs <- matrix(NA, nrow=reps, ncol=k + 1)
      vcov_j_list=list()
      
      b=1
      no_singular_error=c()
      while(b<=reps){
        index <- sample(1:length(clusters), length(clusters), replace=TRUE)
        aa <- clusters[index]
        bb <- table(aa)
        bootdat <- NULL
        dat=cbind(y,X)
        for(b1 in 1:max(bb)){
          cc <- dat[cluster %in% names(bb[bb %in% b1]),]
          for(b2 in 1:b1){
            bootdat <- rbind(bootdat, cc)
          }
        }
        
        bootdatX=matrix(bootdat[,(j+1):ncol(bootdat)],nrow=nrow(bootdat))
        bootdaty=bootdat[,1:j]
        
        sum_expxb = rowSums(exp(bootdatX %*% t(betamat_aug)))
        expxb = exp(bootdatX %*% betamat_aug[i, ])
        G = expxb/sum_expxb
        g = (expxb * sum_expxb - expxb^2)/sum_expxb^2
        X_a = bootdatX * as.vector(sqrt(g^2/(G * (1 - G))))
        A = t(X_a) %*% X_a
        mu = bootdaty[, i] - G
        X_b = bootdatX * as.vector(mu * g/G/(1 - G))
        B = t(X_b) %*% X_b
        
        a_solve_error = tryCatch(solve(A),error=function(e){NULL})
        if(is.null(a_solve_error)){
          no_singular_error=c(no_singular_error,b)
          next
        }
        
        Var_b = solve(A) %*% B %*% solve(A)
        std_b = sqrt(diag(Var_b))
        sterrs[b,]=std_b
        vcov_j_list[[b]]=Var_b
        
        b=b+1
      }
      if(length(no_singular_error)>0){warning(paste('Error in solve.default(A) : Lapack routine dgesv: system is exactly singular: U[28,28] = 0" Appeared',length(no_singular_error),'times within cluster bootstrap for outcome #',i))}
      std_b=apply(sterrs,2,mean)
      vcov[[i]] = Reduce("+", vcov_j_list) / length(vcov_j_list)
      if (i > 1) 
        sigmat[i - 1, ] = std_b
    }
  }else{
    for(i in 1:j){
      sum_expxb = rowSums(exp(X %*% t(betamat_aug))) 
      expxb = exp(X %*% betamat_aug[i,]) 
      G = expxb / sum_expxb 
      g = (expxb * sum_expxb - expxb^2) / sum_expxb^2 
      
      X_a = X * as.vector(sqrt(g^2/(G*(1-G))))
      A = t(X_a) %*% X_a
      
      mu = y[,i] - G
      X_b = X * as.vector(mu * g / G / (1-G))
      B = t(X_b) %*% X_b
      Var_b = solve(A) %*% B %*% solve(A)
      std_b = sqrt(diag(Var_b))
      vcov[[i]] = Var_b
      if(i>1) sigmat[i-1,] = std_b
    }
  }
  
  listmat = list()
  for (i in 1:(j - 1)) {
    tabout = matrix(ncol = 4, nrow = k + 1)
    tabout[, 1:2] = t(rbind(betamat[i, ], sigmat[i, ]))
    tabout[, 3] = tabout[, 1]/tabout[, 2]
    tabout[, 4] = 2 * (1 - pnorm(abs(tabout[, 3])))
    colnames(tabout) = c("estimate", "std", "z", "p-value")
    if (length(Xnames) > 0) 
      rownames(tabout) = Xnames
    listmat[[i]] = tabout
  }
  if (length(ynames) > 0) 
    names(listmat) = ynames[2:j]
  outlist = list()
  outlist$estimates = listmat
  outlist$baseline = ynames[1]
  outlist$likelihood = opt$maximum
  outlist$conv_code = opt$code
  outlist$convergence = paste(opt$type, paste(as.character(opt$iterations), 
                                              "iterations"), opt$message, sep = ",")
  outlist$count = c(Obs = n, Explanatories = k, Choices = j)
  outlist$y = y
  outlist$X = X
  outlist$rowNo = row.remain
  outlist$coefficient = betamat_aug
  names(vcov) = ynames
  outlist$vcov = vcov
  outlist$cluster = cluster
  outlist$reps=ifelse(is.null(cluster),0,reps)
  
  print(paste("Fractional logit model estimation completed. Time:", 
              round(proc.time()[3] - start.time[3], 1), "seconds"))
  return(structure(outlist, class = "fracregmlogit"))
}








#' Plot Marginal or Discrete Effects of Willingness to Pay
#' 
#' Plot marginal or discrete effects of willingness to pay, potentially against another variable.
#' 
#' @param x A "fracregmlogit" object.
#' @param wtp.vec A numeric vector for willingness to pay.
#' @param varlist A string vector which provides the names of variables to plot the effect for.
#'  If missing, all variables in the object will be plotted.
#' @param against A vector with the same length as the number of observations in the model, or the name of a variable. 
#' Serves as the x-axis in the plots.
#' @param mfrow A numeric vector with two elements. Specifies the number of rows and columns in a panel.
#' Similar to par(mfrow=c()). Default to NULL, and the program will choose a square panel. 
#' @param t Number of points to be used for smoothing.
#' @param effect The type of effect ("marginal" or "discrete").
#' @param type Plot type.
#' @param plot.show If TRUE, the plot will be created. Otherwise, the function returns raw data that can be
#' used to create user-specified (custom) plots. 
#' @param ... Additional arguments.
#' @return Panel plots of effects vs. chosen variables.
#' @details 
#' This function provides a visualisation tool for potentially heterogeneous marginal and discrete effects of willingness to pay.
#' The function allows the user to plot marginal effects to detect any patterns in the effects, in itself
#' and against other variables. The plot also allows visualisation of sub-groups in data, which can be
#' very useful to visualise categorical and dummy variables. 
#' 
#' The function takes a \code{fracregmlogit} object, and internally calls \code{wtp} and \code{fracregmlogit.pe} 
#' to compute the willingness to pay at different data points.
#' 
#' Additional parameters include \code{varlist}, a vector of string variable names to be plotted. 
#'  
#' \code{against} allows a different variable to be chosen
#' as the x-axis. \code{against} can supply the column name of a variable in the original dataset to plot against. 
#'  
#' @seealso \code{\link{wtp}}, \code{\link{fracregmlogit.pe}}
#' @examples
#' data("fracreg_spending")
#' X = fracreg_spending[,2:5]
#' y = fracreg_spending[,6:11]
#' results1 = fracregmlogit(y, X)
#' 
#' # Define a willingness to pay vector
#' wtp.vec = c(1, 1, 1, 1, 1)
#' 
#' # Plot WTP for 'popdens'
#' plot(results1, wtp.vec=wtp.vec, varlist="popdens")
#' @export plot.fracregmlogit


plot.fracregmlogit = function(x, wtp.vec=NULL, varlist=NULL, against=NULL, mfrow=NULL, t=500, effect=c("discrete","marginal"), type=NULL, plot.show=TRUE, ...){
  effect = match.arg(effect)
  object <- x
  K = ncol(object$X); j = ncol(object$y); N = nrow(object$X); 
  Xnames = colnames(object$X) ; ynames = colnames(object$y)
  X = object$X; y=object$y
  
  # determine variable list
  var_colNo = which(Xnames %in% varlist)
  k = length(var_colNo)
  
  if(is.null(mfrow)){
    js = ceiling(sqrt(k))
    jr = ifelse(js*(js-1)>=k,js-1,js)
  }else{
    jr = mfrow[1]; js = mfrow[2]
  }
  
  if(!is.null(against)) {
    ag_No = which(Xnames == against)
    if(length(ag_No)==0) stop(paste("The against vector specified,",against,
                                    "is not in the list of explanatory variables. Please check again."))
    ag_min = min(X[,ag_No]); ag_max = max(X[,ag_No])
    ag_vec = seq(ag_min,ag_max,length.out = t)
    wtp_mat = matrix(nrow=t,ncol=k)
    colnames(wtp_mat) = varlist
    for(i in 1:t){
      newdata = colMeans(X[,-K])
      newdata[ag_No] = ag_vec[i]
      wtp_mat[i,] = wtp(fracregmlogit.pe(object,effect=effect,se=F,varlist=varlist,at=newdata),wtp.vec)[[1]]
    }
  }else{
    against="ObsNo"
    ag_vec=1:N
    wtp_mat = matrix(nrow=N,ncol=k)
    colnames(wtp_mat) = varlist
    for(i in 1:N){
      newdata = X[i,-K]
      wtp_mat[i,] = wtp(fracregmlogit.pe(object,effect=effect,se=F,varlist=varlist,at=newdata),wtp.vec)[[1]]
    }
  }
  # plotting
  if(plot.show){
    par(mfrow=c(jr,js))
    if(is.null(type)){type="l"} # default to line plot. 
    for(i in 1:k){
      plot(ag_vec,wtp_mat[,i],xlab=against,ylab=paste(effect,"effect of", varlist[i]),...)
    }}
  return(list(ag_vec,wtp_mat))
}
#' Plot Marginal or Discrete Effects
#' 
#' Plot the desired effect at each observed value for each choice.
#' 
#' @param x A "fracregmlogit.pe" object.
#' @param X A matrix of independent variables.
#' @param y A matrix of dependent variables.
#' @param varlist A string vector which provides the names of variables to plot the effect for.
#'  If missing, all variables in the object will be plotted.
#' @param against A vector with the same length as the number of observations in the model. 
#' Serves as the x-axis in the plots.
#' @param against.x A character string, supply the column name in the X matrix to plot against.
#' @param against.y A character string, supply the column name in the y matrix to plot against.
#' @param group.x A character string. Supply the column name in the X matrix to group upon. 
#' @param group.algebra A character string. Supply additional algebra imposed on the group variable. 
#' @param mfrow A numeric vector with two elements. Specifies the number of rows and columns in a panel.
#' Similar to par(mfrow=c()). Default to NULL, and the program will choose a square panel. 
#' @param ... Additional arguments.
#' @return Panel plots of effects vs. chosen variables.
#' @details 
#' This function provides a visualisation tool for potentially heterogeneous marginal and discrete effects.
#' The function allows the user to plot marginal effects to detect any patterns in the effects, in itself
#' and against other variables. The plot also allows visualisation of sub-groups in data, which can be
#' very useful to visualise categorical and dummy variables. 
#' 
#' The function takes a \code{fracregmlogit.pe} object, created by the \code{fracregmlogit.pe()} function. Note that since 
#' the plotting requires marginal effects for all observations, the object should be created by choosing 
#' \code{marg.type="aveacr"}, the average across method for effects calculation. 
#' 
#' Additional parameters include \code{varlist}, a vector of string variable names to be plotted. \code{X}
#' and \code{y} are the dependent and independent variable matrices in the original regression model. 
#'  
#' \code{against}, \code{against.x}, and \code{against.y} allow different variables to be chosen
#' as the x-axis. \code{against} directly supplies the vector to be plotted against, whereas \code{against.x}
#' and \code{against.y} supply variable names in the original dataset. Note that the user has to provide
#' \code{X} and \code{y} in order to use the column name options, respectively. 
#'  
#' \code{group.x} supplies the column name in the X matrix to group by. The plot will be able to 
#' differentiate different groups by colours. Additionally, the user can supply a string to \code{group.algebra},
#' which provides an algebra operation that will be evaluated on the group vector. For example, choosing 
#' \code{group.x = "a"} and \code{group.algebra = ">0"} will create two groups, one with X$a > 0, and one with X$a <= 0.
#' 
#' @seealso \code{\link{fracregmlogit.pe}}
#' 
#' @examples
#' data("fracreg_spending")
#' X = fracreg_spending[,2:5]
#' y = fracreg_spending[,6:11]
#' results1 = fracregmlogit(y, X)
#' 
#' # Calculate marginal effects with marg.type="aveacr" (no standard errors for speed)
#' effect1 = fracregmlogit.pe(results1, effect="marginal", marg.type="aveacr", se=FALSE)
#' 
#' # Plot effects
#' plot(effect1, X=results1$X, against.x = "houseval", group.x = "popdens", group.algebra = ">10")
#' @export plot.fracregmlogit.pe



plot.fracregmlogit.pe = function(x, varlist=NULL, X=NULL, y=NULL, 
                                against=NULL,against.x=NULL,against.y=NULL,
                                group.x=NULL, group.algebra=NULL,
                                mfrow=NULL, ...){
  object <- x
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("grid", quietly = TRUE)
  
  if(is.null(object[["marg.list"]])) stop("Please choose marg.type=aveacr when calculating effects")
  k = ncol(object$effects); j = nrow(object$effects); N = nrow(object$marg.list[[1]]); 
  Xnames = colnames(object$effects) ; ynames = rownames(object$effects)
  # X = object$X; y=object$y
  
  # determine variable list
  if(length(varlist)==0){
    varlist=Xnames
    var_colNo = 1:k
  }else{
    var_colNo = which(Xnames %in% varlist)
    k = length(var_colNo)
  }
  if(k==0) stop("Variable list not matched. Please check your varlist input.")
  
  # determine panel size
  if(is.null(mfrow)){
    js = ceiling(sqrt(j))
    jr = ifelse(js*(js-1)>j,js-1,js)
  }else{
    jr = mfrow[1]; js = mfrow[2]
  }
  
  # determine plotting x axis. 
  if(is.null(against) & is.null(against.x) & is.null(against.y)){
    M.against=1:N 
    ag.name = "ObsNo"
  }else if(is.null(against.x)==F){
    M.against = X[,against.x]
    if(is.null(M.against)){
      stop("against.x not found in variable list. Please double check your spelling")
    }
    ag.name = against.x
  }else if(is.null(against.y)==F){
    M.against = y[,against.y] 
    ag.name = against.y
  }else{M.against=against}
  
  
  # determine group variables
  if(is.null(group.x) & is.null(group.algebra)) {M.group=NULL; g.name=NULL}
  if(is.null(group.x)==F) {M.group = X[,group.x]; g.name.display <- g.name <- group.x;}
  if(is.null(group.algebra)==F) {
    M.group = eval(parse(text=paste("X[,",'"',group.x,'"',"]",group.algebra,sep="")))
    M.group = ifelse(M.group,"Yes","No")
    g.name = group.x
    g.name.display = paste(group.x,group.algebra,sep="")
    }
  
  for(c in var_colNo){
    ggplot()
    pushViewport(viewport(layout = grid.layout(jr, js)))
    temp.data = cbind(object$marg.list[[c]],M.against)
    temp.data = as.data.frame(temp.data)
    colnames(temp.data) = c(colnames(object$marg.list[[c]]),ag.name)
    if(is.null(M.group)==F){
      temp.data = cbind(temp.data,as.factor(M.group))
      colnames(temp.data)[ncol(temp.data)] = g.name
    }
    for(i in 1:j){
      g <- ggplot(temp.data,aes_string(ag.name,ynames[i],color=g.name)) + geom_point() 
      g <- g + geom_hline(yintercept = 0) + theme_classic() + ggtitle(paste("Effects on", Xnames[c]))
      if(is.null(M.group)==F) g <- g + theme(legend.title = element_text(colour="black"))+
        scale_color_discrete(name=g.name.display)
      print(g,vp = viewport(layout.pos.row = ifelse(i%%jr==0,jr,i%%jr), layout.pos.col = (i-1) %/%js + 1) )
  }
}}
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
#' @seealso \code{\link{fracregmlogit}}
#' @examples
#' data("fracreg_spending")
#' X = fracreg_spending[,2:5]
#' y = fracreg_spending[,6:11]
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
#' @export fitted.fracregmlogit
#' 


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
#' @export residuals.fracregmlogit
#' 
residuals.fracregmlogit <- function(object, ...) {
  yhat = fitted(object)
  return(as.data.frame(object$y-yhat))
}

#' @rdname fitted.fracregmlogit
#' @export predict.fracregmlogit
#' 
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


#' Generate Summary Tables for fracregmlogit Objects
#' 
#' Generate tables of coefficient estimates, partial effects, and willingness to pay from
#' fracregmlogit-type objects. 
#' 
#' @name summary.fracregmlogit
#' @aliases summary.fracregmlogit.pe
#' @aliases summary.fracregmlogit.wtp
#' 
#' @param object an object with class "fracregmlogit", "fracregmlogit.pe", or "fracregmlogit.wtp".
#' @param ... Additional arguments passed to the printCoefmat function.
#' @return Returns the object invisibly.
#' 
#' @details This module provides summary methods for three fracregmlogit objects: \code{fracregmlogit}, \code{fracregmlogit.pe}
#' , and \code{fracregmlogit.wtp}. 
#' 
#' For \code{fracregmlogit} objects, the summary prints the number of observations, log pseudo-likelihood,
#' baseline choice, and the coefficient estimates with standard errors, z-statistics, and p-values
#' for each choice equation. 
#' 
#' For \code{fracregmlogit.pe} objects, it displays the marginal or discrete effects 
#' along with their computed standard errors (if Krinsky-Robb sampling was performed) for each choice.
#' 
#' For \code{fracregmlogit.wtp} objects, it provides a table of the aggregated willingness to pay 
#' along with its standard errors and test statistics.
#' @seealso \code{\link{fracregmlogit}}, \code{\link{fracregmlogit.pe}}
#' @examples
#' data("fracreg_spending")
#' X = fracreg_spending[,2:5]
#' y = fracreg_spending[,6:11]
#' 
#' # generate fracregmlogit summary
#' results1 = fracregmlogit(y, X)
#' summary(results1)
#' 
#' # generate marginal effects summary
#' effects1 = fracregmlogit.pe(results1, effect="marginal", se=FALSE)
#' summary(effects1)
#' 
#' @rdname summary.fracregmlogit
#' @export summary.fracregmlogit
#' @export summary.fracregmlogit.pe
#' @export summary.fracregmlogit.wtp

############
# generate willingness to pay tables
############

summary.fracregmlogit.wtp = function(object, ...) {
  if(!inherits(object, "fracregmlogit.wtp")) stop("Expect an fracregmlogit.wtp object. Wrong object type given.")
  if (is.null(dim(object$wtp))) return(object$wtp)
  if(is.null(colnames(object$wtp)) || colnames(object$wtp)[1]!="estimate") return(object$wtp) # no need to summary.
  
  cat("\nFractional Multinomial Logit Model - Willingness to Pay\n")
  stats::printCoefmat(object$wtp, digits = max(3, getOption("digits") - 2), 
                      signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n")
  invisible(object)
}
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
  j=nrow(object$effects); k=ncol(object$effects)
  Xnames = colnames(object$effects); ynames = rownames(object$effects)
  if(length(varlist)==0){
    varlist=Xnames
    var_colNo = c(1:k)
    k = length(var_colNo)
  }else{
    var_colNo = which(varlist %in% Xnames)
    k = length(var_colNo)
  }
  if(length(wtp.vec)!=j) stop("Wrong length of wtp.vec. Please check specification again.")
  # wtp calcs
  betamat = object$effects[,varlist, drop=FALSE]; semat = object$se[,varlist, drop=FALSE]
  wtp_mean = wtp.vec %*% betamat

  if(object$R>0){ # prevent a bug that does not output R in the fracregmlogit.pe module. 
  wtp_se = sqrt(wtp.vec^2 %*% semat^2)
  # output tables
  tabout = matrix(ncol=4,nrow=k)
  tabout[,1] = wtp_mean
  tabout[,2] = wtp_se
  tabout[,3] = tabout[,1] / tabout[,2]
  tabout[,4] = 2*(1-pnorm(abs(tabout[,3])))
  colnames(tabout) = c("estimate","std","z","p-value")
  rownames(tabout) = varlist
  }else tabout = wtp_mean
  if(indv.obs){
    wtp_mat = matrix(ncol = k, nrow=nrow(object$marg.list[[1]]))
    for(c in var_colNo){
      c1 = which(var_colNo == c)
      wtp_mat[,c1] = as.matrix(object$marg.list[[c1]]) %*% wtp.vec
    }
    colnames(wtp_mat) = varlist
  }
  # output list
  outlist = list()
  outlist$wtp = tabout
  if(indv.obs) outlist$wtp.obs = wtp_mat
  return(structure(outlist,class="fracregmlogit.wtp"))
}
  

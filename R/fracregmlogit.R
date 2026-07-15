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
#' @param na.action A function specifying how to handle missing values, default is \code{stats::na.omit}. If \code{NULL}, no action is taken.
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
                          abstol = 1e-05,cluster=NULL,reps=1000, na.action=stats::na.omit, ...){
  start.time = proc.time()
  
  if (!missing(y) && !missing(X)) {
      args_to_clean <- list(y=y, X=X)
      if (!is.null(cluster)) args_to_clean$cluster <- cluster
      
      args_to_clean$na.action <- na.action
      cleaned <- do.call(fracreg_clean_data, args_to_clean)
      
      y <- cleaned$y
      X <- cleaned$X
      if (!is.null(cluster)) cluster <- cleaned$cluster
  }
  
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
    betas = matrix(betas, nrow = j - 1, byrow = TRUE)
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
    betas = matrix(betas, nrow = j - 1, byrow = TRUE)
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
  betamat = matrix(opt$estimate, ncol = k + 1, byrow = TRUE)
  betamat_aug = rbind(rep(0, k + 1), betamat)
  colnames(betamat_aug) = Xnames
  rownames(betamat_aug) = ynames
  sigmat = matrix(nrow = j - 1, ncol = k + 1)
  vcov = list()
  
  if(is.null(cluster)==FALSE){
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
  wald_list = list()
  for (i in 1:(j - 1)) {
    tabout = matrix(ncol = 4, nrow = k + 1)
    tabout[, 1:2] = t(rbind(betamat[i, ], sigmat[i, ]))
    tabout[, 3] = tabout[, 1]/tabout[, 2]
    tabout[, 4] = 2 * (1 - pnorm(abs(tabout[, 3])))
    colnames(tabout) = c("estimate", "std", "z", "p-value")
    if (length(Xnames) > 0) 
      rownames(tabout) = Xnames
    listmat[[i]] = tabout
    
    if (k > 0) {
      beta_slope = betamat[i, 1:k]
      V_slope = vcov[[i + 1]][1:k, 1:k]
      W = tryCatch(as.numeric(t(beta_slope) %*% solve(V_slope) %*% beta_slope), error = function(e) NA)
      p_W = if(!is.na(W)) 1 - pchisq(W, df = k) else NA
      wald_list[[i]] = list(W = W, p = p_W, df = k)
    } else {
      wald_list[[i]] = list(W = NA, p = NA, df = 0)
    }
  }
  
  y_bar = colMeans(y)
  LL0 = n * sum(y_bar * log(y_bar + 1e-16))
  pseudo_R2 = 1 - opt$maximum / LL0

  if (length(ynames) > 0) {
    names(listmat) = ynames[2:j]
    names(wald_list) = ynames[2:j]
  }
  outlist = list()
  outlist$estimates = listmat
  outlist$wald = wald_list
  outlist$pseudo_R2 = pseudo_R2
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









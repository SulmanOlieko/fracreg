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
#' wtp.vec = c(1, 1, 1, 1, 1, 1)
#' 
#' # Plot WTP for 'popdens'
#' plot(results1, wtp.vec=wtp.vec, varlist="popdens")
#' @exportS3Method plot fracregmlogit
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
      wtp_mat[i,] = wtp(fracregmlogit.pe(object,effect=effect,se=FALSE,varlist=varlist,at=newdata),wtp.vec)[[1]]
    }
  }else{
    against="ObsNo"
    ag_vec=1:N
    wtp_mat = matrix(nrow=N,ncol=k)
    colnames(wtp_mat) = varlist
    for(i in 1:N){
      newdata = X[i,-K]
      wtp_mat[i,] = wtp(fracregmlogit.pe(object,effect=effect,se=FALSE,varlist=varlist,at=newdata),wtp.vec)[[1]]
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
#' \donttest{
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
#' }
#' @exportS3Method plot fracregmlogit.pe
plot.fracregmlogit.pe = function(x, varlist=NULL, X=NULL, y=NULL, 
                                against=NULL,against.x=NULL,against.y=NULL,
                                group.x=NULL, group.algebra=NULL,
                                mfrow=NULL, ...){
  object <- x
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("grid", quietly = TRUE)
  
  if(is.null(object[["marg.list"]])) stop("Please choose marg.type=aveacr when calculating effects")
  k = nrow(object$effects); j = ncol(object$effects); N = nrow(object$marg.list[[1]]); 
  Xnames = rownames(object$effects) ; ynames = colnames(object$effects)
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
  }else if(is.null(against.x)==FALSE){
    M.against = X[,against.x]
    if(is.null(M.against)){
      stop("against.x not found in variable list. Please double check your spelling")
    }
    ag.name = against.x
  }else if(is.null(against.y)==FALSE){
    M.against = y[,against.y] 
    ag.name = against.y
  }else{M.against=against}
  
  
  # determine group variables
  if(is.null(group.x) & is.null(group.algebra)) {M.group=NULL; g.name=NULL}
  if(is.null(group.x)==FALSE) {M.group = X[,group.x]; g.name.display <- g.name <- group.x;}
  if(is.null(group.algebra)==FALSE) {
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
    if(is.null(M.group)==FALSE){
      temp.data = cbind(temp.data,as.factor(M.group))
      colnames(temp.data)[ncol(temp.data)] = g.name
    }
    for(i in 1:j){
      g <- ggplot(temp.data,aes_string(ag.name,ynames[i],color=g.name)) + geom_point() 
      g <- g + geom_hline(yintercept = 0) + theme_classic() + ggtitle(paste("Effects on", Xnames[c]))
      if(is.null(M.group)==FALSE) g <- g + theme(legend.title = element_text(colour="black"))+
        scale_color_discrete(name=g.name.display)
      print(g,vp = viewport(layout.pos.row = ifelse(i%%jr==0,jr,i%%jr), layout.pos.col = (i-1) %/%js + 1) )
  }
}}

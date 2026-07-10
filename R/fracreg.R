fracreg <- function(y,x,x2=x,linkbin,linkfrac,type="1P",inflation=0,intercept=TRUE,table=FALSE,variance=TRUE,var.type="default",var.eim=TRUE,var.cluster,dfc=FALSE,offset=NULL,or=FALSE,level=0.95,...)
{
	cl <- match.call()
	### 1. Error and warning messages

	if(missing(y)) stop("dependent variable is missing")
	if(missing(x)) stop("explanatory variables are missing")

	if(all(type!=c("1P","2Pbin","2Pfrac","2P","3P"))) stop(sQuote(type)," - type not recognised")
	if(any(y>1) | any(y<0)) stop("The dependent variable has values outside the unit interval")
	if(all(inflation!=c(0,1))) stop(inflation," - value not recognised for inflation")
	if(length(inflation)>1) stop(inflation," - only a single value allowed for inflation")
	if(type!="1P" & !any(y==inflation)) stop("The dependent variable has no ",sQuote(inflation)," values")

	if(type=="2Pbin")
	{
		if(missing(linkbin)) stop("linkbin is missing")
		if(!missing(linkfrac)) warning(sQuote(type)," and",sQuote(linkfrac)," - type does not use linkfrac")

		if(all(linkbin[1]!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkbin[1])," - linkbin not recognised")
	}
	if(type=="3P")
	{
		if(missing(linkbin)) stop("linkbin is missing")
		if(missing(linkfrac)) stop("linkfrac is missing")
		if(length(linkbin) == 1) linkbin <- c(linkbin, linkbin)
		if(all(linkbin[1]!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkbin[1])," - linkbin not recognised")
		if(all(linkbin[2]!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkbin[2])," - linkbin not recognised")
		if(all(linkfrac!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkfrac)," - linkfrac not recognised")
	}
	if(type=="2P")
	{
		if(missing(linkbin)) stop("linkbin is missing")
		if(missing(linkfrac)) stop("linkfrac is missing")
		if(all(linkbin[1]!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkbin[1])," - linkbin not recognised")
		if(all(linkfrac!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkfrac)," - linkfrac not recognised")
	}
	if(any(type==c("1P","2Pfrac")))
	{
		if(missing(linkfrac)) stop("linkfrac is missing")
		if(!missing(linkbin)) warning(sQuote(type)," and",sQuote(linkbin)," - type does not use linkbin")

		if(all(linkfrac!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkfrac)," - linkfrac not recognised")
	}

	if(table==TRUE & variance==FALSE)
	{
		variance <- TRUE
		warning("option variance changed from F to TRUE, as required by table=T")
	}

	if(all(var.type!=c("standard","robust","cluster","default"))) stop(sQuote(var.type)," - var.type not recognised")
	if(var.type=="cluster" & missing(var.cluster)) stop("option cluster for covariance matrix but no var.cluster supplied")
	if(missing(var.cluster)) var.cluster <- NULL

	if(!is.logical(intercept)) stop("non-logical value assigned to option intercept")
	if(!is.logical(table)) stop("non-logical value assigned to option table")
	if(!is.logical(variance)) stop("non-logical value assigned to option variance")
	if(!is.logical(var.eim)) stop("non-logical value assigned to option var.eim")
	if(!is.logical(dfc)) stop("non-logical value assigned to option dfc")

	### 2. Data and variables preparation

	if(any(type==c("1P","2Pfrac"))) x2 <- x

	if(is.data.frame(x)) x <- as.matrix(x)
	if(is.data.frame(x2)) x2 <- as.matrix(x2)

	if(!is.matrix(x)) stop("x is not a matrix")
	if(!is.matrix(x2)) stop("x2 is not a matrix")

	x.names <- dimnames(x)[[2]]
	x2.names <- dimnames(x2)[[2]]

	if(is.null(x.names)) stop("x has no column names")
	if(is.null(x2.names)) stop("x2 has no column names")

	if(intercept==TRUE)
	{
		x <- cbind(1,x)
		x2 <- cbind(1,x2)

		x.names <- c("(Intercept)",x.names)
		x2.names <- c("(Intercept)",x2.names)
	}

	if(length(x.names)!=length(unique(x.names))) stop("some covariate names in x are identical")
	if(length(x2.names)!=length(unique(x2.names))) stop("some covariate names in x2 are identical")

	if(length(y)!=nrow(x)) stop("the number of observations for y and x are different")
	if(!is.null(offset) && length(y)!=length(offset)) stop("offset does not have the appropriate dimension")
	if(var.type=="cluster")
	{
		if(length(y)!=length(var.cluster)) stop("var.cluster does not have the appropriate dimension")
	}

	if(any(type==c("2Pbin","2P")))
	{
		if(inflation==0) yb <- y>0
		if(inflation==1) yb <- y==1
	}

	if(any(type==c("2Pfrac","2P")))
	{
		var.cluf <- NULL
		offsetf <- NULL
		if(inflation==0)
		{
			yf <- y[y>0]
			x2f <- x2[y>0,]
			if(var.type=="cluster") var.cluf <- var.cluster[y>0]
			if(!is.null(offset)) offsetf <- offset[y>0]
		}
		if(inflation==1)
		{
			yf <- y[y<1]
			x2f <- x2[y<1,]
			if(var.type=="cluster") var.cluf <- var.cluster[y<1]
			if(!is.null(offset)) offsetf <- offset[y<1]
		}

		if(length(yf)!=nrow(x2f)) stop("the number of observations for y and x2 are different")
	}

	if(type=="3P")
	{
		var.clu1 <- NULL
		var.cluf <- NULL
		offset1 <- NULL
		offsetf <- NULL
		yb0 <- y>0
		y1_subset <- y[y>0]
		x_subset <- x[y>0,]
		yb1 <- y1_subset==1
		
		yf <- y[y>0 & y<1]
		x2f <- x2[y>0 & y<1,]
		
		if(var.type=="cluster") {
			var.clu1 <- var.cluster[y>0]
			var.cluf <- var.cluster[y>0 & y<1]
		}
		if(!is.null(offset)) {
			offset1 <- offset[y>0]
			offsetf <- offset[y>0 & y<1]
		}
	}

	### 3. Estimation

	class <- "fracreg"

	if(any(type==c("2Pbin","2P")))
	{
		if(var.type=="default") var.ty <- "standard"
		else var.ty <- var.type

		method <- "ML"
		results <- fracreg.est(yb,x,linkbin,method,variance,var.ty,var.eim,var.cluster,dfc,offset=offset,...)
		p <- results$p
		if(variance==TRUE) p.var <- results$p.var
		yhat1 <- results$yhat
		xbhat <- results$xbhat
		converged1 <- results$converged
		LL <- results$LL

		table.info <- list(y=yb,yhat=yhat1,p=p,p.var=if(variance) p.var else NA,x.names=x.names,type="2Pbin",link=linkbin,converged=converged1,var.type=var.ty,LL=LL,method=method,var.cluster=if(var.type=="cluster") var.cluster else NULL,dfc=dfc,or=or,level=level)
		if(table==TRUE) do.call(fracreg.table, table.info)

		formula <- yb ~ x - 1
		names(p) <- x.names

		resBIN <- list(call=cl, class=class,formula=formula,type=type,link=linkbin,method=method,p=p,yhat=yhat1,xbhat=xbhat,converged=converged1,LL=LL,x.names=x.names,table.info=table.info)
		if(variance==TRUE)
		{ 
			dimnames(p.var) <- list(x.names,x.names)
			resBIN[["p.var"]] <- p.var
			resBIN[["var.type"]] <- var.ty
			resBIN[["var.eim"]] <- var.eim
			resBIN[["dfc"]] <- dfc
			if(var.type=="cluster") resBIN[["var.cluster"]] <- var.cluster
		}
		if(!is.null(offset)) resBIN[["offset"]] <- offset
	}

	if(any(type==c("1P","2Pfrac","2P")))
	{
		var.clu <- NULL
		off.curr <- NULL
		if(type=="1P")
		{
			yy <- y
			xx2 <- x2
			ty <- "1P"
			if(var.type=="cluster") var.clu <- var.cluster
			if(!is.null(offset)) off.curr <- offset
		}
		else
		{
			yy <- yf
			xx2 <- x2f
			ty <- "2Pfrac"
			if(var.type=="cluster") var.clu <- var.cluf
			if(!is.null(offset)) off.curr <- offsetf
		}

		if(var.type=="default") var.ty <- "robust"
		else var.ty <- var.type

		method <- "QML"
		results <- fracreg.est(yy,xx2,linkfrac,method,variance,var.ty,var.eim,var.clu,dfc,offset=off.curr,...)
		p <- results$p
		if(variance==TRUE) p.var <- results$p.var
		yhat <- results$yhat
		xbhat <- results$xbhat
		converged2 <- results$converged

		table.info <- list(y=yy,yhat=yhat,p=p,p.var=if(variance) p.var else NA,x.names=x2.names,type=ty,link=linkfrac,converged=converged2,var.type=var.ty,LL=results$LL,method=method,var.cluster=var.clu,dfc=dfc,or=or,level=level)
		if(table==TRUE) do.call(fracreg.table, table.info)

		formula <- yy ~ xx2 - 1
		names(p) <- x2.names

		resFRAC <- list(call=cl, class=class,formula=formula,type=type,link=linkfrac,method=method,p=p,yhat=yhat,xbhat=xbhat,converged=converged2,x.names=x2.names,table.info=table.info)
		if(variance==TRUE)
		{ 
			dimnames(p.var) <- list(x2.names,x2.names)
			resFRAC[["p.var"]] <- p.var
			resFRAC[["var.type"]] <- var.ty
			resFRAC[["var.eim"]] <- var.eim
			resFRAC[["dfc"]] <- dfc
			if(var.type=="cluster") resFRAC[["var.cluster"]] <- var.clu
		}
		if(!is.null(off.curr)) resFRAC[["offset"]] <- off.curr
	}

	if(type=="3P")
	{
		if(var.type=="default") var.ty <- "standard"
		else var.ty <- var.type

		# Component 1: y > 0
		method <- "ML"
		results0 <- fracreg.est(yb0,x,linkbin[1],method,variance,var.ty,var.eim,var.cluster,dfc,offset=offset,...)
		p0 <- results0$p
		if(variance==TRUE) p.var0 <- results0$p.var
		yhat1 <- results0$yhat
		converged1 <- results0$converged
		table.info0 <- list(y=yb0,yhat=yhat1,p=p0,p.var=if(variance) p.var0 else NA,x.names=x.names,type="3Pbin0",link=linkbin[1],converged=converged1,var.type=var.ty,LL=results0$LL,method=method,var.cluster=if(var.type=="cluster") var.cluster else NULL,dfc=dfc,or=or,level=level)
		if(table==TRUE) do.call(fracreg.table, table.info0)
		names(p0) <- x.names
		resBIN0 <- list(call=cl, class=class,formula=yb0 ~ x - 1,type="3Pbin0",link=linkbin[1],method=method,p=p0,yhat=yhat1,xbhat=results0$xbhat,converged=converged1,LL=results0$LL,x.names=x.names,table.info=table.info0)
		if(variance==TRUE) { dimnames(p.var0) <- list(x.names,x.names); resBIN0[["p.var"]] <- p.var0; resBIN0[["var.type"]] <- var.ty }
		if(!is.null(offset)) resBIN0[["offset"]] <- offset

		# Component 2: y = 1 | y > 0
		if(var.type=="cluster") var.clu1 <- var.cluster[y>0] else var.clu1 <- NULL
		results1 <- fracreg.est(yb1,x_subset,linkbin[2],method,variance,var.ty,var.eim,var.clu1,dfc,offset=offset1,...)
		p1 <- results1$p
		if(variance==TRUE) p.var1 <- results1$p.var
		converged2 <- results1$converged
		table.info1 <- list(y=yb1,yhat=results1$yhat,p=p1,p.var=if(variance) p.var1 else NA,x.names=x.names,type="3Pbin1",link=linkbin[2],converged=converged2,var.type=var.ty,LL=results1$LL,method=method,var.cluster=var.clu1,dfc=dfc,or=or,level=level)
		if(table==TRUE) do.call(fracreg.table, table.info1)
		names(p1) <- x.names
		resBIN1 <- list(call=cl, class=class,formula=yb1 ~ x_subset - 1,type="3Pbin1",link=linkbin[2],method=method,p=p1,yhat=results1$yhat,xbhat=results1$xbhat,converged=converged2,LL=results1$LL,x.names=x.names,table.info=table.info1)
		if(variance==TRUE) { dimnames(p.var1) <- list(x.names,x.names); resBIN1[["p.var"]] <- p.var1; resBIN1[["var.type"]] <- var.ty }
		if(!is.null(offset1)) resBIN1[["offset"]] <- offset1

		# Component 3: 0 < y < 1
		if(var.type=="default") var.tyF <- "robust" else var.tyF <- var.type
		resultsF <- fracreg.est(yf,x2f,linkfrac,"QML",variance,var.tyF,var.eim,var.cluf,dfc,offset=offsetf,...)
		pF <- resultsF$p
		if(variance==TRUE) p.varF <- resultsF$p.var
		converged3 <- resultsF$converged
		table.infoF <- list(y=yf,yhat=resultsF$yhat,p=pF,p.var=if(variance) p.varF else NA,x.names=x2.names,type="3Pfrac",link=linkfrac,converged=converged3,var.type=var.tyF,LL=resultsF$LL,method="QML",var.cluster=var.cluf,dfc=dfc,or=or,level=level)
		if(table==TRUE) do.call(fracreg.table, table.infoF)
		names(pF) <- x2.names
		resFRAC <- list(call=cl, class=class,formula=yf ~ x2f - 1,type="3Pfrac",link=linkfrac,method="QML",p=pF,yhat=resultsF$yhat,xbhat=resultsF$xbhat,converged=converged3,x.names=x2.names,table.info=table.infoF)
		if(variance==TRUE) { dimnames(p.varF) <- list(x2.names,x2.names); resFRAC[["p.var"]] <- p.varF; resFRAC[["var.type"]] <- var.tyF }
		if(!is.null(offsetf)) resFRAC[["offset"]] <- offsetf
	}

	if(type=="2P")
	{
		off2Pfrac <- if(!is.null(offset)) offset else 0
		yhat2 <- fracreg.links(linkfrac)$linkinv(x2%*%p + off2Pfrac)
		yhat <- yhat1*yhat2

		converged <- converged1*converged2

		table.info <- list(y=y,yhat=yhat,p=NA,p.var=NA,x.names=NA,type=type,link=c(linkbin,linkfrac),converged=converged,var.type="standard",or=or,level=level)
		if(table==TRUE) do.call(fracreg.table, table.info)

		ybase <- y
		x2base <- x2
	}

	if(type=="3P")
	{
		# Overall Expected Value
		off3Pbin2 <- if(!is.null(offset)) offset else 0
		off3Pfrac <- if(!is.null(offset)) offset else 0

		F1 <- yhat1
		F2 <- fracreg.links(linkbin[2])$linkinv(x%*%p1 + off3Pbin2)
		G <- fracreg.links(linkfrac)$linkinv(x2%*%pF + off3Pfrac)
		yhat <- F1 * (F2 + (1 - F2) * G)

		converged <- converged1*converged2*converged3

		table.info <- list(y=y,yhat=yhat,p=NA,p.var=NA,x.names=NA,type=type,link=c(linkbin,linkfrac),converged=converged,var.type="standard",or=or,level=level)
		if(table==TRUE) do.call(fracreg.table, table.info)

		ybase <- y
		x2base <- x2
	}

	### 4. Return results

	if(type=="2Pbin") {
		class(resBIN) <- "fracreg"
		return(invisible(resBIN))
	}
	if(any(type==c("1P","2Pfrac"))) {
		class(resFRAC) <- "fracreg"
		return(invisible(resFRAC))
	}
	if(type=="2P") {
		res <- list(call=cl, resBIN=resBIN,resFRAC=resFRAC,class=class,type=type,ybase=ybase,x2base=x2base,yhat2P=yhat,converged=converged,table.info=table.info)
		class(res) <- "fracreg"
		return(invisible(res))
	}
	if(type=="3P") {
		res <- list(call=cl, resBIN0=resBIN0,resBIN1=resBIN1,resFRAC=resFRAC,class=class,type=type,ybase=ybase,x2base=x2base,yhat3P=yhat,converged=converged,table.info=table.info)
		class(res) <- "fracreg"
		return(invisible(res))
	}
}

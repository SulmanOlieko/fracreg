fracreg.ptest <- function(object1,object2,version="Wald",table=FALSE)
{
	### 1. Error and warning messages

	if(missing(object1)) stop("object1 is missing")
	if(missing(object2)) stop("object2 is missing")

	if(is.null(object1$class)) stop("object1 not the output of an fracreg command")
	if(is.null(object2$class)) stop("object2 not the output of an fracreg command")

	if(object1$class!="fracreg") stop("object1 is not the output of an fracreg command")
	if(object2$class!="fracreg") stop("object2 is not the output of an fracreg command")

	if(any(object1$type==c("1P","2P","3P")) & all(object2$type!=c("1P","2P","3P"))) stop("object1 and object2 cannot be compared")
	if((object1$type=="2Pbin" | object1$type=="3Pbin0" | object1$type=="3Pbin1") & object1$type!=object2$type) stop("object1 and object2 cannot be compared")
	if(object1$type=="2Pfrac" & object2$type!="2Pfrac") stop("object1 and object2 cannot be compared")

	if(object1$converged==0) stop("object1 is not the output of a successful (converged) fracreg command")
	if(object2$converged==0) stop("object2 is not the output of a successful (converged) fracreg command")

	if(all(version!="LM") & all(version!="Wald")) stop("test version not correctly specified")

	if(!is.logical(table)) stop("non-logical value assigned to option table")

	### 2. Variables and other information for the tests

	# 2.1. Model 1

	if(any(object1$type==c("1P","2Pbin","2Pfrac","3Pbin0","3Pbin1","3Pfrac")))
	{
		link1 <- object1$link
		type1 <- object1$type

		mf <- model.frame(object1$formula)
		y1 <- as.vector(model.response(mf))
		x1 <- model.matrix(object1$formula)
		dimnames(x1)[[2]] <- object1$x.names

		yhat1 <- object1$yhat
		xbhat1 <- object1$xbhat

		g1 <- fracreg.links(link1)$mu.eta(xbhat1)

		gx1 <- g1*x1

		if(type1=="1P") title1 <- paste("Fractional",link1,"regression")
		if(type1=="2Pbin") title1 <- paste("Binary",link1,"component of a two-part regression")
		if(type1=="3Pbin0") title1 <- paste("First binary",link1,"component of a three-part regression")
		if(type1=="3Pbin1") title1 <- paste("Second binary",link1,"component of a three-part regression")
		if(type1=="2Pfrac") title1 <- paste("Fractional",link1,"component of a two-part regression")
		if(type1=="3Pfrac") title1 <- paste("Fractional",link1,"component of a three-part regression")

		if(type1=="2Pbin") type.both <- "2Pbin"
		else type.both <- "others"
	}
		if(object1$type=="3P")
	{
		y1 <- object1$ybase
		link1a <- object1$resBIN0$link; type1a <- "3Pbin0"
		x1a <- model.matrix(object1$resBIN0$formula)
		yhat1a <- object1$resBIN0$yhat
		xbhat1a <- object1$resBIN0$xbhat
		g1a <- fracreg.links(link1a)$mu.eta(xbhat1a)
		dimnames(x1a)[[2]] <- object1$resBIN0$x.names

		link1b <- object1$resBIN1$link; type1b <- "3Pbin1"
		x1b <- model.matrix(object1$resBIN1$formula)
		yhat1b <- object1$resBIN1$yhat
		xbhat1b <- object1$resBIN1$xbhat
		g1b <- fracreg.links(link1b)$mu.eta(xbhat1b)
		dimnames(x1b)[[2]] <- object1$resBIN1$x.names

		link1c <- object1$resFRAC$link; type1c <- "3Pfrac"
		x1c <- object1$x2base
		xbhat1c <- x1c%*%object1$resFRAC$p
		yhat1c <- fracreg.links(link1c)$linkinv(xbhat1c)
		g1c <- fracreg.links(link1c)$mu.eta(xbhat1c)
		dimnames(x1c)[[2]] <- object1$resFRAC$x.names

		yhat1 <- object1$yhat3P
		H <- yhat1b + (1-yhat1b)*yhat1c
		g1a_part <- as.vector(g1a * H)
		g1b_part <- as.vector(yhat1a * g1b * (1 - yhat1c))
		g1c_part <- as.vector(yhat1a * (1 - yhat1b) * g1c)

		gx1 <- cbind(g1a_part*x1a, g1b_part*x1b, g1c_part*x1c)
		title1 <- paste("Three-part regression - binary",link1a,", binary",link1b,"+ fractional",link1c)
		type.both <- "others"
	}

	if(object1$type=="2P")
	{
		y1 <- object1$ybase

		link1a <- object1$resBIN$link
		type1a <- "2Pbin"

		x1a <- model.matrix(object1$resBIN$formula)
		yhat1a <- object1$resBIN$yhat
		xbhat1a <- object1$resBIN$xbhat
		g1a <- fracreg.links(link1a)$mu.eta(xbhat1a)
		dimnames(x1a)[[2]] <- object1$resBIN$x.names

		link1b <- object1$resFRAC$link
		type1b <- "2Pfrac"

		x1b <- object1$x2base
		xbhat1b <- x1b%*%object1$resFRAC$p
		yhat1b <- fracreg.links(link1b)$linkinv(xbhat1b)
		g1b <- fracreg.links(link1b)$mu.eta(xbhat1b)
		dimnames(x1b)[[2]] <- object1$resFRAC$x.names

		yhat1 <- object1$yhat2P
		g1ab <- as.vector(g1a*yhat1b)
		g1ba <- as.vector(g1b*yhat1a)

		gx1 <- cbind(g1ab*x1a,g1ba*x1b)

		title1 <- paste("Binary",link1a,"+ Fractional",link1b,"two-part regression")

		type.both <- "others"
	}

	# 2.2. Model 2

	if(any(object2$type==c("1P","2Pbin","2Pfrac","3Pbin0","3Pbin1","3Pfrac")))
	{
		link2 <- object2$link
		type2 <- object2$type

		mf <- model.frame(object2$formula)
		y2 <- as.vector(model.response(mf))
		x2 <- model.matrix(object2$formula)
		dimnames(x2)[[2]] <- object2$x.names

		yhat2 <- object2$yhat
		xbhat2 <- object2$xbhat
		g2 <- fracreg.links(link2)$mu.eta(xbhat2)

		gx2 <- g2*x2

		if(type2=="1P") title2 <- paste("Fractional",link2,"regression")
		if(type2=="2Pbin") title2 <- paste("Binary",link2,"component of a two-part regression")
		if(type2=="3Pbin0") title2 <- paste("First binary",link2,"component of a three-part regression")
		if(type2=="3Pbin1") title2 <- paste("Second binary",link2,"component of a three-part regression")
		if(type2=="2Pfrac") title2 <- paste("Fractional",link2,"component of a two-part regression")
		if(type2=="3Pfrac") title2 <- paste("Fractional",link2,"component of a three-part regression")
	}
		if(object2$type=="3P")
	{
		y2 <- object2$ybase
		link2a <- object2$resBIN0$link; type2a <- "3Pbin0"
		x2a <- model.matrix(object2$resBIN0$formula)
		yhat2a <- object2$resBIN0$yhat
		xbhat2a <- object2$resBIN0$xbhat
		g2a <- fracreg.links(link2a)$mu.eta(xbhat2a)
		dimnames(x2a)[[2]] <- object2$resBIN0$x.names

		link2b <- object2$resBIN1$link; type2b <- "3Pbin1"
		x2b <- model.matrix(object2$resBIN1$formula)
		yhat2b <- object2$resBIN1$yhat
		xbhat2b <- object2$resBIN1$xbhat
		g2b <- fracreg.links(link2b)$mu.eta(xbhat2b)
		dimnames(x2b)[[2]] <- object2$resBIN1$x.names

		link2c <- object2$resFRAC$link; type2c <- "3Pfrac"
		x2c <- object2$x2base
		xbhat2c <- x2c%*%object2$resFRAC$p
		yhat2c <- fracreg.links(link2c)$linkinv(xbhat2c)
		g2c <- fracreg.links(link2c)$mu.eta(xbhat2c)
		dimnames(x2c)[[2]] <- object2$resFRAC$x.names

		yhat2 <- object2$yhat3P
		H <- yhat2b + (1-yhat2b)*yhat2c
		g2a_part <- as.vector(g2a * H)
		g2b_part <- as.vector(yhat2a * g2b * (1 - yhat2c))
		g2c_part <- as.vector(yhat2a * (1 - yhat2b) * g2c)

		gx2 <- cbind(g2a_part*x2a, g2b_part*x2b, g2c_part*x2c)
		title2 <- paste("Three-part regression - binary",link2a,", binary",link2b,"+ fractional",link2c)
	}

	if(object2$type=="2P")
	{
		y2 <- object2$ybase

		link2a <- object2$resBIN$link
		type2a <- "2Pbin"

		x2a <- model.matrix(object2$resBIN$formula)
		yhat2a <- object2$resBIN$yhat
		xbhat2a <- object2$resBIN$xbhat
		g2a <- fracreg.links(link2a)$mu.eta(xbhat2a)
		dimnames(x2a)[[2]] <- object2$resBIN$x.names

		link2b <- object2$resFRAC$link
		type2b <- "2Pfrac"

		x2b <- object2$x2base
		xbhat2b <- x2b%*%object2$resFRAC$p
		yhat2b <- fracreg.links(link2b)$linkinv(xbhat2b)
		g2b <- fracreg.links(link2b)$mu.eta(xbhat2b)
		dimnames(x2b)[[2]] <- object2$resFRAC$x.names

		yhat2 <- object2$yhat2P
		g2ab <- as.vector(g2a*yhat2b)
		g2ba <- as.vector(g2b*yhat2a)

		gx2 <- cbind(g2ab*x2a,g2ba*x2b)

		title2 <- paste("Binary",link2a,"+ Fractional",link2b,"two-part regression")
	}

	### 3. Further error and warning messages

	if(!isTRUE(all.equal(y1,y2))) stop("The dependent variable is not the same in the two models")

	if(any(object1$type==c("2Pbin","2Pfrac")) | (object1$type=="1P" & object2$type=="1P"))
	{
		if(object1$link==object2$link)
		{
			x1.names <- dimnames(x1)[[2]]
			x2.names <- dimnames(x2)[[2]]
			x12.names <- c(x1.names,x2.names)
			x12.names <- unique(x12.names)
			x1.len <- length(x1.names)
			x2.len <- length(x2.names)
			x12.len <- length(x12.names)

			if(identical(x1.names,x2.names)) stop("object 1 and object 2 are based on the same link function and covariates")
			if(x1.len==x12.len & x2.len!=x12.len) stop("object 2 is nested in object 1 - no need to use the P test")
			if(x1.len!=x12.len & x2.len==x12.len) stop("object 1 is nested in object 2 - no need to use the P test")
		}
	}
	if(object1$type=="2P" & object2$type=="2P")
	{
		if(object1$resBIN$link==object2$resBIN$link & object1$resFRAC$link==object2$resFRAC$link)
		{
			x1a.names <- dimnames(x1a)[[2]]
			x2a.names <- dimnames(x2a)[[2]]
			x12a.names <- c(x1a.names,x2a.names)
			x12a.names <- unique(x12a.names)
			x1a.len <- length(x1a.names)
			x2a.len <- length(x2a.names)
			x12a.len <- length(x12a.names)

			x1b.names <- dimnames(x1b)[[2]]
			x2b.names <- dimnames(x2b)[[2]]
			x12b.names <- c(x1b.names,x2b.names)
			x12b.names <- unique(x12b.names)
			x1b.len <- length(x1b.names)
			x2b.len <- length(x2b.names)
			x12b.len <- length(x12b.names)

			if(identical(x1a.names,x2a.names) & identical(x1b.names,x2b.names)) stop("object 1 and object 2 are based on the same link function and covariates")
			if(x1a.len==x12a.len & x2a.len==x12a.len & x1b.len==x12b.len & x2b.len!=x12b.len) stop("object 2 is nested in object 1 - no need to use the P test")
			if(x1a.len==x12a.len & x2a.len==x12a.len & x1b.len!=x12b.len & x2b.len==x12b.len) stop("object 1 is nested in object 2 - no need to use the P test")
			if(x1a.len==x12a.len & x2a.len!=x12a.len & x1b.len==x12b.len & x2b.len==x12b.len) stop("object 2 is nested in object 1 - no need to use the P test")
			if(x1a.len!=x12a.len & x2a.len==x12a.len & x1b.len==x12b.len & x2b.len==x12b.len) stop("object 1 is nested in object 2 - no need to use the P test")
			if(x1a.len==x12a.len & x2a.len!=x12a.len & x1b.len==x12b.len & x2b.len!=x12b.len) stop("object 2 is nested in object 1 - no need to use the P test")
			if(x1a.len!=x12a.len & x2a.len==x12a.len & x1b.len!=x12b.len & x2b.len==x12b.len) stop("object 1 is nested in object 2 - no need to use the P test")
		}
	}

	### 4. Tests

	ver <- NA
	S <- NA
	Sp <- NA

	df <- 1
	for(j in 1:2)
	{
		if(j==1)
		{
			yj <- y1
			yhatj <- yhat1
			gxj <- gx1
			gzj <- yhat2-yhat1
		}
		if(j==2)
		{
			yj <- y2
			yhatj <- yhat2
			gxj <- gx2
			gzj <- yhat1-yhat2
		}

		if(any(version=="LM"))
		{
			ver <- c(ver,"LM")

			results <- fracreg.lm(yj,yhatj,gxj,gzj,type.both)
			Sj <- results$LM
			S <- c(S,Sj)
			Sp <- c(Sp,1-pchisq(Sj,df))
		}
		if(any(version=="Wald"))
		{
			ver <- c(ver,"Wald")

			yt <- yj-yhatj
			gxzj <- cbind(gxj,gzj)
			results <- lm(yt ~ gxzj-1)

			gzj.b <- results$coefficients[ncol(gxzj)]
			dfcc <- nrow(gxzj)/(nrow(gxzj)-ncol(gxzj))
			gzj.var <- dfcc*(solve(t(gxzj)%*%gxzj)%*%t(gxzj)%*%diag(results$residuals^2)%*%gxzj%*%solve(t(gxzj)%*%gxzj))[ncol(gxj)+1,ncol(gxj)+1]

			Sj <- as.vector(gzj.b/(gzj.var^0.5))
			S <- c(S,Sj)
			Spj <- 2*(1-pt(abs(Sj),nrow(gxzj)-ncol(gxzj)))
			Sp <- c(Sp,Spj)
		}
	}

	n.ver <- length(ver[-1])/2 
	table.info <- list(test.which="P",S=S,Sp=Sp,ver=ver,title1=title1,title2=title2,n.ver=n.ver)
	if(table==TRUE) do.call(fracreg.tests.table, table.info)

	### 5. Return results

	statistics=S[-1]
	ver <- ver[-1]
	names(statistics) <- paste(c(rep("H0-obj1",n.ver),rep("H0-obj2",n.ver)),ver,sep="-")

	class(statistics) <- c("fracreg.ptest", "numeric")
	attr(statistics, "table.info") <- table.info

	return(invisible(statistics))
}

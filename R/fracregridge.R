#' Fractional Ridge Regression
#'
#' @description 
#' \code{fracregridge} implements Fractional Ridge Regression (Rokem & Kay, 2020), which is an approach to regularized linear regression.
#' Unlike standard ridge regression where the penalty term \eqn{\alpha} is chosen directly, \code{fracregridge} allows you to specify the desired *fraction* of the unregularized OLS coefficient vector length. The algorithm then automatically determines the corresponding \eqn{\alpha} penalties.
#'
#' @param y A numeric vector or matrix of dependent variables.
#' @param x A numeric matrix of explanatory variables.
#' @param fracs A numeric vector of desired fractions of the unregularized coefficient vector length. Default is \code{seq(0.1, 1.0, by=0.1)}. Must be sorted in ascending order.
#' @param tol A numeric tolerance under which singular values of the \code{x} matrix are considered to be zero. Default is \code{1e-10}.
#' @param intercept logical. If \code{TRUE}, an intercept is included in the model. Default is \code{TRUE}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class \code{fracregridge}, containing:
#' \item{coef}{The estimated ridge coefficients for each requested fraction.}
#' \item{alphas}{The corresponding \eqn{\alpha} penalty values.}
#' \item{fracs}{The grid of fractions.}
#' \item{call}{The matched call.}
#' @references Rokem, A., and Kay, K. (2020), Fractional ridge regression: a fast, 
#'   interpretable reparameterization of ridge regression. GigaScience, 9(12): giaa133.
#' @references Rokem, A., and Kay, K., fracridge: Fractional Ridge Regression. Package repository. <https://github.com/nrdg/fracridge>.
#' @export
#'
#' @examples
#' set.seed(123)
#' y <- rnorm(100)
#' x <- matrix(rnorm(1000), 100, 10)
#' colnames(x) <- paste0("x", 1:10)
#' 
#' # Fit fractional ridge regression
#' mod <- fracregridge(y, x, fracs = c(0.3, 0.5, 0.8))
#' print(mod)
#' summary(mod)
fracregridge <- function(y, x, fracs = seq(0.1, 1.0, by=0.1), tol = 1e-10, intercept = TRUE, ...)
{
	cl <- match.call()

	if(missing(y)) stop("dependent variable is missing")
	if(missing(x)) stop("explanatory variables are missing")

	if(is.data.frame(x)) x <- as.matrix(x)
	if(!is.matrix(x)) stop("x is not a matrix")

	x.names <- dimnames(x)[[2]]
	if(is.null(x.names)) stop("x has no column names")

	if(intercept==TRUE)
	{
		x <- cbind(1,x)
		x.names <- c("(Intercept)",x.names)
	}

	if(length(x.names)!=length(unique(x.names))) stop("some covariate names in x are identical")
	if(length(y)!=nrow(x)) stop("the number of observations for y and x are different")
	
	if(any(diff(fracs) < 0)) stop("The 'fracs' input must be sorted in ascending order.")
	
	nn <- nrow(x)
	pp <- ncol(x)
	
	y_mat <- as.matrix(y)
	bb <- ncol(y_mat)
	ff <- length(fracs)
	
	# SVD Rotation
	if (nn > pp) {
		s <- svd(crossprod(x))
		uu <- s$u
		ss <- s$d
		v_t <- t(s$v)
		selt <- sqrt(pmax(ss, 0))
		
		if (bb >= nn) {
			ynew <- (diag(1/selt) %*% v_t %*% t(x)) %*% y_mat
		} else {
			ynew <- diag(1/selt) %*% v_t %*% (t(x) %*% y_mat)
		}
	} else {
		s <- svd(x)
		uu <- s$u
		selt <- s$d
		v_t <- t(s$v)
		ynew <- t(uu) %*% y_mat
	}
	
	ols_coef <- t(t(ynew) / selt)
	
	isbad <- selt < tol
	if (any(isbad)) {
		warning("Some eigenvalues are being treated as 0")
		ols_coef[isbad, ] <- 0
	}
	
	BIG_BIAS <- 10e3
	SMALL_BIAS <- 10e-3
	BIAS_STEP <- 0.2
	
	val1 <- BIG_BIAS * (selt[1]^2)
	val2 <- SMALL_BIAS * (selt[length(selt)]^2)
	
	alphagrid <- c(0, 10^seq(floor(log10(val2)), ceiling(log10(val1)), by = BIAS_STEP))
	
	seltsq <- selt^2
	
	sclg <- outer(seltsq, alphagrid, function(s, a) s / (s + a))
	sclg_sq <- sclg^2
	
	first_dim <- min(nn, pp)
	
	coef_res <- array(NA, dim = c(first_dim, ff, bb))
	alphas_res <- matrix(NA, nrow = ff, ncol = bb)
	
	for (ii in seq_len(bb)) {
		newlen <- sqrt(as.numeric(t(sclg_sq) %*% (ols_coef[, ii]^2)))
		newlen <- newlen / newlen[1]
		
		# interpolate (log transformed space)
		temp <- stats::approx(x = rev(newlen), y = rev(log(1 + alphagrid)), xout = fracs, rule = 2)$y
		
		targetalphas <- exp(temp) - 1
		alphas_res[, ii] <- targetalphas
		
		sc <- outer(seltsq, targetalphas, function(s, a) s / (s + a))
		
		c_mat <- sc * ols_coef[, ii]
		coef_res[,, ii] <- c_mat
	}
	
	coef_flat <- matrix(coef_res, nrow = first_dim, ncol = ff * bb)
	coef_unrot <- t(v_t) %*% coef_flat
	
	coef_final <- array(coef_unrot, dim = c(pp, ff, bb))
	
	if (bb == 1) {
		coef_final <- matrix(coef_final[,,1], nrow = pp, ncol = ff)
		colnames(coef_final) <- paste0("frac_", fracs)
		rownames(coef_final) <- x.names
		alphas_res <- as.numeric(alphas_res[, 1])
		names(alphas_res) <- paste0("frac_", fracs)
	} else {
		dimnames(coef_final) <- list(x.names, paste0("frac_", fracs), colnames(y_mat))
	}
	
	yhat <- array(NA, dim = c(nn, ff, bb))
	xbhat <- array(NA, dim = c(nn, ff, bb))
	se <- array(NA, dim = c(pp, ff, bb))
	table.info <- list()
	stats.info <- list()
	
	V2 <- t(v_t)^2
	warned_np <- FALSE
	
	for (ii in seq_len(bb)) {
	    for (jj in seq_len(ff)) {
	        if (bb == 1) {
	            beta <- coef_final[, jj]
	        } else {
	            beta <- coef_final[, jj, ii]
	        }
	    
	        y_pred <- x %*% beta
	        yhat[, jj, ii] <- y_pred
	        xbhat[, jj, ii] <- y_pred
	        
	        if (bb == 1) {
	            alpha <- alphas_res[jj]
	        } else {
	            alpha <- alphas_res[jj, ii]
	        }
	        
	        df_alpha <- sum(seltsq / (seltsq + alpha))
	        
	        if (nn > df_alpha) {
	            sigma2 <- sum((y_mat[, ii] - y_pred)^2) / (nn - df_alpha)
	            D <- seltsq / (seltsq + alpha)^2
	            var_diag <- V2 %*% D
	            std_err <- sqrt(as.numeric(sigma2 * var_diag))
	        } else {
	            if (!warned_np) {
	                warning("Degrees of freedom (n - df(alpha)) is less than or equal to 0. Standard errors cannot be computed and will be NA.")
	                warned_np <- TRUE
	            }
	            std_err <- rep(NA, pp)
	        }
	        
	        se[, jj, ii] <- std_err
	        
	        est <- beta
	        zval <- est / std_err
	        pval <- 2 * pnorm(-abs(zval))
	        
	        tab <- cbind(Estimate = est, `Std. Error` = std_err, `z value` = zval, `Pr(>|z|)` = pval)
	        rownames(tab) <- x.names
	        
	        list_name <- paste0("target_", colnames(y_mat)[ii], "_frac_", fracs[jj])
	        if (bb == 1) list_name <- paste0("frac_", fracs[jj])
	        
	        table.info[[list_name]] <- tab
	        stats.info[[list_name]] <- list(n_obs = nn, R2 = cor(y_mat[, ii], y_pred)^2, df_alpha = df_alpha)
	    }
	}
	
	p.var <- se^2
	if (bb == 1) {
	    p.var <- matrix(p.var[,,1], nrow = pp, ncol = ff)
	    colnames(p.var) <- paste0("frac_", fracs)
	    rownames(p.var) <- x.names
	    
	    yhat_out <- matrix(yhat[,,1], nrow = nn, ncol = ff)
	    colnames(yhat_out) <- paste0("frac_", fracs)
	    
	    xbhat_out <- matrix(xbhat[,,1], nrow = nn, ncol = ff)
	    colnames(xbhat_out) <- paste0("frac_", fracs)
	} else {
	    dimnames(p.var) <- list(x.names, paste0("frac_", fracs), colnames(y_mat))
	    dimnames(yhat) <- list(NULL, paste0("frac_", fracs), colnames(y_mat))
	    dimnames(xbhat) <- list(NULL, paste0("frac_", fracs), colnames(y_mat))
	    yhat_out <- yhat
	    xbhat_out <- xbhat
	}
	
	res <- list(call = cl, 
	            coef = coef_final, 
	            alphas = alphas_res, 
	            fracs = fracs, 
	            x.names = x.names,
	            type = "fracregridge",
	            class = "fracregridge",
	            yhat = yhat_out,
	            xbhat = xbhat_out,
	            p.var = p.var,
	            table.info = table.info,
	            stats.info = stats.info)
	class(res) <- "fracregridge"
	
	return(res)
}

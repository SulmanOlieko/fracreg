#' @title Partial Effects for Fractional Panel Data Regression
#'
#' @description
#' Computes Average Partial Effects (APEs) for fractional panel data models fit using \code{fracregpd}. In correlated random effects (CRE) models, the time-averages of the covariates are used merely to control for unobserved heterogeneity, and as such, they are automatically filtered out when computing the structural partial effects. Standard errors are computed using the Delta method and the bootstrapped parameter variance matrix.
#' @param object An object of class \code{fracregpd}.
#' @param APE logical. Compute Average Partial Effects?
#' @param CPE logical. Compute Conditional Partial Effects? (Not currently supported for panel data models).
#' @param at numeric vector. The values at which to evaluate the CPE.
#' @param which.x character vector. Variables for which to compute partial effects. By default, auxiliary CRE parameters (like \code{_mean} and \code{vhat}) are automatically excluded so that partial effects are only computed for the main structural parameters.
#' @param variance logical. Compute standard errors using the Delta method?
#' @param table logical. Print the resulting partial effects table?
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' An object of class \code{fracreg.pe} containing the standard coefficient tables with Average Partial Effects.
#'
#' @author Sulman Olieko Owili <oliekosulman@gmail.com>
#'
#' @seealso
#' \code{\link{fracregpd}}, \code{\link{fracreg.pe}}
#'
#' @examples
#' \donttest{
#' # Simulate Panel Data
#' N <- 50
#' T <- 5
#' id <- rep(1:N, each=T)
#' time <- rep(1:T, N)
#' x1 <- rnorm(N*T)
#' x2 <- rnorm(N*T)
#' z1 <- rnorm(N*T)
#' u <- rnorm(N*T)
#' y <- exp(x1 - x2 + u)/(1 + exp(x1 - x2 + u))
#' X <- cbind(x1 = x1, x2 = x2)
#' Z <- cbind(x1 = x1, z1 = z1)
#' 
#' # Estimate panel data model
#' mod_pd <- fracregpd(id=id, time=time, y=y, x=X, z=Z, 
#'                     x.exogenous=FALSE, type="GMMbgw", link="logit")
#' 
#' # Compute Average Partial Effects
#' pe_res <- fracregpd.pe(mod_pd)
#' summary(pe_res)
#' }
#' @export
fracregpd.pe <- function(object, APE=TRUE, CPE=FALSE, at=NULL, which.x=NULL, variance=TRUE, table=FALSE, ...) {
    if(missing(object)) stop("object is missing")
    if(is.null(object$class)) {
        if(!inherits(object, "fracregpd")) stop("object is not the output of a fracregpd command")
    } else {
        if(object$class != "fracregpd") stop("object is not the output of a fracregpd command")
    }

    if(!is.logical(APE)) stop("non-logical value assigned to option APE")
    if(!is.logical(CPE)) stop("non-logical value assigned to option CPE")
    if(!is.logical(variance)) stop("non-logical value assigned to option variance")
    if(!is.logical(table)) stop("non-logical value assigned to option table")

    if(all(c(APE,CPE)==FALSE)) stop("You must specify at least one option: APE and/or CPE")
    if(CPE==FALSE & !is.null(at)) stop("option at is only required for CPE")

    if(object$converged==0) stop("object is not the output of a successful (converged) fracregpd command")

    x.names <- object$table.info$x.names
    
    # Filter out intercept and CRE parameters
    aux_vars <- grep("_mean$|^vhat$", x.names, value = TRUE)
    struct_vars <- setdiff(x.names, aux_vars)
    if ("(Intercept)" %in% struct_vars) {
        struct_vars <- setdiff(struct_vars, "(Intercept)")
    }
    if (length(struct_vars) == 0) stop("No structural variables found to compute partial effects.")

    if(is.null(which.x)) which.x <- struct_vars
    xw.names <- unique(c(struct_vars,which.x))
    if(!all(which.x %in% x.names)) stop("option which.x contains variables not in the model")

    xvar.names <- struct_vars
    k <- length(xvar.names)
    npar <- length(object$p)
    n <- nrow(object$x)

    linka <- object$link
    pa <- object$p
    pa.var <- object$table.info$p.var
    x <- object$x
    xbhata <- object$xbhat
    
    title <- paste("Panel data fractional", linka, "regression")
    type <- "1P" # Treat structurally as 1P for PE variance computation

    resAPE <- list()
    resCPE <- list()

    if(APE==TRUE) {
        PE.type <- "APE"

        if(any(x.names=="(Intercept)")) p.pe <- matrix(rep(pa[-1],each=n),ncol=npar-1)
        else p.pe <- matrix(rep(pa,each=n),ncol=npar)
        
        # We need to map p.pe columns to x.names (excluding intercept)
        xnames_noint <- setdiff(x.names, "(Intercept)")
        dimnames(p.pe) <- list(NULL, xnames_noint)

        ga <- fracregpd.links(linka)$mu.eta(xbhata)
        PEa.p <- as.matrix(p.pe[,which.x])*ga
        PE.p <- apply(PEa.p,2,mean)

        if(variance==TRUE) {
            # Note: fracreg.pe.var relies on type, p, etc. 
            # We map the inputs directly to fracreg.pe.var expectations
            PE.sd <- fracreg.pe.var(x,npar,which.x,x.names,xnames_noint,"1P",pa,xbhata,ga,linka,pa.var)
        } else {
            PE.sd <- NA
        }

        table.info.APE <- list(PE.p=PE.p,PE.sd=if(variance) PE.sd else NA,PE.type=PE.type,which.x=which.x,xvar.names=xvar.names,title=title,at=at)
        if(table==TRUE) do.call(fracreg.pe.table, table.info.APE)
        resAPE[["table.info"]] <- table.info.APE
    }

    if(CPE==TRUE) {
        PE.type <- "CPE"
        # CPE computation similar to fracreg.pe
        # Since it's panel data with CRE, the user needs to provide values for ALL parameters including _mean.
        # This is complex, so we will throw an error if CPE is requested for now unless we fully support it.
        stop("CPE is currently not supported for fracregpd models due to CRE auxiliary parameters. Use APE.")
    }

    if(APE==TRUE & CPE==TRUE) res <- list(ape=resAPE,cpe=resCPE)
    else if(APE==TRUE & CPE==FALSE) res <- resAPE
    else if(APE==FALSE & CPE==TRUE) res <- resCPE
    
    class(res) <- c("fracregpd.pe", "fracreg.pe")
    return(invisible(res))
}

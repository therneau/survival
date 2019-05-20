# $Id: ridge.S 11166 2008-11-24 22:10:34Z therneau $
ridge <- function(..., theta, df=nvar/2, eps=.1, scale=TRUE) {
    x <- cbind(...)
    nvar <- ncol(x)
    xname <- as.character(parse(text=substitute(cbind(...))))[-1]
    vars <- apply(x, 2, function(z) var(z[!is.na(z)]))
    class(x) <- 'coxph.penalty'

    if (!missing(theta) && !missing(df))
	    stop("Only one of df or theta can be specified")

    if (scale) 
	    pfun <- function(coef,theta, ndead, scale) {
		list(penalty= sum(coef^2 *scale)*theta/2,
		     first  = theta*coef*scale,
		     second = theta*scale,
		     flag=FALSE)
		}
    else
	    pfun <- function(coef,theta, ndead, scale) {
		list(penalty= sum(coef^2)*theta/2,
		     first  = theta*coef,
		     second = theta,
		     flag=FALSE)
		}


    if (!missing(theta)) {
	temp <- list(pfun=pfun,
		     diag=TRUE,
		     cfun=function(parms, iter, history) {
				list(theta=parms$theta, done=TRUE) }, 
		     cparm=list(theta= theta),
		     pparm= vars,
		     varname=paste('ridge(', xname, ')', sep=''))
	}
    else {
	temp <- list(pfun=pfun,
		     diag=TRUE,
		     cfun=frailty.controldf,
		     cargs = 'df',
		     cparm=list(df=df, eps=eps, thetas=0, dfs=nvar,
		         guess=1),
		     pparm= vars,
		     varname=paste('ridge(', xname, ')', sep=''))
	}
	
    attributes(x) <- c(attributes(x), temp)
    x
    }

#  $Id: frailty.gamma.S 11377 2009-12-14 22:59:56Z therneau $
# 
# Defining function for gamma frailty fits
#
frailty.gamma <- function(x, sparse=(nclass >5), theta, df, eps= 1e-5, 
			  method=c("em", "aic", "df", "fixed"), ...) {
    nclass <- length(unique(x[!is.na(x)]))
    if (sparse)	x <-as.numeric(as.factor(x))
    else{
	x <- as.factor(x)
	attr(x,'contrasts') <- contr.treatment(nclass, contrasts=FALSE)
                }
    if (is.R()) class(x) <- c("coxph.penalty",class(x))
    else        oldClass(x) <- "coxph.penalty"

    # Check for consistency of the arguments
    if (missing(method)) {
	if (!missing(theta)) {
	    method <- 'fixed'
	    if (!missing(df)) 
		    stop("Cannot give both a df and theta argument")
	    }
	else if (!missing(df)) method <- 'df'
	}
    method <- match.arg(method)
    if (method=='df' && missing(df)) stop("Method = df but no df argument")
    if (method=='fixed' && missing(theta))
	    stop("Method= fixed but no theta argument")
    if (method!='df' && !missing(df)) 
	    stop("Method is not df, but have a df argument")
    if (method !='fixed' && !missing(theta)) 
	    stop("Method is not 'fixed', but have a theta argument")

    pfun<- function(coef, theta, ndeath){
	if (theta==0) list(recenter=0, penalty=0, flag=TRUE)
	else {
	      recenter <- log(mean(exp(coef)))
	      coef <- coef - recenter
	      nu <- 1/theta
	      list(recenter=recenter,
		   first=   (exp(coef) -1) * nu,
		   second=  exp(coef) * nu,
		   penalty= -sum(coef)*nu,   # The exp part sums to a constant
		   flag=FALSE)
	           }
	   }

    printfun <- function(coef, var, var2, df, history) {
	if (!is.null(history$history)) 
	     theta <- history$history[nrow(history$history),1]
	else theta <- history$theta
	clog  <- history$c.loglik
	
	if (is.matrix(var)) test <- coxph.wtest(var, coef)$test
	else 		    test <- sum(coef^2/var)
	df2 <- max(df, .5)      # Stop silly p-values
	list(coef=c(NA, NA, NA, test, df, 1-pchisq(test, df2)),
		 history=paste("Variance of random effect=", format(theta),
	                       "  I-likelihood =", 
		         format(round(clog,1), digits=10)))
	}
    # The final coxph object will contain a copy of printfun.  Stop it from
    #   also containing huge unnecessary variables, e.g. 'x', known at this 
    #   point in time.  Not an issue for pfun, which does not get saved.
    # Setting to globalenv() won't work since coxph.wtest is not visible 
    #   outside the survival library's name space.
    if (is.R()) environment(printfun) <- asNamespace('survival')

    if (method=='fixed') {
	temp <- list(pfun=pfun,
		     printfun=printfun,
		     diag =TRUE,
		     sparse= sparse,
		     cargs = c("x", "status", "loglik"),		 
		     cfun = frailty.controlgam,
		     cparm= list(theta=theta, ...))
        }
    else if (method=='em'){
	temp <- list(pfun=pfun,
		     printfun=printfun,
		     diag =TRUE,
		     sparse= sparse,
		     cargs = c("x", "status", "loglik"),		 
		     cfun = frailty.controlgam,
		     cparm= c(list(eps=eps), ...))
	}
    
    else if (method=='aic') {
	temp <- list(pfun=pfun,
		     printfun=printfun,
		     diag =TRUE,
		     sparse= sparse,
		     cargs = c("x", "status", "loglik", "neff","df", "plik"),
		     cparm=list(eps=eps, lower=0, init=c(.1, 1), ...),
		     cfun =function(opt, iter, old, group, status, loglik,...){
			 temp <- frailty.controlaic(opt, iter, old, ...)
			 if (iter >0) {
			     #compute correction to the loglik
			     if (old$theta==0) correct <- 0
			     else {
				 if (is.matrix(group)) 
					 group <-c(group %*% 1:ncol(group))
				 d <- tapply(status,group,sum)
				 correct <- frailty.gammacon(d, 1/old$theta)
				 }
			     temp$c.loglik <- loglik + correct
			     }
			 temp
			 })
	}
    else {  #df method
	# The initial guess is based on the observation that theta=1 often
	#   gives about df= (#groups)/3
	if (missing(eps)) eps <- .1
	temp <- list(pfun=pfun,
		     printfun=printfun,
		     diag =TRUE,
		     sparse= sparse,
		     cargs= c('df', "x", "status", "loglik"),
		     cparm=list(df=df, thetas=0, dfs=0, eps=eps,
		                guess=3*df/length(unclass(x)), ...),
		     cfun =function(opt, iter, old, df, group, status, loglik){
			 temp <- frailty.controldf(opt, iter, old, df)
			 if (iter >0) {
			     #compute correction to the loglik
			     if (old$theta==0) correct <- 0
			     else {
				 if (is.matrix(group)) 
					 group <-c(group %*% 1:ncol(group))
				 d <- tapply(status,group,sum)
				 correct <- frailty.gammacon(d, 1/old$theta)
				 }
			     temp$c.loglik <- loglik + correct
			     }

			 temp
		         })
	}

    # If not sparse, give shorter names to the coefficients, so that any
    #   printout of them is readable.
    if (!sparse) {
	vname <- paste("gamma", levels(x), sep=':')
	temp <- c(temp, list(varname=vname))
	}
    attributes(x) <- c(attributes(x), temp)
    x
    }

			  
			   
			   

# $Id: frailty.t.S 11377 2009-12-14 22:59:56Z therneau $
# 
# Defining function for t-distribution frailty fits
#
frailty.t <- function(x, sparse=(nclass>5), theta, df, eps= 1e-5,  tdf=5,
			  method=c("aic", "df", "fixed"), ...) {
    nclass <- length(unique(x[!is.na(x)]))
    if (sparse){
	x <-as.numeric(as.factor(x))
	if (is.R()) class(x) <- "coxph.penalty"
	else        oldClass(x) <- "coxph.penalty"
        }
    else{
	x <- as.factor(x)
	if (is.R()) class(x) <- c("coxph.penalty",class(x))
	else        oldClass(x) <- "coxph.penalty"
	attr(x,'contrasts') <- contr.treatment(nclass, contrasts=FALSE)
        }

    if (tdf <=2) stop("Cannot have df <3 for the t-frailty")
    # Check for consistency of the arguments
    if (missing(method)) {
	if (!missing(theta)) {
	    method <- 'fixed'
	    if (!missing(df)) 
		    stop("Cannot give both a df and theta argument")

	    }
	else if (!missing(df)) {
	    if (df==0) method <- 'aic'
	    else       method <- 'df'
	    }
	}
    method <- match.arg(method)
    if (method=='df' && missing(df)) stop("Method = df but no df argument")
    if (method=='fixed' && missing(theta))
	    stop("Method= fixed but no theta argument")
    if (method !='fixed' && !missing(theta)) 
	    stop("Method is not 'fixed', but have a theta argument")

    pfun<- function(coef, theta, ndead, tdf){
	if (theta==0) list(recenter=0, penalty=0, flag=TRUE)
	else {
	    sig <- theta* (tdf-2)/tdf  #scale contant^2 in density formula
	    #
	    # Find the centering constant, using 1 NR step
	    #
	    temp  <- 1 + coef^2/(tdf*sig)
	    temp1 <- coef/temp
	    temp2 <- 1/temp - (2/(tdf*sig))*coef^2/temp^2
	    recenter <- sum(temp1)/sum(temp2)  #NR step towards MLE

	    coef <- coef - recenter
	    const <- (tdf+1)/(tdf*sig)
	    temp  <- 1 + coef^2/(tdf*sig)
	    list(recenter=recenter, 
		 first=   const*coef/temp,
		 second=  const*(1/temp - (2/(tdf*sig))*coef^2/temp^2),
		 penalty= sum(.5*log(pi*tdf*sig) + ((tdf+1)/2)*log(temp) +
		                lgamma(tdf/2) - lgamma((tdf+1)/2)),
		 flag=FALSE)
	    }
	}

    printfun <- function(coef, var, var2, df, history) {
	if (!is.null(history$history)) 
	     theta <- history$history[nrow(history$history),1]
	else theta <- history$theta
	
	if (is.matrix(var)) test <- coxph.wtest(var, coef)$test
	else 		    test <- sum(coef^2/var)
	df2 <- max(df, .5)      # Stop silly p-values
	list(coef=c(NA, NA, NA, test, df, 1-pchisq(test, df2)),
		 history=paste("Variance of random effect=", format(theta)))
	}
    # The final coxph object will contain a copy of printfun.  Stop it from
    #   also containing huge unnecessary variables, e.g. 'x', known at this 
    #   point in time.  Not an issue for pfun, which does not get saved.
    # The reason for using the survival namespace instead of globalenv() is 
    # that we call coxph.wtest, which may not be visible outside the name space
    if (is.R()) environment(printfun) <- asNamespace('survival')

    if (method=='fixed') {
	temp <- list(pfun=pfun, pparm=tdf,
		     printfun=printfun,
		     diag =TRUE,
		     sparse= sparse,
		     cfun = function(parms, iter, old){
		          list(theta=parms$theta, done=TRUE)},
		     cparm= list(theta=theta, ...))
        }
    
    else if (method=='aic') {
	temp <- list(pfun=pfun, pparm=tdf,
		     printfun=printfun,
		     diag =TRUE,
		     sparse= sparse,
		     cargs = c("neff", "df", "plik"),	
		     cparm=list(lower=0, init=c(.1,1), eps=eps, ...),
		     cfun = frailty.controlaic)
	}
    else {  #df method
	if (missing(eps)) eps <- .1
	temp <- list(pfun=pfun, pparm=tdf,
		     printfun=printfun,
		     diag =TRUE,
		     sparse= sparse,
		     cargs= c('df'),
		     cparm=list(df=df, eps=eps, thetas=0, dfs=0,
		                guess=3*df/length(unclass(x)), ...),
                     cfun = frailty.controldf)
	}

    # If not sparse, give shorter names to the coefficients, so that any
    #   printout of them is readable.
    if (!sparse) {
	vname <- paste("t", levels(x), sep=':')
	temp <- c(temp, list(varname=vname))
	}
    attributes(x) <- c(attributes(x), temp)
    x
    }

			  
			   
			   

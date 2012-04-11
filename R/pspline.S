#
# the p-spline function for a Cox model
#
pspline <- function(x, df=4, theta, nterm=2.5*df, degree=3, eps=0.1, 
		    method, Boundary.knots=range(x), 
                    intercept=FALSE, penalty=TRUE, ...) {
    if (!missing(theta)) {
	method <- 'fixed'
	if (theta <=0 || theta >=1) stop("Invalid value for theta")
	}
    else if (df ==0 || (!missing(method) && method=='aic')) {
	method <- 'aic'
	nterm <- 15    #will be ok for up to 6-8 df
	if (missing(eps)) eps <- 1e-5
	}
    else {
	method <- 'df'
	if (df <=1) stop ('Too few degrees of freedom')
	# The below used to say "df+1 > nterm", but we need some scope for
	#  the smoother parameter to avoid strange conditions
        if (df > nterm) stop("`nterm' too small for df=",df)
    }

    xname <- deparse(substitute(x))
    keepx <- !is.na(x)
    if (!all(keepx)) x <- x[keepx] #this is done before any reference to 
                                   # Boundary.knots, so the default works
    nterm <- round(nterm)
    if (nterm < 3) stop("Too few basis functions")
    
    if (!missing(Boundary.knots)) {
        if (!is.numeric(Boundary.knots) || length(Boundary.knots) !=2 ||
            Boundary.knots[1] >= Boundary.knots[2])
            stop("Invalid values for Boundary.knots")
            
        # Check for data values outside the knot range
        outl <- (x < Boundary.knots[1])
        outr<- (x > Boundary.knots[2])
        outside <- outl | outr
    }
    else outside <- FALSE

    # Set up the evenly spaced knots
    dx <- (Boundary.knots[2] - Boundary.knots[1])/nterm
    knots <- c(Boundary.knots[1] + dx*((-degree):(nterm-1)), 
               Boundary.knots[2]+ dx*(0:degree))

    # Set up the basis.  Inside the boundary knots we use spline.des.
    # Outside of them we use  f(edge) + (x-edge)* f'(edge)
    if (any(outside)) {
        newx <- matrix(0., length(x), nterm + degree)
        if (any(outl)) {
            tt <- spline.des(knots, Boundary.knots[c(1,1)], degree+1, 0:1)
            newx[outl,] <- cbind(1, x[outl] - Boundary.knots[1]) %*% tt$design
        }
        if (any(outr)) {
            tt <- spline.des(knots, Boundary.knots[c(2,2)], degree+1, 0:1)
            newx[outr,] <- cbind(1, x[outr] - Boundary.knots[2]) %*% tt$design
        }
        if (any(inside <- !outside)) 
            newx[inside,] <-  spline.des(knots, x[inside], degree+1)$design
    }
    else newx <- spline.des(knots, x, degree+1, outer.ok=TRUE)$design

    # put missings back in so that the number of rows is right
    if (!all(keepx)) {
	temp <- matrix(NA, length(keepx), ncol(newx))
	temp[keepx,] <- newx
        newx <- temp
        }

    nvar <- ncol(newx)   #should be nterm + degree
    dmat <- diag(nvar)
    dmat <- apply(dmat, 2, diff, 1, 2) 
    dmat <- t(dmat) %*% dmat

    if (intercept) xnames <-paste('ps(', xname, ')', 1:nvar, sep='')
    else {
        newx <- newx[,-1, drop=FALSE]
        dmat <- dmat[-1,-1, drop=FALSE]    # rows corresponding to the 0 coef
        xnames <-paste('ps(', xname, ')', 1+ 2:nvar, sep='')
    }

    if (!penalty) {
        attributes(newx) <- c(attributes(newx), list(intercept=intercept,
                                          nterm=nterm,
                                          Boundary.knots=Boundary.knots))
        if (is.R()) class(newx) <- "pspline"
        else        oldClass(newx) <- "pspline"
        return(newx)
    }

    pfun <- function(coef, theta, n, dmat) {
	if (theta >=1) list(penalty= 100*(1-theta), flag=TRUE)
	else {
	    if (theta <= 0) lambda <- 0 
	    else lambda <- theta / (1-theta)
	    list(penalty= c(coef %*% dmat %*% coef) * lambda/2,
		 first  = c(dmat %*% coef) * lambda ,
		 second = c(dmat * lambda),
		 flag=FALSE
		 )
	    }
        }	

    printfun <- function(coef, var, var2, df, history, cbase) {
	test1 <- coxph.wtest(var, coef)$test
	# cbase contains the centers of the basis functions
	#   do a weighted regression of these on the coefs to get a slope
	xmat <- cbind(1, cbase)
	xsig <- coxph.wtest(var, xmat)$solve   # V X , where V = g-inverse(var)
	# [X' V X]^{-1} X' V
	cmat <- coxph.wtest(t(xmat)%*% xsig, t(xsig))$solve[2,]  
        linear <- sum(cmat * coef)
	lvar1  <- c(cmat %*% var %*% cmat)
	lvar2  <- c(cmat %*% var2%*% cmat)
	test2 <- linear^2 / lvar1
	# the "max(.5, df-1)" below stops silly (small) p-values for a
	#  chisq of 0 on 0 df, when using AIC gives theta near 1
	cmat <- rbind(c(linear, sqrt(lvar1), sqrt(lvar2), 
			test2, 1, 1-pchisq(test2, 1)),
		      c(NA, NA, NA, test1-test2, df-1, 
			1-pchisq(test1-test2, max(.5,df-1))))
	dimnames(cmat) <- list(c("linear", "nonlin"), NULL)
	nn <- nrow(history$thetas)
	if (length(nn)) theta <- history$thetas[nn,1]
	else  theta <- history$theta
	list(coef=cmat, history=paste("Theta=", format(theta)))
	}

    if (is.R()) {
	# The printfun needs to remember the spline's knots,
	#  but I don't need (or want) to carry around the entire upteen 
	#  variables defined here as an environment
	# Remove cbase from the arg list, and make it the environment
	formals(printfun) <- alist(coef=, var=, var2=, df=, history=)
	tempenv <- new.env(parent=asNamespace('survival'))
	assign('cbase',  knots[2:nvar] + (Boundary.knots[1] -knots[1]), 
               envir=tempenv)
	environment(printfun) <- tempenv
	}
    else {
	# Somewhat simpler in Splus, but because it depends on the 
	#  undocumented manner in which functions are stored, it might
	#  stop working one day
	printfun[[6]] <- knots[2:nvar] + (Boundary.knots[1] - knots[1])
	}
	
    if (method=='fixed') {
	temp <- list(pfun=pfun,
		     printfun=printfun,
		     pparm=dmat,
		     diag =FALSE,
		     cparm=list(theta=theta),
		     varname=xnames,
		     cfun = function(parms, iter, old)
			         list(theta=parms$theta, done=TRUE))
	}
    else if (method=='df') {
	temp <- list(pfun=pfun,
		     printfun=printfun,
		     diag =FALSE,
		     cargs=('df'),
		     cparm=list(df=df, eps=eps, thetas=c(1,0),
		                dfs=c(1, nterm), guess=1 - df/nterm, ...),
		     pparm= dmat,
		     varname=xnames,
		     cfun = frailty.controldf)
	}

    else { # use AIC
	temp <- list(pfun=pfun,
		     printfun=printfun,
		     pparm=dmat,
		     diag =FALSE,
		     cargs = c('neff', 'df', 'plik'),
		     cparm=list(eps=eps, init=c(.5, .95), 
		                lower=0, upper=1, ...),
		     varname=xnames,
		     cfun = frailty.controlaic)
	}
    
    attributes(newx) <- c(attributes(newx), temp,
                          list(intercept=intercept, nterm=nterm,
                          Boundary.knots=Boundary.knots))
    if (is.R()) class(newx) <- c("pspline", 'coxph.penalty')
    else        oldClass(newx) <- 'coxph.penalty'
    newx
    }

makepredictcall.pspline <- function(var, call) {
    if (call[[1]] != as.name("pspline")) return(call)  #wrong phone number
    newcall <- call[1:2]  #don't let the user override anything
    at <- attributes(var)[c("nterm", "intercept", "Boundary.knots")]
    newcall[names(at)] <- at
    newcall
}
    
predict.pspline <- function(object, newx, ...) {
    if (missing(newx)) return(object)
    a <- c(list(x=newx, penalty=FALSE), 
           attributes(object)[c("intercept, Boundary.knots")])
    do.call("pspline", a)
}

# Given a pspline basis, recover x
psplineinverse <- function(x) {
    if (!inherits(x, "pspline")) 
        stop("Argment must be the result of a call to pspline")
    intercept <- attr(x, "intercept")
    knots <- attr(x, "knots")
    nknot <- length(knots)
    if (!intercept) {
        indx <- 1:(ncol(x)+1) + (nknot- (ncol(x) +1))/2
        as.vector(cbind(1-rowSums(x), x) %*% knots[indx])
    }
    else {
        indx <- 1:ncol(x) + (nknot - ncol(x))/2
        as.vector(x %*% knots)
    }
}

    

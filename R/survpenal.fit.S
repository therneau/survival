# 
# fit a penalized parametric model
#
survpenal.fit<- function(x, y, weights, offset, init, controlvals, dist, 
		       scale=0, nstrat=1, strata, pcols, pattr, assign,
			 parms=NULL) {

    iter.max <- controlvals$iter.max
    outer.max <- controlvals$outer.max
    eps <- controlvals$rel.tolerance
    toler.chol <- controlvals$toler.chol

    if (!is.matrix(x)) stop("Invalid X matrix ")
    n <- nrow(x)
    nvar <- ncol(x) 
    ny <- ncol(y)
    if (is.null(offset)) offset <- rep(0,n)
    if (missing(weights)|| is.null(weights)) weights<- rep(1.0,n)
    else if (any(weights<=0)) stop("Invalid weights, must be >0")

    # The strata() term in survreg signals one scale parameter is
    #  to be fit per strata.  Here strata contains the strata level of each
    #  subject (variable not needed for only one strata), nstrat= # of strata.
    # Set nstrat2 = the number of coefficients I need to fit (which is 0
    #  if the scale is pre-fixed).
    if (scale <0) stop("Invalid scale")
    if (scale >0 && nstrat >1) 
	    stop("Cannot have both a fixed scale and strata")
    if (nstrat>1 && (missing(strata) || length(strata)!= n))
	    stop("Invalid strata variable")
    if (nstrat==1) strata <- rep(1,n)
    if (scale >0)  nstrat2 <- 0
    else           nstrat2 <- nstrat

    if (is.character(dist)) {
	sd <- survreg.distributions[[dist]]
	if (is.null(sd)) stop ("Unrecognized distribution")
	}
    else sd <- dist
    dnum <- match(sd$name, c("Extreme value", "Logistic", "Gaussian"))
    if (is.na(dnum)) {
	# Not one of the three distributions built in to the C code
        #  We need to set up a callback routine
        #  This returns the 5 number distribution summary (see the density
        #  functions in survreg.distributions).  Interval censored obs require
        #  2 evals and all others 1, so the call to the routine will have n2
        #  values.
	dnum <- 4  # flag for the C routine
	n2 <- n + sum(y[,ny]==3)  
        
	#
        # Create an expression that will be evaluated by the C-code,
        #   but with knowledge of some current variables
        # In the R doc, this would be "body(function(z) {"
	#  in Splus (Chambers book):  "functionBody(function(z)"
	#  same action, different name.  Luckily 'quote' exists in both
	# We make very sure the result is the right type and length here,
	#  rather than in the C code, for simplicity.
	fdensity <- quote({
	    if (length(parms)) temp <- sd$density(z, parms)
            else               temp <- sd$density(z)
	    
	    if (!is.matrix(temp) || any(dim(temp) != c(n2,5)) ||
                !is.numeric(temp))
		    stop("Density function returned an invalid matrix")
	    as.vector(as.double(temp))
	    })
        }
    else {
        fdensity <-1  #dummy value for the .Call
        n2 <- n       #a dummy value for inclusion in rho
        }

    # This is a subset of residuals.survreg: define the first and second
    #   derivatives at z=0 for the 4 censoring types
    #   Used below for starting estimates
    derfun <- function(y, eta, sigma, density, parms) {
	ny <- ncol(y)
	status <- y[,ny]
	z <- (y[,1] - eta)/sigma
	dmat <- density(z,parms)
	dtemp<- dmat[,3] * dmat[,4]    #f'
	if (any(status==3)) {
	    z2 <- (y[,2] - eta)/sigma
	    dmat2 <- density(z2)
	    }
	else {
	    dmat2 <- matrix(0,1,5)   #dummy values
	    z2 <- 0
	    }
	tdenom <- ((status==0) * dmat[,2]) +
		  ((status==1) * 1 )       +
		  ((status==2) * dmat[,1]) +
		  ((status==3) * ifelse(z>0, dmat[,2]-dmat2[,2], 
		                             dmat2[,1] - dmat[,1]))
	tdenom <- 1/(tdenom* sigma)
	dg <- -tdenom   *(((status==0) * (0-dmat[,3])) +
			  ((status==1) * dmat[,4]) + 
			  ((status==2) * dmat[,3]) +
			  ((status==3) * (dmat2[,3]- dmat[,3])))

	ddg <- (tdenom/sigma)*(((status==0) * (0- dtemp)) +
			       ((status==1) * dmat[,5]) +
			       ((status==2) * dtemp) +
			       ((status==3) * (dmat2[,3]*dmat2[,4] - dtemp))) 
	list(dg = dg, ddg = ddg - dg^2)
	}
    status <- y[,ny]

    #
    # are there any sparse frailty terms?
    # 
    npenal <- length(pattr)  #total number of penalized terms
    if (npenal == 0 || length(pcols) != npenal)
	    stop("Invalid pcols or pattr arg")
    sparse <- sapply(pattr, function(x) !is.null(x$sparse) &&  x$sparse)
    if (sum(sparse) >1) stop("Only one sparse penalty term allowed")
        
    #
    # Create a marking vector for the terms, the same length as assign
    #    with pterms == 0=ordinary term, 1=penalized, 2=sparse,
    #    pindex = length of pcols = position in pterms
    # 
    # Make sure that pcols is a strict subset of assign, so that the
    #   df computation (and printing) can unambiguously decide which cols of
    #   X are penalized and which are not when doing "terms" like actions.
    # To make some downstream things easier, order pcols and pattr to be
    #   in the same relative order as the terms in 'assign' 
    #
    pterms <- rep(0, length(assign))
    names(pterms) <- names(assign)
    pindex <- rep(0, npenal)
    for (i in 1:npenal) {
	temp <- unlist(lapply(assign, function(x,y) (length(x) == length(y) &&
					     all(x==y)), pcols[[i]]))
	if (sparse[i]) pterms[temp] <- 2
	else pterms[temp] <- 1
	pindex[i] <- (seq(along.with=temp))[temp]
	}
    if ((sum(pterms==2) != sum(sparse)) || (sum(pterms>0) != npenal))
	    stop("pcols and assign arguments disagree")
    if (any(pindex != sort(pindex))) {
	temp <- order(pindex)
	pindex <- pindex[temp]
	pcols <- pcols[temp]
	pattr <- pattr[temp]
	}
    # ptype= 1 or 3 if a sparse term exists, 2 or 3 if a non-sparse exists
    ptype <- any(sparse) + 2*(any(!sparse))
        
    if (any(sparse)) {
	sparse.attr <- (pattr[sparse])[[1]]  #can't use [[sparse]] directly
	                                     # if 'sparse' is a T/F vector
	fcol <- unlist(pcols[sparse])
	if (length(fcol) > 1) stop("Sparse term must be single column")

	# Remove the sparse term from the X matrix
	frailx <- x[, fcol]
	x <- x[, -fcol, drop=FALSE]
	for (i in 1:length(assign)){
	    j <- assign[[i]]
	    if (j[1] > fcol) assign[[i]] <- j-1
	    }
	for (i in 1:npenal) {
	    j <- pcols[[i]]
	    if (j[1] > fcol) pcol[[i]] <- j-1
	    }

	frailx <- match(frailx, sort(unique(frailx)))
	nfrail <- max(frailx)
	nvar <- nvar - 1
        
	#Set up the callback for the sparse frailty term
	#  (At most one sparse term is allowed).  The calling code will
	#  first set 'coef1' to the current value of the sparse coefficients,
	#  then call the expression below.  It uses a separate context (Splus
        #  frame or R environment), so there is no conflict between that
        #  variable name and the rest of the code.  Thus, think of the below as
	#  a funcion of the temporary variable coef1 (current value found
	#  in the calling C code), theta1 (current value in the S code
	#  below, using calls to cfun), and fixed known values of pfun1 etc.
        # The expression will constantly replace components of "coxlist1". By
        #  creating it first, we assure the order of the components, again
        #  to make it simpler for the C code (it can grab the first component
        #  and know that that is 'coef', etc).
	#
	pfun1 <- sparse.attr$pfun
        coxlist1 <- list(coef=0, first=0, second=0, penalty=0, flag=F)
	f.expr1 <- quote({
	    if (is.null(extra1)) temp <- pfun1(coef1, theta1, n.eff)
	    else  temp <- pfun1(coef1, theta1, n.eff, extra1)

	    if (!is.null(temp$recenter)) 
		    coxlist1$coef <- coef1 - as.double(temp$recenter)
	    else    coxlist1$coef <- coef1
	    if (!temp$flag) {
		coxlist1$first <- -as.double(temp$first)
		coxlist1$second <- as.double(temp$second)
	        }
	    else {
		coxlist1$first <- double(nfrail)
		coxlist1$second <- double(nfrail)
		}
	    coxlist1$penalty <- -as.double(temp$penalty)
	    coxlist1$flag   <- as.logical(temp$flag)

	    # Make sure the list has exactly the right structure, so
	    #  the the C code can be simple.  The first line below is
            #  probably unnecessary (belt AND suspenders); the second is
            #  checking a possibly user-supplied penaly function
	    if (any(names(coxlist1) != c('coef', 'first', 'second', 'penalty',
			                 'flag'))) 
		    stop("Invalid coxlist1")
	    if (any(sapply(coxlist1, length) != c(rep(nfrail,3), 1, 1)))
		    stop("Incorrect length in coxlist1")
	    coxlist1
            })
	}
    else {   # no sparse terms
	frailx <- 0
	nfrail <- 0
	f.expr1 <- NULL  #dummy value
        pfun1 <- NULL    #dummy
	coxlist1 <- NULL # "
	}
    nvar2 <- nvar + nstrat2
    if (nvar2 ==0) {
        # There are no non-sparse coefficients, and no scale parameters
        #   A strange model, leading to an hmat with 0 columns.  The
        #   underlying C code will choke, since this case is not built in.
        stop("Cannot fit a model with no coefficients other than sparse ones")
        }

    # Now the non-sparse penalties
    #  There can be multiple penalized terms
    if (sum(!sparse) >0) {
	full.imat <- !all(unlist(lapply(pattr, function(x) x$diag)))
	ipenal <- (1:length(pattr))[!sparse]   #index for non-sparse terms
        if (full.imat) {
            coxlist2 <- list(coef=double(nvar), first=double(nvar), 
                    second= double(nvar^2), penalty=0.0, 
                             flag=rep(FALSE,nvar))
            length2 <- c(nvar, nvar, nvar*nvar, 1, nvar)
            }  
        else {
            coxlist2 <- list(coef=double(nvar), first=double(nvar),
                    second=double(nvar), penalty= 0.0, flag=rep(FALSE,nvar))
            length2 <- c(nvar, nvar, nvar, 1, nvar)
            }
        
	# The C code will set the variable coef2, containing the concatonation
	#  of all the non-sparse penalized coefs.  Think of the below as
	#  a function of coef (from the C code), thetalist (set further
	#  below), and unchanging variables such as pattr.
	f.expr2 <- quote({
	    pentot <- 0
	    newcoef <- coef2

	    for (i in ipenal) {
		pen.col <- pcols[[i]]
		tcoef <- coef2[pen.col]
		if (is.null(extralist[[i]]))
			temp <- ((pattr[[i]])$pfun)(tcoef, thetalist[[i]],
                                                    n.eff)
		else    temp <- ((pattr[[i]])$pfun)(tcoef, thetalist[[i]],
						n.eff,extralist[[i]])
		if (!is.null(temp$recenter))
		    newcoef[pen.col] <- tcoef -  temp$recenter

		if (temp$flag) coxlist2$flag[pen.col] <- TRUE
		else {
		    coxlist2$flag[pen.col] <- FALSE
		    coxlist2$first[pen.col] <- -temp$first
		    if (full.imat) {
			tmat <- matrix(coxlist2$second, nvar, nvar)
			tmat[pen.col,pen.col] <- temp$second
			coxlist2$second <- c(tmat)
		        }
		    else coxlist2$second[pen.col] <- temp$second
		    }
		pentot <- pentot - temp$penalty
	        }
	    coxlist2$penalty <- as.double(pentot)
	    coxlist2$coef <- newcoef
	    if (any(sapply(coxlist2, length) != length2)) 
		    stop("Length error in coxlist2")
 	    coxlist2
	    })
        }
    else {
	full.imat <- FALSE  # no non-sparse penalties
	length2 <- 0  #dummy value
	f.expr2 <- NULL
	coxlist2 <- NULL
	ipenal <- NULL
	}
    

    # Create the frame for penalized evaluation
    # In R new.env inherits everything, in Splus new.frame only has
    #   what I specify at this time.  The variable thetalist will
    #   be iterated below, so we need to remember to update it within
    #   the Splus rho each time we do!
    # The variables parms, sd, and n2 are used for fdensity evaluation
    #
    if (is.R()) rho <- new.env()  
    #Splus    else rho <- new.frame(list(pfun1=pfun1, theta1=NULL, extra1=NULL,
    #			   nfrail=nfrail, pcols=pcols, pattr=pattr,
    #			   length2=length2, full.imat=full.imat,
    #			   ipenal = ipenal, nvar=nvar,
    #			   coxlist1=coxlist1, coxlist2=coxlist2, 
    #                       sd=sd, parms=parms, n2=n2))

    #
    # "Unpack" the passed in paramter list, 
    #   and make the initial call to each of the external routines
    #
    cfun <- lapply(pattr, function(x) x$cfun)
    parmlist <- lapply(pattr, function(x,eps) c(x$cparm, eps2=eps), sqrt(eps))
    extralist<- lapply(pattr, function(x) x$pparm)
    iterlist <- vector('list', length(cfun))
    thetalist <- vector('list', length(cfun))
    printfun  <- lapply(pattr, function(x) x$printfun)
    extra1 <- NULL
    theta1 <- NULL
    for (i in 1:length(cfun)) {
	temp <- (cfun[[i]])(parmlist[[i]], iter=0)
	if (sparse[i]) {
	    assign('theta1', temp$theta, rho)
	    assign('extra1', extralist[[i]], rho)
	    }
	thetalist[[i]] <- temp$theta
	iterlist[[i]] <- temp
	}
#    if (!is.R()) {  # Splus support
#	assign('thetalist', thetalist, frame=rho)
#    	assign('extralist', extralist, frame=rho)
#	}

    #
    # Manufacture the list of calls to cfun, with appropriate arguments
    #
    temp1 <- c('x', 'coef', 'plik', 'loglik', 'status', 'neff',  'df', 'trH')
    temp2 <- c('frailx', 'fcoef', 'fit$loglik-fit$penalty',  'fit$loglik', 
               'status', 'n.eff')
    temp3 <- c('x[,pen.col]', 'coef[pen.col]', 'fit$loglik-fit$penalty',
	       'fit$loglik', 'status', 'n.eff')
    calls <- vector('expression', length(cfun))
    cargs <- lapply(pattr, function(x) x$cargs)
    for (i in 1:length(cfun)) {
	tempchar <- paste("(cfun[[", i, "]])(parmlist[[", i, "]], iter,",
			  "iterlist[[", i, "]]")
	temp2b <- c(temp2, paste('pdf[', i, ']'), paste('trH[', i, ']'))
	temp3b <- c(temp3, paste('pdf[', i, ']'), paste('trH[', i, ']'))
	if (length(cargs[[i]])==0) 
	    calls[i] <- parse(text=paste(tempchar, ")"))
	else {
	    temp <- match(cargs[[i]], temp1)
	    if (any(is.na(temp))) stop(paste((cargs[[i]])[is.na(temp)],
					    "not matched"))
	    if (sparse[i]) temp4 <- paste(temp2b[temp], collapse=',')
	    else           temp4 <- paste(temp3b[temp], collapse=',')
	    
	    calls[i] <- parse(text=paste(paste(tempchar,temp4,sep=','),')'))
	    }
        }
    need.df <- any(!is.na(match(c('df', 'trH'), unlist(cargs))))#do any use df?

    #
    # Last of the setup: create the vector of variable names
    #
    varnames <- dimnames(x)[[2]]
    for (i in 1:npenal) {
	if (!is.null(pattr[[i]]$varname))
		varnames[pcols[[i]]] <- pattr[[i]]$varname
        }
        
    nvar2 <- nvar + nstrat2
    nvar3 <- nvar2 + nfrail
    #
    # A good initial value of the scale turns out to be critical for successful
    #   iteration, in a surprisingly large number of data sets.
    # The best way we've found to get one is to fit a model with only the
    #   mean and the scale.  We also the loglik of the mean-only model in the
    #   result
    # Even this model needs starting guesses...
	yy <- ifelse(status !=3, y[,1], (y[,1]+y[,2])/2 )
    coef <- sd$init(yy, weights,parms)
    # We sometimes get into trouble with a small initial estimate of sigma,
    #  (the surface isn't SPD), but never with a large one.  Double it.
    if (scale >0) vars <- log(scale)
	else vars <- log(4*coef[2])/2  # init gives \sigma^2, I need log(sigma)
    coef <- c(coef[1], rep(vars, nstrat))
    # get a better initial value for the mean using the "glim" trick
    deriv <- derfun(y, yy, exp(vars), sd$density, parms)
    wt <-  -1*deriv$ddg*weights
    coef[1] <- sum(weights*deriv$dg + wt*(yy -offset)) / sum(wt)

    fit0 <- .Call(Csurvreg6,
		  iter = as.integer(20),
		  nvar = as.integer(1),
		  as.double(y),
		  as.integer(ny),
		  x = as.double(rep(1.0, n)),
		  as.double(weights),
		  as.double(offset),
		  coef= as.double(coef),
		  as.integer(nstrat2),
		  as.integer(strata),
		  as.double(eps),
		  as.double(toler.chol), 
		  as.integer(dnum),
		  fdensity,
		  rho)

    # The "effective n" of the model
    temp <-  mean(exp(fit0$coef[-1]))   #overall sd
    n.eff <- sd$variance(temp^2) * (solve(matrix(fit0$var,1+nstrat2)))[1,1]
    #if (!is.R()) assign('n.eff', n.eff, frame=rho)

    #
    # Fit the model with all covariates
    #   Start with initial values
    #
    if (is.numeric(init)) {
	if (length(init) == nvar) {
	    if (scale >0) init <- c(init, log(scale))
	    else          init <-c(rep(0, nfrail), init, fit0$coef[-1])
	    }
	else if (length(init) == nvar2)  init <- c(rep(0,nfrail), init)
	else if (length(init) != nvar3) 
		stop("Wrong length for inital values")
	if (scale >0) init <- c(init, log(scale))
	}
    else  {
	# The algebra behind the 'glim' trick just doesn't work here
	#  Use the intercept fit + zeros
	#    coef order = frailty, intercept, other covariates, sigmas
	init <- c(rep(0, nfrail), fit0$coef[1], rep(0, nvar-1), fit0$coef[-1])
	}

    #
    # Tack on the sigmas to "assign", so that the df component includes
    #   the sigmas
    if (nstrat2 >0) assign <- c(assign, list(sigma=(1+nvar):nvar2))
    iter2 <- 0
    iterfail <- NULL
    thetasave <- unlist(thetalist)
    for (iterx in 1:outer.max) {
	fit <- .Call(Csurvreg7,
		   iter = as.integer(iter.max),
		   as.integer(nvar),
		   as.double(y),
		   as.integer(ny),
		   as.double(x),
	           as.double(weights),
		   as.double(offset),
		   coef= as.double(init),
	           as.integer(nstrat2),
	           as.integer(strata),
		   as.double(eps),
	           as.double(toler.chol), 
		   as.integer(dnum),
                   fdensity, 
                   rho,
	           as.integer(ptype),
		   as.integer(full.imat),
		   as.integer(nfrail),
		   as.integer(frailx),
		   f.expr1, f.expr2)

	iter <- iterx
	iter2 <- iter2 + fit$iter
	if (fit$flag == 1000) iterfail <- c(iterfail, iter)

	if (nfrail >0) {
	    fcoef <- fit$coef[1:nfrail]
	    coef  <- fit$coef[nfrail + 1:nvar2]
	    }
	else coef <- fit$coef[1:nvar2]

        # We need to fetch back some of the results from the
        #  evaluation area of f.expr1 and f.expr2
	if (is.R()) {
	    if (nfrail >0) coxlist1 <- get('coxlist1', envir=rho)
	    if (ptype >1 ) coxlist2 <- get('coxlist2', envir=rho)
	    }
	#else {
	#    if (nfrail >0) coxlist1 <- get('coxlist1', frame=rho)
	#    if (ptype >1 ) coxlist2 <- get('coxlist2', frame=rho)
	#    }

	# If any penalties were infinite, the C code has made hdiag=1 out
	#  of self-preservation (avoid zero divides).  But such coefs are 
	#  guarranteed to be zero so the variance should be too.
	temp <- rep(FALSE, nvar2+nfrail)
	if (nfrail>0) temp[1:nfrail] <- coxlist1$flag
	if (ptype >1) temp[nfrail+ 1:nvar] <- coxlist2$flag
	hdiag <- ifelse(temp, 0, fit$hdiag)

	if (need.df) {
            #get the penalty portion of the second derive matrix
	    if (nfrail>0) temp1 <- coxlist1$second
	    else 	  temp1 <- 0
	    if (ptype>1)  {
		if (full.imat) {
		    temp2 <- matrix(0., nvar2, nvar2)
		    temp2[1:nvar, 1:nvar] <- coxlist2$second
		    }
		else  temp2 <- diag(c(coxlist2$second, rep(0, nstrat2)))
		}
	    else          temp2 <- 0
					
	    dftemp <-coxpenal.df(matrix(fit$hmat, ncol=nvar2),  
			         matrix(fit$hinv, ncol=nvar2), hdiag, 
				 assign, ptype, nvar2,
		                 temp1, temp2, pindex[sparse])
	    df <- dftemp$df
	    var  <- dftemp$var
	    var2 <- dftemp$var2
	    pdf <- df[pterms>0]	          # df's for penalized terms
	    trH <- dftemp$trH[pterms>0]   # trace H 
	    }

	#
	# Call the control function(s)
	#
	done <- TRUE
	for (i in 1:length(cfun)) {
	    pen.col <- pcols[[i]]
	    temp <- eval(calls[i])
	    if (sparse[i]) assign('theta1', temp$theta, rho)
	    thetalist[[i]] <- temp$theta
	    iterlist[[i]] <- temp
	    done <- done & temp$done
    	    }
	if (done) break
	#if (!is.R()) assign('thetalist', thetalist, frame=rho)

	# 
	# Choose starting estimates for the next iteration
	#
	if (iter==1) {
	    init <- coefsave <- fit$coef
	    thetasave <- cbind(thetasave, unlist(thetalist))
	    }
	else {
	    temp <- unlist(thetalist)
	    coefsave <- cbind(coefsave, fit$coef)
	    # temp = next guess for theta
	    # *save = prior thetas and the resultant fits
	    # choose as initial values the result for the closest old theta
	    howclose <- apply((thetasave-temp)^2,2, sum)
	    which <- min((1:iter)[howclose==min(howclose)])
	    init <- coefsave[,which]
	    thetasave <- cbind(thetasave, temp)
	    }
        }   #end of the iteration loop
        
   if (!need.df) {  #didn't need it iteration by iteration, but do it now
        #get the penalty portion of the second derive matrix
	if (nfrail>0) temp1 <- coxlist1$second
	else 	      temp1 <- 0
	if (ptype>1)  {
		if (full.imat) {
		    temp2 <- matrix(0., nvar2, nvar2)
		    temp2[1:nvar, 1:nvar] <- coxlist2$second
		    }
		else  temp2 <- diag(c(coxlist2$second, rep(0, nstrat2)))
		}
	else  temp2 <- 0
					
	dftemp <-coxpenal.df(matrix(fit$hmat,ncol=nvar2),  
			     matrix(fit$hinv,ncol=nvar2),  hdiag, 
		             assign, ptype, nvar2, 
		             temp1, temp2, pindex[sparse])
	df <- dftemp$df
	trH <- dftemp$trH
	var <- dftemp$var
	var2  <- dftemp$var2
        }
        
    if (iter.max >1 && length(iterfail)>0)
	    warning(paste("Inner loop failed to coverge for iterations", 
			  paste(iterfail, collapse=' ')))
    which.sing <- (hdiag[nfrail + 1:nvar] ==0)
    coef[which.sing] <- NA

    names(iterlist) <- names(pterms[pterms>0])
    cname <- varnames
    cname <- c(cname, rep("Log(scale)", nstrat2))
    dimnames(var) <- list(cname, cname)
    names(coef) <- cname

    if (nfrail >0) {
	lp <- offset + fcoef[frailx]
	lp <- lp + x %*%coef[1:nvar] 
	list(coefficients  = coef,
	     icoef = fit0$coef,
	     var    = var,
	     var2   = var2,
	     loglik = c(fit0$loglik, fit$loglik- fit$penalty),
	     iter   = c(iter, iter2),
	     linear.predictors = as.vector(lp),
	     frail = fcoef,
	     fvar  = dftemp$fvar,
	     df = df, 
	     penalty= c(fit0$penalty, -fit$penalty),
	     pterms = pterms, assign2=assign,
	     history= iterlist,
	     printfun=printfun,
	     score = fit$u)
	}
    else {  #no sparse terms
	list(coefficients  = coef,
	     icoef = fit0$coef,
	     var    = var,
	     var2   = var2,
	     loglik = c(fit0$loglik, fit$loglik- fit$penalty),
	     iter   = c(iter, iter2),
	     linear.predictors = as.vector(x%*%coef[1:nvar]),
	     df = df, df2=dftemp$df2,
	     penalty= c(0, -fit$penalty),
	     pterms = pterms, assign2=assign,
	     history= iterlist,
	     printfun= printfun,
	     score = fit$u)
	}
    }

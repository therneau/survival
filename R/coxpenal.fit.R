#
# General penalized likelihood
#
coxpenal.fit <- function(x, y, strata, offset, init, control,
			weights, method, rownames, 
			pcols, pattr, assign) {
    eps <- control$eps
    n <-  nrow(y)
    if (is.matrix(x)) nvar <- ncol(x)
    else  if (length(x)==0) stop("Must have an X variable")
    else nvar <-1

    if (missing(offset) || is.null(offset)) offset <- rep(0,n)
    if (missing(weights)|| is.null(weights))weights<- rep(1,n)
    else {
	if (any(weights<=0)) stop("Invalid weights, must be >0")
	}

    # Get the list of sort indices, but don't sort the data itself
    if (ncol(y) ==3) {
	if (length(strata) ==0) {
	    sorted <- cbind(order(-y[,2], y[,3]), 
			    order(-y[,1]))
	    newstrat <- n
	    }
	else {
	    sorted <- cbind(order(strata, -y[,2], y[,3]),
			    order(strata, -y[,1]))
	    newstrat  <- cumsum(table(strata))
	    }
	status <- y[,3]
	andersen <- TRUE
#	routines <- paste('agfit5', c('a', 'b', 'c'), sep='_')
        routines <- list(Cagfit5a, Cagfit5b, Cagfit5c)
        }
    else {
	if (length(strata) ==0) {
	    sorted <- order(-y[,1], y[,2])
	    newstrat <- n
	    }
	else {
	    sorted <- order(strata, -y[,1], y[,2])
	    strata <- (as.numeric(strata))[sorted]
	    newstrat <-  cumsum(table(strata))
	    }
	status <- y[,2]
	andersen <- FALSE
#	routines <- paste('coxfit5', c('a', 'b', 'c'), sep='_')
        routines <- list(Ccoxfit5a, Ccoxfit5b, Ccoxfit5c)
        }

    n.eff <- sum(y[,ncol(y)])  #effective n for a Cox model is #events
    #
    # are there any sparse frailty terms?
    # 
    npenal <- length(pattr)
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
    ## can't compute assign attribute in R without terms
    ## if (missing(assign)) assign <- attr(x, 'assign')[-1]
    ##Remove 'intercept'
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

    ## Make sure these get defined <TSL>
    f.expr1<-function(coef) NULL
    f.expr2<-function(coef) NULL


    if (any(sparse)) {
	sparse.attr <- (pattr[sparse])[[1]]  #can't use [[sparse]] directly
	                                     # if 'sparse' is a T/F vector
	fcol <- unlist(pcols[sparse])
	if (length(fcol) > 1) stop("Sparse term must be single column")

	# Remove the sparse term from the X matrix
	xx <- x[, -fcol, drop=FALSE]
	for (i in 1:length(assign)){
	    j <- assign[[i]]
	    if (j[1] > fcol) assign[[i]] <- j-1
	    }
	for (i in 1:npenal) {
	    j <- pcols[[i]]
	    if (j[1] > fcol) pcols[[i]] <- j-1
	    }

	frailx <- x[, fcol]
	frailx <- match(frailx, sort(unique(frailx)))
	nfrail <- max(frailx)
	nvar <- nvar - 1

	#Set up the callback for the sparse frailty term
	pfun1 <- sparse.attr$pfun
        ### In R we use a function and eval() it, not an expression
	f.expr1 <- function(coef){
	    coxlist1$coef <- coef 
	    if (is.null(extra1)) temp <- pfun1(coef, theta1, n.eff)
	    else  temp <- pfun1(coef, theta1, n.eff, extra1)

	    if (!is.null(temp$recenter)) 
		    coxlist1$coef <- coxlist1$coef - as.double(temp$recenter)
	    if (!temp$flag) {
		coxlist1$first <- -as.double(temp$first)
		coxlist1$second <- as.double(temp$second)
	        }
	    coxlist1$penalty <- -as.double(temp$penalty)
	    coxlist1$flag   <- as.logical(temp$flag)
	    if (any(sapply(coxlist1, length) != c(rep(nfrail,3), 1, 1)))
		    stop("Incorrect length in coxlist1")
	    coxlist1
        }
        if (!is.null(getOption("survdebug"))) debug(f.expr1)
        
	coxlist1 <- list(coef=double(nfrail), first=double(nfrail), 
			 second=double(nfrail), penalty=0.0, flag=FALSE)
        ## we pass f.expr1 in as an argument in R
	##.C("init_coxcall1", as.integer(sys.nframe()), expr1)
    }
    else {
	xx <- x
	frailx <- 0
	nfrail <- 0
	}

    # Now the non-sparse penalties
    if (sum(!sparse) >0) {
	full.imat <- !all(unlist(lapply(pattr, function(x) x$diag)))
	ipenal <- (1:length(pattr))[!sparse]   #index for non-sparse terms
	f.expr2 <- function(coef){
            coxlist2$coef<-coef ##<TSL>
	    pentot <- 0
	    for (i in ipenal) {
		pen.col <- pcols[[i]]
		coef <- coxlist2$coef[pen.col]
		if (is.null(extralist[[i]]))
			temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]],n.eff)
		else    temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]],
						n.eff,extralist[[i]])
		if (!is.null(temp$recenter))
		    coxlist2$coef[pen.col] <- coxlist2$coef[pen.col]- 
			                               temp$recenter
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
	    if (any(sapply(coxlist2, length) != length2)) 
		    stop("Length error in coxlist2")
	    coxlist2
        }
        if (!is.null(getOption("survdebug")))
            debug(f.expr2)
	if (full.imat) {
	    coxlist2 <- list(coef=double(nvar), first=double(nvar), 
		    second= double(nvar*nvar), penalty=0.0, flag=rep(FALSE,nvar))
	    length2 <- c(nvar, nvar, nvar*nvar, 1, nvar)
	    }  
	else {
	    coxlist2 <- list(coef=double(nvar), first=double(nvar),
		    second=double(nvar), penalty= 0.0, flag=rep(FALSE,nvar))
	    length2 <- c(nvar, nvar, nvar, 1, nvar)
	    }
        ## in R, f.expr2 is passed as an argument later
	##.C("init_coxcall2", as.integer(sys.nframe()), expr2)
        }
    else full.imat <- FALSE

    #
    # Set up initial values for the coefficients
    #  If there are no sparse terms, finit is set to a vector of length 1
    #  rather than length 0, just to stop some "zero length" errors for
    #  later statements where fcoef is saved (but not used)
    #
    if (nfrail >0) finit <- rep(0,nfrail)
    else finit <- 0
    if (!missing(init) && !is.null(init)) {
	if (length(init) != nvar) {
	    if (length(init) == (nvar+nfrail)) {
		finit <- init[-(1:nvar)]
		init  <- init[1:nvar]
		}
	    else stop("Wrong length for inital values")
	    }
	}
    else init <- double(nvar)

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
    for (i in 1:length(cfun)) {
	temp <- (cfun[[i]])(parmlist[[i]], iter=0)
	if (sparse[i]) {
	    theta1 <- temp$theta
	    extra1 <- extralist[[i]]
	    }
	thetalist[[i]] <- temp$theta
	iterlist[[i]] <- temp
	}

    #
    # Manufacture the list of calls to cfun, with appropriate arguments
    #
    ## Amazingly, all this works in R, so I don't need to understand it.
    ##
    temp1 <- c('x', 'coef', 'plik', 'loglik', 'status', 'neff', 'df', 'trH')
    temp2 <- c('frailx', 'coxfit$fcoef', 'loglik1',  'coxfit$loglik', 'status',
	       'n.eff')
    temp3 <- c('xx[,pen.col]', 'coxfit$coef[pen.col]','loglik1',
	       'coxfit$loglik', 'status', 'n.eff')
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
    varnames <- dimnames(xx)[[2]]
    for (i in 1:length(cfun)) {
	if (!is.null(pattr[[i]]$varname))
		varnames[pcols[[i]]] <- pattr[[i]]$varname
        }

    ## need the current environment for callbacks
    rho<-environment()
    
    #
    # Have C store the data, and get the loglik for beta=initial, frailty=0
    #
    coxfit <- .C(routines[[1]],
                       as.integer(n),
                       as.integer(nvar), 
                       as.double(y),
                       x= as.double(xx) ,
                       as.double(offset),
                       as.double(weights),
		       as.integer(newstrat),
		       as.integer(sorted-1),
                       means= double(nvar),
                       coef= as.double(init),
                       u = double(nvar),
		       loglik=double(1),
		       as.integer(method=='efron'),
		       as.integer(ptype),
		       as.integer(full.imat),
		       as.integer(nfrail),
		       as.integer(frailx),
                 #R callback additions
                 f.expr1,f.expr2,rho)

    loglik0 <- coxfit$loglik
    means   <- coxfit$means

    #
    #  Now for the actual fit
    #
    iter2 <- 0
    iterfail <- NULL
    thetasave <- unlist(thetalist)
    for (outer in 1:control$outer.max) {
        coxfit <- .C(routines[[2]], 
		        iter=as.integer(control$iter.max),
			as.integer(n),
			as.integer(nvar),
		        as.integer(newstrat),
			coef = as.double(init),
		        u    = double(nvar+nfrail),
			hmat = double(nvar*(nvar+nfrail)),
			hinv = double(nvar*(nvar+nfrail)),
			loglik = double(1),
			flag = integer(1),
			as.double(control$eps),
		        as.double(control$toler.chol),
			as.integer(method=='efron'),
			as.integer(nfrail),
		        fcoef = as.double(finit),
			fdiag = double(nfrail+nvar),
                     ## R additions
                     f.expr1,f.expr2,rho)

	iter <- outer
	iter2 <- iter2 + coxfit$iter
	if (coxfit$iter >=control$iter.max) iterfail <- c(iterfail, iter)

	# If any penalties were infinite, the C code has made fdiag=1 out
	#  of self-preservation (0 divides).  But such coefs are guarranteed
	#  zero so the variance should be too.)
	temp <- rep(FALSE, nvar+nfrail)
	if (nfrail>0) temp[1:nfrail] <- coxlist1$flag
	if (ptype >1) temp[nfrail+ 1:nvar] <- coxlist2$flag
	fdiag <- ifelse(temp, 0, coxfit$fdiag)

	if (need.df) {
            #get the penalty portion of the second derive matrix
	    if (nfrail>0) temp1 <- coxlist1$second
	    else 	  temp1 <- 0
	    if (ptype>1)  temp2 <- coxlist2$second
	    else          temp2 <- 0
					
	    dftemp <-coxpenal.df(matrix(coxfit$hmat, ncol=nvar),  
			         matrix(coxfit$hinv, ncol=nvar), fdiag, 
				 assign, ptype, nvar,
		                 temp1, temp2, pindex[sparse])
	    df <- dftemp$df
	    var  <- dftemp$var
	    var2 <- dftemp$var2
	    pdf <- df[pterms>0]	          # df's for penalized terms
	    trH <- dftemp$trH[pterms>0]   # trace H 
	    }

	if (nfrail >0)  penalty <- -coxlist1$penalty
	else            penalty <- 0
	if (ptype >1) penalty <- penalty - coxlist2$penalty
	loglik1 <- coxfit$loglik + penalty  #C code returns PL - penalty
	if (iter==1) penalty0 <- penalty

	#
	# Call the control function(s)
	#
	done <- TRUE
	for (i in 1:length(cfun)) {
	    pen.col <- pcols[[i]]
	    temp <- eval(calls[i])
	    if (sparse[i]) theta1 <- temp$theta
	    thetalist[[i]] <- temp$theta
	    iterlist[[i]] <- temp
	    done <- done & temp$done
    	    }
	if (done) break

	# 
	# Choose starting estimates for the next iteration
	#
	if (iter==1) {
	    init <- coefsave <- coxfit$coef
	    finit <- fsave   <- coxfit$fcoef
	    thetasave <- cbind(thetasave, unlist(thetalist))
	    }
	else {
	    # the "as.vector" removes names, dodging a bug in Splus5.1
	    temp <- as.vector(unlist(thetalist))
	    coefsave <- cbind(coefsave, coxfit$coef)
	    fsave    <- cbind(fsave, coxfit$fcoef)
	    # temp = next guess for theta
	    # *save = prior thetas and the resultant fits
	    # choose as initial values the result for the closest old theta
	    howclose <- apply((thetasave-temp)^2,2, sum)
	    which <- min((1:iter)[howclose==min(howclose)])
	    if (nvar>0)   init <- coefsave[,which]
	    if (nfrail>0) finit<- fsave[,which]
	    thetasave <- cbind(thetasave, temp)
	    }
        }

    # release the memory
    expect <- .C(routines[[3]], as.integer(n),
		             as.integer(nvar),
		             as.integer(newstrat),
		             as.integer(method=='efron'),
		             expect= double(n))$expect

    if (!need.df) {  #didn't need it iteration by iteration, but do it now
        #get the penalty portion of the second derive matrix
	if (nfrail>0) temp1 <- coxlist1$second
	else 	      temp1 <- 0
	if (ptype>1)  temp2 <- coxlist2$second
	else          temp2 <- 0
					
	dftemp <-coxpenal.df(matrix(coxfit$hmat,ncol=nvar),  
			     matrix(coxfit$hinv,ncol=nvar),  fdiag, 
		             assign, ptype, nvar, 
		             temp1, temp2, pindex[sparse])
	df <- dftemp$df
	trH <- dftemp$trH
	var <- dftemp$var
	var2  <- dftemp$var2
        }

    if (control$iter.max >1 && length(iterfail)>0)
	    warning(paste("Inner loop failed to coverge for iterations", 
			  paste(iterfail, collapse=' ')))
    which.sing <- (fdiag[nfrail + 1:nvar] ==0)
    
    coef <- coxfit$coef
    names(coef) <- varnames
    coef[which.sing] <- NA
    resid <- double(n)
    resid <- status - expect
    names(resid) <- rownames

    names(iterlist) <- names(pterms[pterms>0])
    if (nfrail >0) {
	lp <- offset + coxfit$fcoef[x[,fcol]]
	if (nvar >0) {   #sparse frailties and covariates
	    lp <- lp + x[,-fcol,drop=FALSE] %*%coxfit$coef - 
                sum(means*coxfit$coef)
	    list(coefficients  = coef,
		 var    = var,
		 var2   = var2,
		 loglik = c(loglik0, loglik1),
		 iter   = c(iter, iter2),
		 linear.predictors = as.vector(lp),
		 residuals = resid,
		 means = means,
                 concordance= survConcordance.fit(y, lp, strata, weights),
		 method= c('coxph.penal', 'coxph'),
		 frail = coxfit$fcoef,
		 fvar  = dftemp$fvar,
		 df = df, df2=dftemp$df2, 
		 penalty= c(penalty0, penalty),
		 pterms = pterms, assign2=assign,
		 history= iterlist,
		 coxlist1=coxlist1, 
		 printfun=printfun)
	    }
	else {  #sparse frailties only
	    list( loglik = c(loglik0, loglik1),
		 iter   = c(iter, iter2),
		 linear.predictors = as.vector(lp),
		 residuals = resid,
		 means = means,
                 concordance= survConcordance.fit(y, lp, strata, weights),
		 method= c('coxph.penal', 'coxph'),
		 frail = coxfit$fcoef,
		 fvar  = dftemp$fvar,
		 df = df, df2=dftemp$df2, 
		 penalty = c(penalty0, penalty), 
		 pterms = pterms, assign2=assign,
		 history= iterlist,
		 printfun=printfun)
	    }
         }
    else {  #no sparse terms
        lp <- offset + as.vector(x%*%coxfit$coef) - sum(means*coxfit$coef)
	list(coefficients  = coef,
	     var    = var,
	     var2   = var2,
	     loglik = c(loglik0, loglik1),
	     iter   = c(iter, iter2),
	     linear.predictors = lp,
	     residuals = resid,
	     means = means,
             concordance= survConcordance.fit(y, lp, strata, weights),
	     method= c('coxph.penal', 'coxph'),
	     df = df, df2=dftemp$df2,
	     penalty= c(penalty0, penalty), 
	     pterms = pterms, assign2=assign,
	     history= iterlist,
	     coxlist2=coxlist2,
	     printfun= printfun)
	}
    }

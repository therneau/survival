#
# This is common code for survpenal.fit and coxpenal.fit.
#  It's all the bookkeeping to set up the penalized callbacks
# This code is in development, and not yet used by anything, the if(FALSE)
#  keeps it out of the distribution
if(FALSE){
survcallback <- function(pcols, pattr, assign, x) {
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
	#  first set 'coef' to the current value of the sparse coefficients,
	#  then call the expression below.  It uses a separate context (Splus
        #  frame or R environment), so there is no conflict between that
        #  variable name and the rest of the code.  Thus, think of the below as
	#  a funcion of the temporary variable coef (current value found
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
        
	# The C code will set the variable coef, containing the concatonation
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
    
    list(f.expr1=f.expr1,   f.expr2=f.expr2, 
	 coxlist1=coxlist1, coxlist2=coxlist2,
	 full.imat=full.imat, ipenal=ipenal, length2=length2,
	 pfun1=pfun1, pindex=pindex, pcols=pcols, pattr=pattr,
	 sparse=sparse, frailx=frailx, nfrail=nfrail,
	 nvar=nvar)
    }
}

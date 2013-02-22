# Aalen's additive regression model
#  Originally, this tried to call coxph with certain options.
#  But we found the passing ... to a model method just doesn't work (for
#   optional things like weights).  So the first portion of this is
#   essentially coxph, to set up for coxph.detail.
# For distribution, the "variance" test is omitted. Not all aspects are
#   yet supported by the downstream printing.
#
aareg <- function(formula, data, weights, subset, na.action,
                  qrtol=1e-7, nmin, dfbeta=FALSE, taper=1,
		  test = c('aalen', 'variance', 'nrisk'),
		  model=FALSE, x=FALSE, y=FALSE) {
    call <- match.call()

    m <- match.call(expand.dots=FALSE)
    temp <- c("", "formula", "data", "weights", "subset", "na.action")
    m <- m[ match(temp, names(m), nomatch=0)]
    special <- c("strata", "cluster")
    Terms <- if(missing(data)) terms(formula, special)
	     else              terms(formula, special, data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())

    test <- match.arg(test)     #check for legal argument

    # Now grab the terms that we need
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    offset<- attr(Terms, "offset")
    tt <- length(offset)
    offset <- if(tt == 0)
		    rep(0, nrow(Y))
	      else if(tt == 1)
		      m[[offset]]
	      else {                 #multiple offset terms!  add them
		    ff <- m[[offset[1]]]
		    for(i in 2:tt)
			    ff <- ff + m[[offset[i]]]
		    ff
		    }

    # Create an X matrix, to feed to the coxdetail routine.
    attr(Terms,"intercept")<- 1  # force no intercept
    strats <- attr(Terms, "specials")$strata
    cluster<- attr(Terms, "specials")$cluster
    dropx <- NULL
    if (length(cluster)) {
	dfbeta <- TRUE
	tempc <- untangle.specials(Terms, 'cluster', 1:10)
	ord <- attr(Terms, 'order')[tempc$terms]
	if (any(ord>1)) stop ("Cluster can not be used in an interaction")
	cluster <- strata(m[,tempc$vars], shortlabel=TRUE)  #allow multiples
	cluster <- as.numeric(cluster) #labels don't matter, and processing
	                               # is a bit faster without them
	dropx <- tempc$terms
	}
    else cluster <- 1:nrow(m)

    # Adding strata, when there is a coefficent per death, is identical
    #  to doing a totally separate fit per group.
    # Using "factor(group)" to get multiple baselines is likely what the
    #  user wants.  However, because we have not processed the strata
    #  statement (taken it out of X, and created the 'newstrat' of coxph)
    #  it will act just like a factor.
    # I've changed my mind multiple times on commenting out the line below.
    #  Computationally identical to factor() -- is an error message or not
    #  an error message the greater source of confusion to a user?
    if (length(strats)) {
       stop("Strata terms not allowed")
       }

    if (length(dropx)) X <- model.matrix(Terms[-dropx], m)[,-1,drop=FALSE]
    else               X <- model.matrix(Terms, m)[,-1,drop=FALSE]
    nvar <- ncol(X)
    nused<- nrow(X)
    weights <- model.extract(m, 'weights')
    if (length(weights)==0) weights <- rep(1.0, nused)

    type <- attr(Y, "type")
    if (type!='right' && type!='counting')
	stop(paste("Aalen model doesn't support \"", type,
			  "\" survival data", sep=''))
   
    # Get the peices that I need from the coxdetail routine
    #  1. It expects a "counting process" type of Y
    if (ncol(Y)==2) {
	mintime <- min(Y[,1])
	if (mintime < 0) Y <- cbind( 2*mintime -1, Y)
	else 	Y <- cbind(-1,Y)
	}
    # Because there are no strata, the number of unique death times is the
    #  number that will be in the output structures
    times <- as.vector(Y[,2])  # toss the labels away
    status<- as.vector(Y[,3])
    ndeath <- length(unique(times[status==1]))

    # Sort everything
    ord <- order(times, -status)
    times <- times[ord]
    status <- status[ord]
    weights <- weights[ord]
    if (x) saveX <- X
    X <- X[ord,,drop=FALSE]

    storage.mode(Y) <- 'double'
    ff <- .C(Ccoxdetail, as.integer(nused),
			  as.integer(nvar),
			  ndeath= as.integer(ndeath),
			  y = Y[ord,],
			  as.double(X),
			  index = as.integer(rep(0,nused)),
			  event2 = rep(1.0, nused),
			  weights = as.double(weights),
			  means= c(0., double(ndeath*nvar-1)),
			  u = double(ndeath*nvar),
			  i = double(ndeath*nvar*nvar),
	                  rmat = integer(ndeath*nused),
	                  nrisk2 = double(ndeath),
			  double(nvar*(3 + 2*nvar)) )
    # riskmat is an nused by ndeath 0/1 matrix showing who is present
    riskmat <- matrix(ff$rmat, nused, ndeath) 

    # Note that imat, as returned by coxdetail, is Var(X) * nevents.
    dt <- list(means= (matrix(ff$means,ndeath, nvar)),
	       var = aareg.taper(taper, array(ff$i, c(nvar, nvar, ndeath)), 
                                     ff$event2[1:ndeath]),
	       time = times[ff$index[1:ndeath]],
	       nrisk= ff$nrisk2,            #weighted # at risk
	       nevent=ff$event2[1:ndeath])  #weighted number of events

    # Set the number of deaths that will be used in the analysis
    #  This may be smaller than the curren "ndeath", due to small nrisk
    # The number of times may even be smaller, if imat is singular at that
    #   time point.
    if (missing(nmin)) nmin <- 3*nvar
    if (nvar==1)
        ndeath   <- sum(dt$nrisk>= nmin & c(dt$var)>0)
    else {
        ndeath <- sum(dt$nrisk >= nmin)
        if (ndeath >0) {
            while (1) {  #we expect very few iterations of this loop
                qri <- qr(dt$var[,,ndeath], tol=qrtol)
                if (qri$rank >= nvar) break   #not singular
                ndeath <- ndeath -1
                }
            }
        }
    if (ndeath<=1) 
	    stop("The threshold 'nmin' is too high, no model can be fit")
    
    # This matches the death times in the data set to the
    #  sorted list of unique death times.  "0" = not a death 
    index <- match(times, dt$time[1:ndeath], nomatch=0) * status
    deaths <- (status==1 & index >0)
    dindex <- index[deaths]  #for each death, a pointer into dt objects
    nevent <- length(dindex) #total number of events (ndeath = #unique times)

    if (length(cluster)) ncluster <- length(unique(cluster))
    else                 ncluster <- nused

    if (dfbeta) {
	dmat <- array(0.0, dim=c(ncluster, nvar+1, ndeath))
	# the resid marix has a row for each death, and nused cols
	#  each row has a "1" in it at the position of the death
	#  the yhat part is subtracted later
        resid <- rep(0., nevent*nused)
        resid[nevent*((1:nused)[deaths]-1) + 1:nevent] <- 1.0 
        resid <- matrix(resid, ncol=nused)  
	}

    #   Coefficient is the step in Aalen's plots
    #   If we keep one row of "coefficent" per death, then Aalen's
    # variance is coef *% t(coef), treating coef as a col vector.
    # If we kept one row per death, then the ndeath nvar by nvar variance
    # matrices would need to be kept too.  So keep 1 row per event.
    # Things like plot will end up accumlating.
    #   There is no such cheat for dfbeta: it is kept as "# unique deaths"
    # p by p matrices.
    #
    if (nvar==1)  { # special case of only 1 covariate
	means <- dt$means[dindex]
	nrisk <- dt$nrisk[dindex]
	xx <- (X[deaths] - means) * weights[deaths]

	v.inverse <- 1/dt$var[dindex] #for all time points
	twt  <- nrisk* 1/cbind(1+ means^2*v.inverse, v.inverse)
	coefficient <- v.inverse * xx / nrisk 
	# Note that ybar is always w_i/nrisk, since we are doing the 
        #  regressions one event at a time.
	b0 <- weights[deaths]/nrisk - means*coefficient

	if (dfbeta) {
            # We first create the nused * #events matrix, and then
            #  collapse it to be ncluster by n-unique-death-times
	    xx <- c(X) * riskmat[,dindex] # X repeated in each col, if at risk
	    predicted <- coefficient * t(xx) + b0*t(riskmat[,dindex]) 
	    resid <- resid - predicted  #nused cols, nvevent rows
            temp1 <- (resid * (t(xx) -means)/(nrisk*dt$var[dindex])) *
                       rep(weights, rep(nevent, nused)) 

	    # temp1[i,j] is the change in alpha at time i for subject j
	    # the "intercept dfbeta" is resid*wt/sum(wt) - xbar*temp1
	    temp0 <-  resid * outer(1/nrisk, weights) -  temp1 * means

	    # get the matrix, nused by 2, which is the influence of each 
	    #   subject on the test statistic.
	    #   This is a bit easier before collapsing
	    if (test=='nrisk') {
		test.dfbeta <- cbind(apply(temp0*nrisk, 2, sum),
				     apply(temp1*nrisk, 2, sum))
		}
	    else {
		test.dfbeta <- cbind(apply(temp0*twt[,1], 2, sum),
				     apply(temp1*twt[,2], 2, sum))
		}
		
            # Now collapse dfbeta, first on the deaths, and then on the cluster
            if (nevent > ndeath) {
                temp1 <- rowsum(temp1, times[deaths], reorder=FALSE)
                temp0 <- rowsum(temp0, times[deaths], reorder=FALSE)
                }
            dmat[,1,] <- rowsum(t(temp0), cluster[ord], reorder= FALSE)
            dmat[,2,] <- rowsum(t(temp1), cluster[ord], reorder =FALSE)
	    }

	# Compute the test statistic, including the intercept term
	# (Much of the code above was a litte easier to write without
	#  the intercept term in thec coef matrix, that below is easier
	#  with it in).
	coefficient <- cbind(b0,coefficient)
	if (test=='nrisk') {
	    temp <- coefficient*nrisk
	    test.statistic <- apply(temp,2,sum)
	    test.var  <- matrix(0.,2,2)
	    diag(test.var) <- apply(temp^2, 2, sum)
	    test.var[1,2] <- test.var[2,1] <- sum(temp[,1]*temp[,2])
	    }
	else {  # full V^{-1} and diag(V){-1} variance (Aalen) are the same
	    temp <- coefficient* twt
	    test.statistic <- apply(temp,2,sum)
	    test.var  <- matrix(0.,2,2)
	    diag(test.var) <- apply(temp^2, 2, sum)
	    test.var[1,2] <- test.var[2,1] <- sum(temp[,1]*temp[,2])
	    }
	}
    
    else { # 2 or more covariates
	coefficient <- matrix(0,nevent, nvar)
	twt <-  matrix(0, nevent, nvar+1)
	means <- dt$means[dindex,]  # vector of means, at each deatj
	nrisk <- dt$nrisk[dindex]
        dindex2 <- (1:nused)[deaths]  # row number of each death
	ybar <- weights[deaths]/nrisk
	test.var <- matrix(0.0, nvar, nvar)
	if (dfbeta) test.dfbeta <- matrix(0., nused, nvar+1)

	for (i in 1:nevent) {	    
            who <- riskmat[,dindex[i]] # 0/1 vector of who is at risk
	    wt <- weights* who
            xx <- who* (X- rep(means[i,], rep(nused, nvar))) # (X-Xbar)

	    # solve, and check for singularity
	    # Note that the increment to imat, as returned by
	    #   the coxph.detail function, is Var(X) * #events
	    #   and qri is intended to be the qr of V-inverse
            if (i==1 || dindex[i] != dindex[i-1]) {  #don't redo qr for ties
                qri <- qr(dt$var[,,dindex[i]], tol=qrtol)
		vmat <- qr.coef(qri, diag(nvar))
                twt[i,] <- nrisk[i] /c(1+ means[i,] %*% vmat %*% means[i,],
				       diag(vmat))
                }
            else twt[i,] <- twt[i-1,]
	    j <- dindex2[i]
	    coefficient[i,] <-qr.coef(qri, wt[j]*xx[j,]) / nrisk[i]
	    if (test=='variance') {
		temp <- wt[j]*xx[j,] 
		test.var <- test.var + outer(temp,temp)
		}

	    if (dfbeta) {
		resid[i, ] <- resid[i,]- c(ybar[i] + xx %*% c(coefficient[i,]))
		temp1 <- t(qr.coef(qri, t(resid[i,]* wt *xx)))/ nrisk[i]
		temp0 <- resid[i,]*wt/nrisk[i] - temp1%*% means[i,]

		if (test=='aalen') 
			test.dfbeta <- test.dfbeta + 
				       cbind(temp0, temp1) %*% diag(twt[i,]) 
		else if (test=='nrisk')
			test.dfbeta <- test.dfbeta + 
				        cbind(temp0, temp1)* nrisk[i] 
		else {
		    test.dfbeta[,-1] <- test.dfbeta[,-1] + resid[i,]* wt *xx
		    # There really isn't a definition for what weight to
		    #  put on the intercept in the "variance" weighting
		    #  (and who really cares about "testing the intercept"
		    #  anyway).  So use the twt one
		    test.dfbeta[,1] <- test.dfbeta[,1] + temp0*twt[i,1]
		    }
		dmat[,-1,dindex[i]] <- dmat[,-1, dindex[i]] +
                                       rowsum(temp1, cluster[ord], reorder=FALSE)
		dmat[,1, dindex[i]] <- dmat[,1,dindex[i]] +
			              rowsum(temp0, cluster[ord], reorder=FALSE)
		}
	    }
        temp <- apply(means*coefficient, 1, sum) # xbar * coef at time t
	b0 <- weights[deaths]/nrisk - temp
	coefficient <- cbind(b0,coefficient)

	# Note - the intercept is a part of the test statistic, even
	#   though it will always be ignored in the overall chisquare test
	if (test=='aalen') {
	    temp <- twt* coefficient
	    test.statistic <- colSums(temp)
	    test.var <- t(temp) %*% temp
	    }
	else if (test=='nrisk') {
	    temp <- coefficient * nrisk
	    test.statistic <- colSums(temp)
	    test.var <- t(temp) %*% temp
	    }
	else  {
	    xx <- weights[deaths]*(X[deaths,] - means[dindex,]) 
	    test.statistic <- apply(xx, 2, sum)
	    }
	}

    if (dfbeta) {
	# The model variance is sum( term[i]^2), i ranging over times,
	#   and each term an n by p matrix (one row per person)
	# The dfbeta one is essentially [sum(term[i])]^2
	#   the test.dfbeta matrix contains this sum over death times
	temp <- rowsum(test.dfbeta, cluster, reorder=FALSE)
	test.var2 <- t(temp) %*% temp
	}

    dimnames(coefficient) <- list(times[deaths], 
				c("Intercept", dimnames(X)[[2]]))
    names(test.statistic) <- c("Intercept", dimnames(X)[[2]])
    dimnames(twt) <- NULL

    ans <- list(n= c(nused, ndeath, length(dt$time)), times=times[deaths], 
		nrisk=dt$nrisk[dindex], 
		coefficient=coefficient, 
		test.statistic=test.statistic, test.var=test.var, test=test,
		tweight = twt, call=call) 

    if (dfbeta) {
	ans$dfbeta <- dmat
	ans$test.var2 <- test.var2
	}
    if (any(weights!=1)) ans$weights <- weights
#    if (ncluster < nused) ans$cluster <- as.numeric(cluster)
    na.action <- attr(m, "na.action")
    if (length(na.action)) ans$na.action <- na.action
    if (model) ans$model <- m
    else {
	if (x) ans$x <- saveX
	if (y) ans$y <- Y
	}

    if (is.R()) class(ans) <- 'aareg'
    else oldClass(ans) <- 'aareg'
    ans
    }

"[.aareg" <- function(x, ..., drop=FALSE) {
    if (!inherits(x, 'aareg')) stop ("Must be an aareg object")

    # There is a bug in Splus6, the "=FALSE" on drop is ignored as a
    #   default.  Add temporary hack
    drop <- FALSE

    i <- ..1
    if (is.matrix(x$coefficient)) {
	x$coefficient <- x$coefficient[,i, drop=drop]
	x$tweight     <- x$tweight[,i,drop=drop]
	}
    else stop("Subsripting impossible, coefficient component not a matrix")

    if (!is.null(x$dfbeta)){
	x$dfbeta <- x$dfbeta[,i,,drop=drop]
	x$test.var2 <- x$test.var2[i,i,drop=drop]
	}
    x$test.statistic <- x$test.statistic[i, drop=drop]
    x$test.var <- x$test.var[i,i,drop=drop]
    x
    }


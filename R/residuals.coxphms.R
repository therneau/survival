#
# Much of this function is a close copy of residuals.coxph.  The difference
#  is that we need to set up the X matrix, strata, and id differently, using
#  stacker. The model component of the object, or result of model.frame(object),
#  will be in the same row order as the original data; and if collapse=FALSE
#  that is exactly the order we want for the residuals. But the analysis
#  used a stacked data frame, which we need to recreate for any residual that
#  needs X or Y.  Then at the end put results in the right order.
residuals.coxphms <- function(object, type=c("martingale","score",
                                             "schoenfeld",
                                             "dfbeta", "dfbetas", "scaledsch"),
                              collapse=FALSE, 
                              weighted=(type %in% c("dfbeta", "dfbetas")),
                              na.action, ...) {

    omit <- object$na.action
    if (!is.null(omit) && !missing(na.action)){
        if (collapse && (na.action=="na.omit" || inherits(omit, "exclude")))
            stop("collapse and expansion of missing are mutually exclusive")
        if (na.action=="na.omit") class(omit) <- "omit"
        else if (na.action=="na.exclude") class(omit) <- "exclude"
        else stop("invalid na.action argument")
    }

    type <- match.arg(type)
    if (type=="deviance") 
        stop("deviance residuals not supported for multi-state")
    otype <- type
    if (type=='dfbeta' || type=='dfbetas') {
        otype <- type   # used for error messge
	type <- 'score'
	if (missing(weighted))
            weighted <- TRUE  # different default for this case
    }
    if (type=='scaledsch') type<-'schoenfeld'
    method <- object$method
    if (method=='exact' && (type=='score' || type=='schoenfeld'))
	stop(paste(otype, 'residuals are not available for the exact method'))

    # if type = martingale and collapse is false the computation is quick
    # it is also a common use case
    if (type=='martingale' && is.logical(collapse) && !collapse &&
          is.logical(weighted) && !weighted) {
        rr <- matrix(0, nrow=object$n, ncol= ncol(object$cmap))
        rr[object$rmap] <- object$residuals   # matrix subscript
        colnames(rr) <- colnames(object$cmap)
        if (!is.null(omit)) return(naresid(omit, rr))
        else return(rr)
    }

    # For all other cases we need the model frame
    mf <- model.frame(object)  # grab the model frame
    id <- model.extract(mf, "id")  # required for a coxphms models
    weights <- model.weights(mf)
    offset <- model.offset(mf)

    # and we will refer to these multiple times
    smap <- object$smap
    cmap <- object$cmap
    utran <- colnames(smap)[unique(smap[1,])]  # the unique transitions

    if (!is.logical(collapse)) {
        # they gave us a collapse vector
        nrow <- nrow(mf)
        if (length(collapse) != nrow) 
            stop("collapse vector not the same length as the model frame")
        cluster <- collapse
        collapse <- TRUE
    } else if (collapse) {
        cluster <- model.extract(mf, "cluster")
        if (is.null(cluster)) cluster <- id
    }        
    if (collapse & !inherits(cluster, "factor")) 
        cluster <- factor(cluster, unique(cluster)) # remember original order

    if (type=="martingale") { # finish up the martingale case
        rr <- matrix(0, nrow=object$n, ncol= ncol(object$cmap))
        rr[object$rmap] <- object$residuals   #matix subscript
        colnames(rr) <- colnames(object$cmap)
        if (weighted) rr <- weights*rr
        if (collapse) return(rowsum(rr, cluster, reorder=FALSE))
        else if (!is.na(omit)) return(naresid(omit, rr))
        else return(rr)
    }
            
    # For all other types we need to create the updated X and Y matrices using
    #  stacker.  First get X, Y, and strata
    # istrat will be the expanded strata that includes states
    Terms <- terms(object)
    strats <- attr(Terms, "specials")$strata

    Y <- model.response(mf)
    X <- model.matrix(object, mf)
    istrat <- attr(X, "strata")
    istate <- model.extract(mf, "istate")
    mcheck <- survcheck2(Y, id, istate)
    transitions <- mcheck$transitions
    istate <- mcheck$istate
    states <- mcheck$states

    xstack <- stacker(object$cmap, object$smap, as.integer(istate), 
                      X, Y, mf=mf, states=states)

    rkeep <- unique(xstack$rindex)
    mcheck<- survcheck2(Y[rkeep,], id[rkeep], istate[rkeep])

    x <- xstack$X  # the code in residuals.coxph uses lower case
    y <- xstack$Y
    n <- nrow(y)
    istrat <- xstack$strata
    if (collapse) cluster <- cluster[xstack$rindex]
    istate <- mcheck$istate[xstack$rindex]  #initial state
    if (length(offset)) offset <- offset[xstack$rindex]
    else offset <- double(n)  # all zero
    if (length(weights)) weights <- weights[xstack$rindex]
    else weights <- rep(1, n)
    
    ny <- ncol(y)
    status <- y[,ny,drop=TRUE]
    vv <- drop(object$naive.var)
    if (is.null(vv)) vv <- drop(object$var)

    nvar <- ncol(x)
    ord <- order(istrat, y[,ny-1], -status)
    newstrat <- c(diff(as.numeric(istrat[ord]))!=0 ,1)
    newstrat[n] <- 1  #C routines use a 'mark last obs of each strata' var

    # sort the data
    x <- x[ord,]
    y <- y[ord,]
    score <- exp(object$linear.predictors)[ord]
    istrat <- istrat[ord]
    if (ny==3) sort1 <- order(istrat, y[,1])
 
    #
    # Now I have gotton the data that I need-- do the work
    # 
    # Schoenfeld and scaled Schoenfeld are a residual per event
    #  neither collapse nor na.omit are relevant
    if (type=='schoenfeld') {
  	if (ny==2) {
 	    mintime <- min(y[,1])
 	    if (mintime < 0) y <- cbind(2*mintime -1, y)
 	    else             y <- cbind(-1,y)
        }
	temp <- .C(Ccoxscho, n=as.integer(n),
			    as.integer(nvar),
			    as.double(y),
			    resid=  as.double(x),
			    as.double(score * weights[ord]),
			    as.integer(newstrat),
			    as.integer(method=='efron'),
			    double(3*nvar) )

	deaths <- y[,3]==1

	if (nvar==1) rr <- temp$resid[deaths]
	else         rr <- matrix(temp$resid[deaths], ncol=nvar) #pick rows 
	if (weighted) rr <- rr * weights[deaths]

 	if (otype=='scaledsch') {
	    ndead <- sum(deaths)
	    coef <- ifelse(is.na(object$coefficients), 0, object$coefficients)
            if (nvar==1) rr <- rr * vv * ndead + coef
	    else rr <- drop(rr %*% vv *ndead + rep(coef, each=nrow(rr)))
        }
        
        # We can't return an array of dimension (time, covariate, state), since
        #  each transition has unique event times.  So transitions are strung
        #  out one after the other. Strata, if present, are within transition.
        #  (I've never found a use for the transition or strata attributes,
        #  but include them as information.)
        #  A single observation might affect multiple transitions as a part of
        #  the risk set, but an event belongs to only one.
  	time <- c(y[deaths,2])  # 'c' kills all of the attributes
        dimnames(rr) <- list(time, names(object$coef))

        tran <- object$rmap[,2] # the observed transitions
        attr(rr, "transition") <- 
            factor((tran[ord])[deaths], unique(tran), utran)
        if (length(strats))
            attr(rr, "strata") <- (xstack$strata[ord])[deaths]
	return(rr)
    }   

    if (type=='score') {
        storage.mode(y) <- storage.mode(x) <- "double"
        storage.mode(newstrat) <- "integer"
        storage.mode(score) <- storage.mode(weights) <- "double"
 	if (ny==2) {
	    resid <- .Call(Ccoxscore2, 
                           y, 
                           x, 
                           istrat,
                           score,
                           weights[ord],
                           as.integer(method=='efron'))
        }
	else {
	    resid<- .Call(Cagscore3,
                           y, 
                           x, 
                           istrat,
                           score,
                           weights[ord],
                           as.integer(method=='efron'), 
                           sort1 -1L)
        }
        
	if (nvar >1) {
	    rr <- matrix(0, n, nvar)
	    rr[ord,] <- resid
            dimnames(rr) <- list(names(object$residuals), 
                                     names(object$coefficients))
        }
	else rr[ord] <- resid

	if      (otype=='dfbeta') {
	    if (is.matrix(rr)) rr <- rr %*% vv
	    else               rr <- rr * vv
        }
	else if (otype=='dfbetas') {
	    if (is.matrix(rr))  rr <- (rr %*% vv) %*% diag(sqrt(1/diag(vv)))
	    else                rr <- rr * sqrt(vv)
        }
    }
    
    # The result has rows for transition 1, then rows for transition 2,
    #   etc. It also will have separate columns for sex_1:2 and sex_1:3 for 
    #   instance.  
    # Our target is an array of (observation, coefficient)
    # Weights are the simplest, do them first
    if (weighted) rr <- rr * weights[object$rmap[,1]]

    # Since each row of the current rr goes with exactly 1 transition,
    #  only one of sex_1:2, sex_1:3 etc in a given row will be nonzero.
    #  Because of this we can use the unadorned rowsum function (which is fast)
    if (collapse){
        newr <- rowsum(rr, cluster, reorder=TRUE)
        dimnames(newr) <- list(levels(cluster), names(object$coefficients))
    } else {
        if (TRUE) newr <- rowsum(rr, xstack$rindex, reorder=TRUE)
        else{ # is placing each resid by hand faster than rowsum?
            newr <- matrix(0, nrow(mf), length(object$coefficients))
            obs <- object$rmap[,1] # where each row in rr comes from
            for (i in 1:ncol(cmap)) {
                keep <- object$rmap[,2]==i # this row applies to this trans
                j <- cmap[cmap[,i]>0, i] # coefficient ind for this transition
                newr[obs[keep], j] <- rr[keep, j]
            }
        }
        dimnames(newr) <- list(rownames(mf), names(object$coefficients))
        if (class(omit)== 'exclude')
            newr <- naresid(omit, newr)
    }
    newr
}


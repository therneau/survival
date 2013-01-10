predict.survreg <-
    function(object, newdata, type=c('response', "link", 'lp', 'linear',
				     'terms', 'quantile','uquantile'),
				se.fit=FALSE,  terms=NULL,
	                        p=c(.1, .9), na.action=na.pass, ...) {
    #
    # What do I need to do predictions ?
    #   
    #  linear predictor: exists
    #           +se    : X matrix
    #          newdata : new X matrix
    #
    #  response -- same as lp, +transform, from distribution
    #  
    #  p --  density function from distribution
    #          scale(s) -- if multiple I need the strata
    #          +se : variance matrix
    #	   newdata: new X

    type <-match.arg(type)
    if (type=='link') type<- 'lp'  #true until their are link functions
    if (type=='linear') type<- 'lp'
    n <- length(object$linear.predictors)
    Terms <- object$terms
    if(!inherits(Terms, "terms"))
	    stop("invalid terms component of  object")

    strata <- attr(Terms, 'specials')$strata
    Terms <- delete.response(Terms)
    coef <- object$coefficients
    intercept <- attr(Terms, "intercept")
    nvar <- length(object$coefficients)
    vv <- object$var[1:nvar, 1:nvar]
    fixedscale <- (nvar == ncol(object$var))

    if (missing(newdata) && (type=='terms' || se.fit)) need.x <- TRUE
    else  need.x <- FALSE

    if (!missing(newdata)){
        newframe <- model.frame(Terms, data=newdata, na.action= na.action,
                                xlev=object$xlevels)
        na.action.used <- attr(newframe, 'na.action')
        }
    else na.action.used <- object$na.action

    if (length(strata) && (type=='quantile' || type=='uquantile') &&
	      !fixedscale) {
	#
	# We need to reconstruct the original "strata" variable
	#
	mf <- model.frame(object)
	temp <- untangle.specials(Terms, 'strata', 1)
	dropx <- temp$terms
	if (length(temp$vars)==1) strata.keep <- mf[[temp$vars]]
	else strata.keep <- strata(mf[,temp$vars], shortlabel=TRUE)
	strata <- as.numeric(strata.keep)
	nstrata <- max(strata)
	    
	if (missing(newdata) && need.x){  #need the old x
	    x <- object[['x']] 
	    if (is.null(x)) x <- model.matrix(object, mf)
	    }

	else if (!missing(newdata)) {  #need the new x
	    if (length(temp$vars)==1) newstrat <- newframe[[temp$vars]]
	    else newstrat <- strata(newframe[,temp$vars], shortlabel=TRUE)
	    strata <- match(newstrat, levels(strata.keep))
	    x <- model.matrix(object, newframe)
	    offset <- model.offset(newframe)
	    }
	}

    else {  # per subject strata not needed
        nstrata <- 1
        if (missing(newdata)) {
            strata <- rep(1L, n)
            if (need.x) x <- model.matrix(object)
        }
	
        else {
            x <- model.matrix(object, newframe)
            strata <- rep(1L, nrow(x))
            offset <- 0
        }
    }
    scale <- object$scale[strata]
    #center x if terms are to be computed
    if(type=='p' || (type == "terms" && intercept)) 
	    x <- sweep(x, 2, object$means)

    #
    # Grab the distribution
    #
    if (is.character(object$dist)) dd <- survreg.distributions[[object$dist]]
    else dd <- object$dist
    if (is.null(dd$itrans)) {
	itrans <- function(x) x  # identity transformation
        dtrans <- function (x) 1 # derivative of the transformation
	}
    else {
	itrans <- dd$itrans
	dtrans <- dd$dtrans
	}
    if (!is.null(dd$dist)) dd <- survreg.distributions[[dd$dist]]

    #
    # Now, lay out the code one case at a time.
    #  There is some repetition this way, but otherwise the code just gets
    #    too complicated.
    #
    if (type=='lp' || type=='response') {
	if (missing(newdata)) {
 	    pred <- object$linear.predictors
#	    names(pred) <- names(object$residuals)
	    }
	else  pred <- drop(x %*% coef)  + offset
	if (se.fit) se <- sqrt(diag(x %*% vv %*% t(x)))

	if (type=='response') {
	    pred <- itrans(pred)
	    if (se.fit) se <- se/ dtrans(pred)
	    }
	}
    else if (type=='quantile' || type=='uquantile') {
	if (missing(newdata)) pred <- object$linear.predictors
	else  pred <- x %*% coef 
	# "pred" is the mean of the distribution,
	#   now add quantiles and then invert
	qq <- dd$quantile(p, object$parm)
	if (length(qq)==1 || length(pred)==1) {
	    pred <- pred + qq*scale
	    if (se.fit && fixedscale) {
		var <- ((x %*% vv) * x) %*% rep(1., ncol(x))
		se <- rep(sqrt(drop(var)), length(qq))
		}
	    else if (se.fit) {
		x.strata <- outer(strata, 1:nstrata, 
				  function(x,y) 1*(x==y))
		se <- matrix(0, ncol=length(qq), nrow=nrow(x))
		for (i in 1:(length(qq))) {
		    temp <- cbind(x, (qq[i]*scale)* x.strata)
		    var <- ((temp %*% object$var) *temp) %*% rep(1, ncol(temp))
		    se[,i] <- sqrt(drop(var))
		    }
		se <- drop(se)
		}
	    }
	else {
	    pred <- c(pred) + outer(scale, qq)
	    if (se.fit && fixedscale) {
		var <- ((x %*% vv) * x) %*% rep(1., ncol(x))
		if (length(qq) >1) {
		    se <- rep(sqrt(drop(var)), length(qq))
		    se <- matrix(se, ncol=length(qq))
		    }
		else se <- sqrt(drop(var))
		}
	    else if (se.fit) {
		x.strata <- outer(strata, 1:nstrata, 
				  function(x,y) 1*(x==y))
		se <- pred
		nc <- rep(1., ncol(object$var))
		for (i in 1:length(qq)) {
		    temp <- cbind(x, (qq[i]*scale)*x.strata)
		    var <- ((temp %*% object$var)* temp) %*% nc
		    se[,i] <- sqrt(drop(var))
		    }
		se <- drop(se)
		}
	    }
	pred <- drop(pred)
	if (type == 'quantile') {
	    pred <- itrans(pred)
	    if (se.fit) se <- se/dtrans(pred)
	    }
	}

    else {  #terms
	if (is.R()) {
	    # In S we can use Build.terms, in R we have to do it ourselves
	    asgn <- attrassign(x,Terms)
	    hasintercept<-attr(Terms,"intercept")>0
	    if (hasintercept)
		    asgn$"(Intercept)"<-NULL
	    nterms<-length(asgn)
	    pred<-matrix(ncol=nterms,nrow=NROW(x))
	    dimnames(pred)<-list(rownames(x),names(asgn))
	    if (se.fit){
		se<-matrix(ncol=nterms,nrow=NROW(x))
		dimnames(se)<-list(rownames(x),names(asgn))
		R<-object$var
		ip <- double(NROW(x))
		}
	    for (i in 1:nterms){
		ii<-asgn[[i]]
		pred[,i]<-x[,ii,drop=FALSE]%*%(coef[ii])
		if (se.fit){
		    for(j in (1:NROW(x))){
			xi<-x[j,ii,drop=FALSE]*(coef[ii])
			vci<-R[ii,ii]
			se[j,i]<-sqrt(sum(xi%*% vci %*%t( xi)))
			}
		    }
		}
	    if (!is.null(terms)){
		pred<-pred[,terms,drop=FALSE]
		if (se.fit)
			se<-se[,terms,drop=FALSE]
		}
	    if (hasintercept)
		    const <- coef(object)['(Intercept']
		else    const <- 0
	    }

#  Splus code, commented out to stop a warning from R CMD check
#	else {
#	    # Splus: use Build.terms to do the work
#	    asgn <- attr(x, 'assign')
#	    attr(x, 'constant') <- object$means
#	    terms <- match.arg(Terms, labels.lm(object))
#	    asgn <- asgn[terms]
#
#	    if (se.fit) {
#		temp <- Build.terms(x, coef, vv, asgn, FALSE)
#		pred <- temp$fit
#		se   <- temp$se.fit
#		}
#	    else pred<- Build.terms(x, coef, NULL, asgn, FALSE)
#	    const<- attr(pred, 'constant')
#	    }
        }

    #Expand out the missing values in the result
    # 
    if (!is.null(na.action.used)) {
	pred <- naresid(na.action.used, pred)
	if(se.fit) se <- naresid(na.action.used, se)

	
	}
    if (se.fit) list(fit=pred, se.fit=se)
    else pred
    }

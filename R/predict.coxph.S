#$Id: predict.coxph.S 11370 2009-09-30 18:02:09Z therneau $
#What do I need to do predictions --
#
#linear predictor:  exists
#        +se     :  X matrix
#        +newdata:  means of old X matrix, new X matrix, new offset
#
#risk -- same as lp
#
#expected --    cumulative hazard for subject= baseline haz + time + risk
#        +se :  sqrt(expected)
#      +new  :  baseline hazard function, new time, new x, means of old X,
#                        new offset, new strata
#
#terms -- : X matrix and the means
#    +se  :  ""  + I matrix
#   +new  : new X matrix and the old means + I matrix
predict.coxph <- function(object, newdata, 
		       type=c("lp", "risk", "expected", "terms"),
		       se.fit=FALSE, na.action=na.pass,
		       terms=names(object$assign), collapse, ...)

    {
    type <-match.arg(type)
    n <- object$n
    Terms <-  delete.response(object$terms)
    strata <- attr(Terms, 'specials')$strata
    dropx <- NULL
    if (length(strata)) {
	   temp <- untangle.specials(Terms, 'strata', 1)
	   dropx <- temp$terms
	   }
    if (length(attr(Terms, 'specials')$cluster)) {
	temp <- untangle.specials(Terms, 'cluster', 1)
	dropx <- c(dropx, temp$terms)
	}
    if (length(dropx)) Terms2 <- Terms[-dropx]
    else  Terms2 <- Terms

    offset <- attr(Terms, "offset")
    resp <- attr(Terms, "variables")[attr(Terms, "response")]

    if (missing(newdata)) {
        na.action.used <- object$na.action
	if (type=='terms' || (se.fit && (type=='lp' || type=='risk'))) {
	    x <- object[['x']]  #don't accidentally grab 'xlevels'
	    if (is.null(x)) {
		x <- model.matrix(Terms2, model.frame(object),
                                  contr=object$contrasts)[,-1,drop=FALSE]
		}
	    x <- sweep(x, 2, object$means)
	    }
	else if (type=='expected') {
	    y <- object$y
	    if (is.null(y)) {
		m <- model.frame(object)
		y <- model.extract(m, 'response')
		}
	    }
	}
    else {
	if (type=='expected'){
	     m <- model.frame(Terms, data=newdata, xlev=object$xlevels,
                              na.action=na.action)
             x <- model.matrix(Terms2, m,
                               contr=object$contrasts)[,-1,drop=FALSE]
         }
	else {
            m <- model.frame(Terms2, data=newdata, xlev=object$xlevels,
                             na.action=na.action)
            x <- model.matrix(delete.response(Terms2), m,
                              contr=object$contrasts)[,-1,drop=FALSE]
        }
        na.action.used <- attr(m, 'na.action')

	x <- sweep(x, 2, object$means)
	if (length(offset)) offset<- model.offset(m)
	else offset <- 0
	}

    #
    # Now, lay out the code one case at a time.
    #  There is some repetition this way, but otherwise the code just gets
    #    too complicated.
    if (is.null(object$coefficients))
        coef<-numeric(0)
    else {
	# Replace any NA coefs with 0, to stop NA in the linear predictor
        coef <- ifelse(is.na(object$coefficients), 0, object$coefficients)
        }
    if (type=='lp' || type=='risk') {
	if (missing(newdata)) {
	    pred <- object$linear.predictors
	    names(pred) <- names(object$residuals)
	    }
	else                  pred <- x %*% coef  + offset
	if (se.fit) se <- sqrt(diag(x %*% object$var %*% t(x)))

	if (type=='risk') {
	    pred <- exp(pred)
	    if (se.fit) se <- se * sqrt(pred)
	    }
	}

    else if (type=='expected') {
	if (missing(newdata)) pred <- y[,ncol(y)] - object$residuals
	else  stop("Method not yet finished")
	se   <- sqrt(pred)
	}

    else {  # type = terms, which is different for R and S.  In S, the
	# function Build.terms does all the work --- R doesn't have it.
	if (is.R()) {
	    asgn <- object$assign
	    nterms<-length(terms)
	    pred<-matrix(ncol=nterms,nrow=NROW(x))
	    if (is.character(terms))
		    termnames<-terms
	    else
		    termnames<-names(object$assign)[terms]
	    dimnames(pred)<-list(rownames(x),termnames)
	    if (se.fit){
		se<-matrix(ncol=nterms,nrow=NROW(x))
		dimnames(se)<-list(rownames(x),termnames)
		R<-object$var
		ip <- real(NROW(x))
		}
	    for (i in 1:nterms){
		ii<-asgn[[terms[i] ]]
		pred[,i]<-x[,ii,drop=FALSE]%*%(coef[ii])
		if (se.fit){
		    for(j in (1:NROW(x))){
			xi<-x[j,ii,drop=FALSE]
			vci<-R[ii,ii]
			se[j,i]<-sqrt(sum(xi%*% vci %*%t( xi)))
			}
		    }
		}
	    }
	else {  #terms for Splus
	    attr(x, "constant") <- rep(0, ncol(x))
	    asgn <- object$assign
	    terms <- match.arg(Terms2, labels.lm(object))
	    asgn <- asgn[terms]
	    if (se.fit) {
		temp <- Build.terms(x, coef, object$var, asgn, FALSE)
		pred <- temp$fit
		se   <- temp$se.fit
		}
	    else pred<- Build.terms(x, coef, NULL, asgn, FALSE)
	    }
	}
 
        
    # Terms should always return a matrix, or R termplot() gets unhappy.
    #   So the next two lines are commented out.
    # if (se.fit) se <- drop(se)
    # pred <- drop(pred)

    # Expand out the missing values in the result
    if (!is.null(na.action.used)) {
	pred <- naresid(na.action.used, pred)
        if (is.matrix(pred)) n <- nrow(pred)
	else               n <- length(pred)
	if(se.fit) se <- naresid(na.action.used, se)
	}

    # Collapse over subjects, if requested
    if (!missing(collapse)) {
	if (length(collapse) != n) stop("Collapse vector is the wrong length")
	pred <- rowsum(pred, collapse)  # in R, rowsum is a matrix, always
	if (se.fit) se <- sqrt(rowsum(se^2, collapse))
        if (type != 'terms') {
            pred <- drop(pred)
            if (se.fit) se <- drop(se)
            }
	}

    if (se.fit) list(fit=pred, se.fit=se)
    else pred
    }

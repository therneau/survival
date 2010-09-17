# $Date: 2006-08-28 18:12:46 $ $Id: survfit.formula.S 11430 2010-07-29 22:34:23Z therneau $
survfit.formula <- function(formula, data, weights, subset, 
			    na.action, etype, id, ...) {
    Call <- match.call()
    Call[[1]] <- as.name('survfit')  #make nicer printout for the user
    # create a copy of the call that has only the arguments we want,
    #  and use it to call model.frame()
    mfnames <- c('formula', 'data', 'weights', 'subset','na.action',
                 'etype', 'id')  #legal args for model.frame
    temp <- Call[c(1, match(mfnames, names(Call), nomatch=0))]
    temp[[1]] <- as.name("model.frame")
    if (is.R()) m <- eval.parent(temp)
    else        m <- eval(temp, sys.parent())
    
    Terms <- terms(formula, 'strata')
    ord <- attr(Terms, 'order')
    if (length(ord) & any(ord !=1))
	    stop("Interaction terms are not valid for this function")

    n <- nrow(m)
    Y <- model.extract(m, 'response')
    if (!is.Surv(Y)) stop("Response must be a survival object")

    casewt <- model.extract(m, "weights")
    if (is.null(casewt)) casewt <- rep(1,n)

    if (!is.null(attr(Terms, 'offset'))) warning("Offset term ignored")

    ll <- attr(Terms, 'term.labels')
    if (length(ll) == 0) X <- factor(rep(1,n))  # ~1 on the right
    else X <- strata(m[ll])
    
    etype <- model.extract(m, 'etype')
    id    <- model.extract(m, 'id')

    if (!is.Surv(Y)) stop("y must be a Surv object")

    # Round the survival times.  Otherwise unique() and
    #  table() can give a different list of unique times, and both are
    #  used in the subsequent routines.  A test case with times of
    #  sqrt(2)^2 and 2 revealed this issue.  Normally double.digits is
    #  53 in base 2, which leads to 15 decimal digits.
    digits <- floor((.Machine$double.digits -1) * 
                    logb(.Machine$double.base,10)) #base 10 digits
    Y[,1] <- signif(Y[,1], digits)
    if (ncol(Y)==3) Y[,2] <- signif(Y[,2],digits)

    if (length(etype) >0)
         temp <- survfitCI(X, Y, weights=casewt, etype=etype, id=id,  ...)
    else if (attr(Y, 'type') != 'right' && attr(Y, 'type') != 'counting')
         temp <- survfitTurnbull(X, Y, casewt, ...)
    else temp <- survfitKM(X, Y, casewt, ...)

    if (is.R()) class(temp) <- 'survfit'
    else        oldClass(temp) <- "survfit"
    if (!is.null(attr(m, 'na.action')))
	    temp$na.action <- attr(m, 'na.action')

    temp$call <- Call
    temp
    }










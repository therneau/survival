# $Id: coxph.rvar.S 11166 2008-11-24 22:10:34Z therneau $
coxph.rvar <- function(fit, collapse) {
    rcall <- match.call()
    if (class(fit) != 'coxph')
	stop ("First argument must be a fitted Cox model")

    if (missing(collapse)) temp <- residuals.coxph(fit, type='dfbeta')
    else temp <- residuals.coxph(fit, type='dfbeta', collapse=collapse)
    if (any(is.na(temp)))
       if (ncol(temp)==1) temp<- temp[!is.na(temp),,drop=FALSE]
       else               temp <- temp[!is.na(temp %*% rep(1,ncol(temp))),]
    fit$robust.var <- t(temp) %*% temp
    fit$rcall <- rcall
    fit
    }

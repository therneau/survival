points.survfit <- function(x, fun, censor=FALSE,
                           col=1, pch, noplot="(s0)", cumhaz=FALSE, ...) {
    x <- survfit0(x, x$start.time)

    # Organize data into stime, ssurv, supper, slower
    stime <- x$time
    std   <- NULL
    yzero <- FALSE   # a marker that we have an "ordinary survival curve" with min 0
    smat <- function(x) {
        # the rest of the routine is simpler if everything is a matrix
        dd <- dim(x)
        if (is.null(dd)) as.matrix(x)
        else if (length(dd) ==2) x
        else matrix(x, nrow=dd[1])
    }

    if (is.numeric(cumhaz)) { # plot the cumulative hazard
        if (!inherits(x, "survfitms") && any(cumhaz != 1))
            stop("numeric cumhaz argument only applies to multi-state")
        dd <- dim(x$cumhaz)
        if (is.null(dd)) nhazard <- 1
        else nhazard <- prod(dd[-1])

        if (!all(cumhaz == floor(cumhaz))) stop("cumhaz argument is not integer")
        if (any(cumhaz < 1 | cumhaz > nhazard)) stop("subscript out of range")
        ssurv <- smat(x$cumhaz)[,cumhaz, drop=FALSE]
        if (!is.null(x$std.chaz)) std <- smat(x$std.chaz)[,cumhaz, drop=FALSE]
        cumhaz <- TRUE # for the rest of the code
    } else if (cumhaz) {
        if (is.null(x$cumhaz)) 
            stop("survfit object does not contain a cumulative hazard")
        ssurv <- smat(x$cumhaz)
        if (!is.null(x$std.chaz)) std <- smat(x$std.chaz)
    }
    else if (inherits(x, "survfitms")) {
        if (!missing(cumprob) && !(is.logical(cumprob) && !cumprob)) {
            dd <- dim(x)
            j <- match("states", names(dd), nomatch=0)
            if (j==0) stop("survfitms object with no states dimension")
            
            #  cumprob is T/F or a vector of integers
            if (is.logical(cumprob)) cumprob <- 1:dd[j]
            else if (!is.numeric(cumprob) || any(cumprob <1 | cumprob > dd[j])
                     || any(cumprob != floor(cumprob)))
                stop("cumprob contains an invalid numeric")
                
            if (dd[j] ==1) {
                # nothing to do, user subscripted to only 1 state
                ssurv <- x$pstate
            } else {
                # reorder the states, pstate has dimension 2 or 3,
                #  time/strata is first, data (if present), then states
                #  (dd is the dimension from the user's point of view, of
                #  strata, data, state)
                if (length(dim(x$pstate))==2) {
                    # drop = FALSE for the rare case of a single time point
                    ssurv <- t(apply(x$pstate[,cumprob, drop=FALSE],1,cumsum))
                } else {
                    temp <- apply(x$pstate[,,cumprob, drop=FALSE],1:2, cumsum)
                    ssurv <- smat(aperm(temp, c(2,3,1)))
                }
                cumprob <- TRUE  # for the lastx line
            }
        } else {
            i <- !(x$states %in% noplot)
            if (all(i) || !any(i)) {
                # the !any is a failsafe, in case none are kept we ignore noplot
                ssurv <- smat(x$pstate)
                if (!is.null(x$std.err)) std <- smat(x$std.err)
                if (!is.null(x$lower)) {
                    slower <- smat(x$lower)
                    supper <- smat(x$upper)
                }
            }
            else {
                i <- which(i)  # the states to keep
                # we have to be careful about subscripting
                if (length(dim(x$pstate)) ==3) {
                    ssurv <- smat(x$pstate[,,i, drop=FALSE])
                    if (!is.null(x$std.err))
                        std <- smat(x$std.err[,,i, drop=FALSE])
                    if (!is.null(x$lower)) {
                        slower <- smat(x$lower[,,i, drop=FALSE])
                        supper <- smat(x$upper[,,i, drop=FALSE])
                    }
                }
                else {
                    ssurv <- x$pstate[,i, drop=FALSE]
                    if (!is.null(x$std.err)) std <- x$std.err[,i, drop=FALSE]
                    if (!is.null(x$lower)) {
                        slower <- smat(x$lower[,i, drop=FALSE])
                        supper <- smat(x$upper[,i, drop=FALSE])
                    }
                }
            }
        }
    }
    else {
        yzero <- TRUE
        ssurv <- as.matrix(x$surv)   # x$surv will have one column
        if (!is.null(x$std.err)) std <- as.matrix(x$std.err)
        # The fun argument usually applies to single state survfit objects
        #  First deal with the special case of fun='cumhaz', which is here for
        #  backwards compatability; people should use the cumhaz argument
        if (!missing(fun) && is.character(fun) && fun=="cumhaz") {
            cumhaz <- TRUE
            if (!is.null(x$cumhaz)) {
                ssurv <- as.matrix(x$cumhaz)
                if (!is.null(x$std.chaz)) std <- as.matrix(x$std.chaz)
            } 
            else {
                ssurv <- as.matrix(-log(x$surv))
                if (!is.null(x$std.err)) {
                    if (x$logse) std <- as.matrix(x$std.err)
                    else std <- as.matrix(x$std.err/x$surv)
                }
             }
        }
    }

    # set up strata
    if (is.null(x$strata)) {
        nstrat <- 1
        stemp <- rep(1, length(x$time)) # same length as stime
    }
    else {
        nstrat <- length(x$strata)
        stemp <- rep(1:nstrat, x$strata) # same length as stime
    }
    ncurve <- nstrat * ncol(ssurv)
    if (!missing(fun)){
        if (is.character(fun)) {
            if (cumhaz) {
                tfun <- switch(tolower(fun),
                               'log' = function(x) x,
                               'cumhaz'=function(x) x,
                               'identity'= function(x) x,
                               stop("Invalid function argument")
                               )
            } else if (inherits(x, "survfitms")) {
                tfun <-switch(tolower(fun),
                              'log' = function(x) log(x),
                              'event'=function(x) x,
                              'cloglog'=function(x) log(-log(1-x)),
                              'cumhaz' = function(x) x,
                              'pct' = function(x) x*100,
                              'identity'= function(x) x,
                              stop("Invalid function argument")
                              )
            } else {
                yzero <- FALSE
                tfun <- switch(tolower(fun),
                           'log' = function(x) x,
                           'event'=function(x) 1-x,
                           'cumhaz'=function(x) x,
                           'cloglog'=function(x) log(-log(x)),
                           'pct' = function(x) x*100,
                           'logpct'= function(x) 100*x,  #special case further below
                           'identity'= function(x) x,
                           'f' = function(x) 1-x,
                           's' = function(x) x,
                           'surv' = function(x) x,
                           stop("Unrecognized function argument")
                           )
            }
        }
        else if (is.function(fun)) tfun <- fun
        else stop("Invalid 'fun' argument")
        
        ssurv <- tfun(ssurv )
        if (!is.null(supper)) {
            supper <- tfun(supper)
            slower <- tfun(slower)
        }
    }
    
    if (ncurve==1 || (length(col)==1 && missing(pch))) {
        if (censor) points(stime, ssurv, ...)
        else points(stime[x$n.event>0], ssurv[x$n.event>0], ...)
    }
    else {
        c2 <- 1  #cycles through the colors and characters
        col <- rep(col, length=ncurve)
        if (!missing(pch)) {
            if (length(pch)==1)
                pch2 <- rep(strsplit(pch, '')[[1]], length=ncurve)
            else pch2 <- rep(pch, length=ncurve)
        }
        for (j in 1:ncol(ssurv)) {
            for (i in unique(stemp)) {
                if (censor) who <- which(stemp==i)
                else who <- which(stemp==i & x$n.event >0)
                if (missing(pch))
                    points(stime[who], ssurv[who,j], col=col[c2], ...)
                else
                    points(stime[who], ssurv[who,j], col=col[c2], 
                           pch=pch2[c2], ...) 
                c2 <- c2+1
            }
        }
    }
}

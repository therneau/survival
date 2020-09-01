#
# The redistribute-to-the-right algorithm
#  Useful for teaching, and for data checks
#  Input
#    formula: usually Surv(time, status) ~1.  Any variables on the
#      right hand side are not used
#    data, weights, subset, na.action: usual formula arguments
#    times: vector of reporting times.  If omitted, the largest event time
#      is used
#   
#  Output: a vector (one output time) or matrix of weights, one per observation
#
rttright <- function(formula, data, weights, subset, na.action, time,
                     id, timefix=TRUE) {
    Call <- match.call()  # save a copy of the call

    # The first chunk of this function is essentially a copy of survfit.
    # This is more fussy than needed, but the extra complexity save us from
    #   rare edge cases.  User's do odd things. (I expect
    #   someone to just paste in one of their coxph calls). Because
    #   we have it worked out, we use it.  
    #
    indx <- match(c('formula', 'data', 'weights', 'subset','na.action',
                    'istate', 'id', 'cluster', "etype"), names(Call), nomatch=0)
    #It's very hard to get the next error message other than malice
    #  eg survfit(wt=Surv(time, status) ~1) 
    if (indx[1]==0) stop("a formula argument is required")
    temp <- Call[c(1, indx)]
    temp[[1L]] <- quote(stats::model.frame)
    mf <- eval.parent(temp)

    Terms <- terms(formula, c("strata", "cluster"))
    ord <- attr(Terms, 'order')
    if (length(ord) & any(ord !=1))
            stop("Interaction terms are not valid for this function")

    n <- nrow(mf)
    Y <- model.response(mf)
    ny <- ncol(Y)
    if (!is.Surv(Y)) stop("response must be a Surv object")
    if (!attr(Y, 'type') == "right" && !attr(Y, 'type')== "counting")
         stop("invalid survival type")
    
    casewt <- model.extract(mf, "weights")
    if (is.null(casewt)) casewt <- rep(1, n)
    else {
        if (!is.numeric(casewt)) stop("weights must be numeric")
        if (any(!is.finite(casewt))) stop("weights must be finite") 
        if (any(casewt <0)) stop("weights must be non-negative")
        casewt <- as.numeric(casewt)  # transform integer to numeric
    }

    if (!is.null(attr(Terms, 'offset'))) warning("Offset term ignored")

    temp <- untangle.specials(Terms, "cluster")
    if (length(temp$vars)>0) {
        Terms <- Terms[-temp$terms]
    }

    ll <- attr(Terms, 'term.labels')
    if (length(ll) == 0) X <- factor(rep(1,n))  # ~1 on the right
    else X <- strata(mf[ll])

    # Deal with the near-ties problem
    if (!is.logical(timefix) || length(timefix) > 1)
        stop("invalid value for timefix option")
    if (timefix) Y <- aeqSurv(Y) 
        
    id <- model.extract(mf, "id")

    # Finally, it's time to do the actual work.  
    # The compact algorithm depends on G(t), the censoring distribution
    # Two issues:  1. G(t) is left continuous, the KM is right continuous
    #              2. For tied times, censors happen after deaths
    # For an id with multiple observations, only the last can be a 
    #   true censoring
    if (is.null(id)) censor <- (Y[,ny] ==0)
    else {
        ord <- order(id)
        ctemp <- (!duplicated(id[ord]) & Y[ord, ny]==0)
        censor <- vector("logical", n)
        censor[ord] <- ctemp
    }
    ctime <- unique(Y[censor,ny-1])  # the unique censoring times
    event <- (Y[,ny] > 0)
    etime <- unique(Y[event, ny-1])  # unique event times
    if (any(duplicated(c(etime, ctime)))) {
        eps <- min(diff(sort(unique(c(ctime, etime)))))/2  # an offset
        y2 <- Y
        y2[,ny] <- ifelse(censor, 1, 0)
        y2[censor, ny-1] <- y2[censor, ny-1] + eps  # break the ties
    }  
    else {
        y2 <- Y
        y2[,ny] <- ifelse(censor, 1, 0)
        eps <- 0
    }
    G <- survfitKM(X, y2, casewt, se.fit=FALSE)
    # read off separately for each stratum, and for each time
    istrat <- as.numeric(X)
    nstrat <- max(istrat)
    if (nstrat>1) gstrat <- rep(1:nstrat, G$strata)
    else gstrat <- rep(1, length(G$time))

    gread <- function(y, weight, gtime, gsurv, cut) {
        if (length(gtime) ==0) new <- weight # no censorings
        else {
            gindx <- findInterval(pmin(cut, y[,ny-1]), gtime, left.open=TRUE)
            new <- ifelse(y[,ny]==0 & y[,ny-1] < cut, 0,
                   weight/c(1,gsurv)[1+gindx])
        }
        unname(new)
    } 
            
    if (missing(time)) time <- 2*max(abs(G$time)) +1  # past the last
    wtmat <- matrix(casewt, n, length(time))
    for (i in 1:nstrat) {
        keep <- (gstrat==i & G$n.event > 0)
        if (any(keep)) {
            gtime <- G$time[keep] - eps
            for (j in 1:length(time)) {
                wtmat[istrat==i, j] <- gread(Y[istrat==i,], casewt[istrat==i],
                                             gtime, G$surv[keep], time[j])
            }       
        }
    }   
    colnames(wtmat) <- time
    drop(wtmat)  
}


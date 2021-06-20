#
# The redistribute-to-the-right algorithm
#  Useful for teaching, and for data checks
#  Input
#    formula: usually Surv(time, status) ~1.  
#    data, weights, subset, na.action: usual formula arguments
#    times: vector of reporting times.  If omitted, the largest event time
#      is used
#   
#  Output: a vector (one output time) or matrix of weights, one per observation
#
rttright <- function(formula, data, weights, subset, na.action, times,
                     id, timefix=TRUE) {
    Call <- match.call()  # save a copy of the call

    # The first chunk of this function is essentially a copy of the code
    #  from the survfit function, the exception being how covariates are handled
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
    type <- attr(Y, "type")
    if (!(type %in% c("right", "mright", "counting", "mcounting")))
        stop("response must be right censored")
    
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
    if (ny==3) {
        if (is.null(id)) stop("id is required for start-stop data")
    }
    if (!is.null(id)) {
        
        if (is.null(attr(Y, 'states'))) {
            ytemp <- Y
            attr(ytemp, 'states') <- 'fail'  # survcheck2 wants a states attr
            check <- survcheck2(ytemp, id)
        }
        else check <- survcheck2(Y, id)

        if (any(check$flag > 0)) 
                stop("one or more flags are >0 in survcheck")
        n.startstate <- sum(check$transitions[,1] >1)
        if (ny ==2) samestart=TRUE
        else {
           etemp  <- tapply(Y[,1], id, min)
           samestart <- all(temp==temp[1])
        }    
    } else check <- NULL

    # Finally, it's time to do the actual work. 
    # For simple survival or competing risks data,
    #  we can use the a compact algorithm that needs only one call to survfit
    # The compact algorithm depends on G(t), the censoring distribution
    if (is.null(check) || (n.startstate==1 & samestart)) {
        # Compute the censoring distribution $G$
        # 1. G(t) is left continuous, the KM is right continuous
        # 2. For tied times, censors happen after deaths
        # 3. For an id with multiple observations, censorings in the middle
        #     are not true censoring
        #
        # Mark the actual censoring times, then create a y2 with the censoring
        #  times shifted just a bit.
        if (ny==2) censor <- ifelse(Y[,2]== 0, 1, 0)
        else {
            ord <- order(id, Y[,2])  # time within subject
            ltemp <- !duplicated(id[ord], fromLast=TRUE)  # last for each id
            last <- ltemp[order(ord)]  # marks the last obs, in data order
            censor <- ifelse(last & Y[,2]==0, 1, 0)
        }

        delta <- min(diff(sort(unique(c(Y[,-ny]))))) /2
        y2 <- Y
        y2[,ny] <- censor
        y2[censor==1, ny-1] <- y2[censor==1, ny-1] + delta
        G <- survfitKM(X, y2, casewt, se.fit=FALSE)

        # read off separately for each stratum, and for each time
        istrat <- as.numeric(X)
        nstrat <- max(istrat)
        if (nstrat>1) gstrat <- rep(1:nstrat, G$strata)
        else gstrat <- rep(1, length(G$time))

        # Grab the weights at a given time
        gread <- function(y, weight, gtime, gsurv, cut) {
            if (length(gtime) ==0) new <- weight # no censorings
            else {
                gindx <- findInterval(pmin(cut, y[,ny-1]), gtime, left.open=TRUE)
                new <- ifelse(y[,ny]==0 & y[,ny-1] < cut, 0,
                              weight/c(1,gsurv)[1+gindx])
            }
            unname(new)
        } 
            
        if (missing(times)) times <- 2*max(abs(G$time)) +1  # past the last
        wtmat <- matrix(casewt, n, length(times))
        for (i in 1:nstrat) {
            keep <- (gstrat==i & G$n.event > 0)
            if (any(keep)) {
                gtime <- G$time[keep]
                ikeep <- (istrat==i)  # longer than keep if data has tied times
                for (j in 1:length(times)) {
                    wtmat[ikeep, j] <- gread(Y[istrat==i,], wtmat[ikeep, j],
                                             gtime, G$surv[keep], times[j])
                }
            }
        }
    }       
    else { 
        # The more difficult case, where there are multiple states or delayed
        #   entry
        stop("Code not yet complete for multistate")
    }
    colnames(wtmat) <- times
    drop(wtmat)  
}


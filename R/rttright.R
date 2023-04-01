#
# The redistribute-to-the-right algorithm
#  Useful for teaching, and for data checks
#  Input
#    formula: usually Surv(time, status) ~1.  
#    data, weights, subset, na.action: usual formula arguments
#    times: vector of reporting times.  If omitted, the largest event time
#      is used
#   
#  Output: a vector or matrix of weights, one row per observation
#
rttright <- function(formula, data, weights, subset, na.action, times,
                     id, timefix=TRUE, renorm=TRUE) {
    Call <- match.call()  # save a copy of the call
    if (missing(times)) times <- NULL

    # The first chunk of this function is essentially a copy of the code
    #  from the survfit function, the exception being how covariates are handled
    # This is more fussy than needed, but the extra complexity save us from
    #   rare edge cases.  User's do odd things. (I expect
    #   someone to just paste in one of their coxph calls). Because
    #   we have it worked out, we use it.  
    #
    indx <- match(c('formula', 'data', 'weights', 'subset','na.action',
                    'istate', 'id'), names(Call), nomatch=0)
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
    if (ny==3 && (is.null(id))) stop("id is required for start-stop data")

    if (!is.null(id)) {
        if (is.null(attr(Y, 'states'))) {
            ytemp <- Y
            attr(ytemp, 'states') <- 'event'  # survcheck2 wants a states attr
            check <- survcheck2(ytemp, id)
        }
        else check <- survcheck2(Y, id)

        if (any(check$flag > 0)) 
                stop("one or more flags are >0 in survcheck")
        n.startstate <- sum(check$transitions[,1] >1)
        if (ny ==2) samestart=TRUE   # everyone starts at the same time
        else {
           etemp  <- tapply(Y[,1], id, min)
           samestart <- all(etemp== etemp[1])
        }  

        temp <- tapply(casewt, id, function(x) diff(range(x)))
        if (any(temp > 0)) stop("there are subjects with multiple weights")
    } else check <- NULL

    # Finally, it's time to do the actual work. 
    # For simple survival or competing risks data,
    #  we can use the a compact algorithm that needs only one call to survfit
    # The compact algorithm depends on G(t), the censoring distribution
    if (is.null(check) || (n.startstate==1 & samestart)) {
        istrat <- as.numeric(X)   # X will be a factor
        nstrat <- max(istrat)
        if (renorm) {
            for (i in unique(istrat)) {
                j <- which(istrat==i)
                if (is.null(id)) uj <- j
                else uj <- j[!duplicated(id[j])]  # don't count someone twice
                casewt[j] <- casewt[j]/sum(casewt[uj])
            }
        }

        # Compute the censoring distribution $G$
        # 1. G(t) is left continuous, the KM is right continuous
        # 2. For tied times, censors happen after deaths
        # 3. For an id with multiple observations, censorings in the middle
        #     are not true censoring
        #
        # Mark the actual censoring times.
        last <- NULL
        if (ny==2) censor <- ifelse(Y[,2]== 0, 1, 0)
        else {
            ord <- order(id, Y[,2])  # time within subject
            ltemp <- !duplicated(id[ord], fromLast=TRUE)  # last for each id
            last <- ltemp[order(ord)]  # marks the last obs, in data order
            censor <- ifelse(last & Y[,3]==0, 1, 0)
        }

        # Compute the censoring distribution by creating y2, which is Y with
        #  the censors shifted right by just a bit.  If the type of Y is
        #  mright or mcounting, change it to right or counting
        delta <- min(diff(sort(unique(c(Y[,-ny], times))))) /2
        y2 <- Y
        y2[,ny] <- censor
        y2[censor==1, ny-1] <- y2[censor==1, ny-1] + delta
        attr(y2, "type") <- if (ny==2) "right" else "counting"
        G <- survfitKM(X, y2, casewt, se.fit=FALSE)
        class(G) <- "survfit"

        # The ouput matrix has one column for each requested time point
        if (is.null(times)) { # single time point
            grabwt <- function(Y, GG, casewt, last, times) {
                # See the longer discusson below
                index <- findInterval(Y[,ny-1], GG$time, left.open=TRUE)
                gwt <- c(1,GG$surv)[1+ index]  # gwt at event time
                if (ncol(Y)==2) haswt <- (Y[,2] > 0)
                else haswt <- (last & Y[,3]>0)
                unname(ifelse (haswt, casewt/gwt, 0))
            }
            times <- 0 # a dummy value
        }
        else {
            # for (time1, time2) data the weight is passed forward from row to
            #  row for a person.  The last row, in not censored, retains the
            #  weight.
            # Result is a matrix with one row for each Y, one col for each
            #  time.   Column j of this with be for times[j]
            # If Y[,1] >= times[j] the weight is 0.  The weight for a set of
            #  rows belonging to a given id gets handed off like a baton in a 
            #  relay race, and this row hasn't yet received it.
            # If Y[,2] < times[j] the weight is 0 as well, unless this is the
            #  last, uncensored row for a subject.  That weight preserves: there
            #  is no hand off.  Censors "hand off" by redistribution. 
            # Otherwise, the weight = casewt/gwt[col]: not a last row for an id)
            #                         casewt/max(gwt[col], gwt2[row]): last row
            #
            #  For (time, status) data this is simpler: every row is a last row
            grabwt <- function(Y, GG, casewt, last, times) {
                # G(t) at censoring times
                indx <- findInterval(times, GG$time, left.open=TRUE)
                gwt <- c(1, GG$surv)[1+indx]
                indx2 <- findInterval(Y[,ny-1], GG$time, left.open=TRUE)
                gwt2 <- c(1, GG$surv)[1+indx2]
                
                ntime <- length(times)
                n <- nrow(Y)

                 if (ncol(Y)==2) { 
                    wtmat <- ifelse(outer(Y[,1], times, ">="), 
                                    outer(casewt, gwt, "/"), 0)
                    indx3 <- (Y[,2] > 0)  # events
                    wtmat[indx3,] <- casewt[indx3]/outer(gwt2[indx3], gwt, pmax)
                }
                else {
                    inner <- outer(Y[,1], times, "<") & 
                        outer(Y[,2], times, ">=")
                    wtmat <- ifelse(inner, outer(casewt, gwt, "/"), 0)
                    indx3 <- (Y[,2]> 0 & last)
                    wtmat[indx3,] <-
                        casewt[indx3] / outer(gwt2[indx3], gwt, pmax)
                }
                wtmat
            }
        }   

       if (nstrat ==1) { # only one stratum
           wtmat <- grabwt(Y, G, casewt, last, times)
       } else {
           wtmat <- matrix(0, n, length(times))
           for (i in 1:nstrat) {
               who <- (istrat==i)
               wtmat[who,] <- grabwt(Y[who,], G[i], casewt[who], last[who],
                                     times)
           }  
       }
    } else { 
        # The more difficult case, where there are multiple states or delayed
        #   entry
        stop("function not defined for delayed entry or multistate data")
    }   

    if (length(times)> 1) {
        dimnames(wtmat) <- list(NULL, times)
        wtmat
    }
    else drop(wtmat)  
}

#
# Here is the "cheating" version of the function
#  I know that the weights, done correctly, will reproduce the Aalen-Johansen
#  So compute the AJ and work backwards
rttright2 <- function(formula, data, weights, subset, na.action, times,
                     id, timefix=TRUE, renorm=TRUE) {
    Call <- match.call()  # save a copy of the call

    indx <- match(c('formula', 'data', 'weights', 'subset','na.action',
                    'istate', 'id'), names(Call), nomatch=0)

    if (indx[1]==0) stop("a formula argument is required")
    temp <- Call[c(1, indx)]
    temp[[1L]] <- quote(survival::survfit)
    temp$model <- TRUE
    sfit <- eval.parent(temp) # a KM or Aalen-Johansen

    X <- model.matrix(sfit$model)  # get the grouping variable for each obs
    Y <- model.response(sfit$model)
    n <- nrow(Y)
    if (!renorm) {
        wt <-model.weights(sfit$model)
        if (is.null(wt)) wt <- rep(1, n)
    }

    dd <- dim(sfit)
}    
    
